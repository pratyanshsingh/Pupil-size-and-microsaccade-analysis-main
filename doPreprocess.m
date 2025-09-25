function doPreprocess(projectDir)
%%
rawDir = dir([projectDir '\derivatives\raw\']);
rawDir = rawDir(3:end);
precDir = [projectDir '\derivatives\prec\'];
numSubj = length(rawDir);

for nSubj = 1:numSubj
    %%
    fprintf('\rProgress: %d / %d (%.1f%%)', nSubj, numSubj, nSubj/numSubj*100);
    subIdx = sprintf('\\sub-%02d',nSubj);
    subjDir = [precDir subIdx '\\'];
    if ~exist(subjDir, "dir")
        mkdir(subjDir)
    end

    subRawDir = dir([rawDir(nSubj).folder subIdx '\*.mat']);
    for i = 1:length(subRawDir)
        load([subRawDir(i).folder '\\' subRawDir(i).name])
    end
    subjFileName = sprintf('sub-%02d_task-%s_',nSubj,subInfo.task);

    % Prepare struct
    nData  = struct();
    rmFlag = struct();
    nData.rawData = pupData.trials;

    % Upsample to 1000 Hz
    nData.upSampledData = upSample(nData.rawData);

    % Convert into diameter
    nData.cvtData = cvtDia(nData.upSampledData);

    % Remove blink and invalid interval (Kret & Sjak-Shie, 2018)
    [nData.rbkData, rmFlag.pup_l] = rmBlink(nData.cvtData,'interp','left');
    [nData.rbkData, rmFlag.pup_r] = rmBlink(nData.rbkData,'interp','right');

    % Epoch pupil data
    evtList = pupData.targeton;
    tw      = [-1000 4000];
    subInfo.epoch = tw;
    [l_epData, r_epData, rmFlag.ep] = doEpoch(nData.rbkData,evtList,tw);

    % Baseline correction
    tw = [-200 0];
    tw = tw - subInfo.epoch(1);
    l_epbcData = l_epData - mean(l_epData(tw(1):tw(2),:));
    r_epbcData = r_epData - mean(r_epData(tw(1):tw(2),:));

    % To ddetect sac, remove blink and invalid interval without interpolating
    [nData.nanData, ~] = rmBlink(nData.cvtData,'nan','left');
    [nData.nanData, ~] = rmBlink(nData.nanData,'nan','right');

    % Detect binocular candidates for saccades (Engbert & Mergenthaler,
    % 2006). The position data were smoothed using a 4th-order low-pass
    % Butterworth filter with a cut-off frequency of 200 Hz.
    [sacData, rmFlag.sac] = detectSac(nData.nanData,evtList);

    % Remove saccades those exceed specific criteria
    tw   = [-1000 3999]; % Use whole window, we can epoch small window later
    thr  = [0.1 2]; % Angle amplitude criteria
    vthr = [0 200]; % Velocity criteria
    dur  = false; % Duration of microsaccades
    intv = false;
    Sac = rmSac(sacData, tw, thr, vthr, dur, intv);

    % Extract feature of microsaccads and epoch
    tw   = [-1000 4000];
    [RateMatx, AmpMatx, VelMatx, AngMatx] = getSacMTX(Sac,tw);
%%
    save([subjDir subjFileName 'beh.mat'],        "behData")
    save([subjDir subjFileName 'info.mat'],       "subInfo")
    save([subjDir subjFileName 'l_epData.mat'],   "l_epData")
    save([subjDir subjFileName 'r_epData.mat'],   "r_epData")
    save([subjDir subjFileName 'l_epbcData.mat'], "l_epbcData")
    save([subjDir subjFileName 'r_epbcData.mat'], "r_epbcData")
    save([subjDir subjFileName 'Sac.mat'],        "Sac")
    save([subjDir subjFileName 'RateMatx.mat'],   "RateMatx")
    save([subjDir subjFileName 'AmpMatx.mat'],    "AmpMatx")
    save([subjDir subjFileName 'VelMatx.mat'],    "VelMatx")
    save([subjDir subjFileName 'AngMatx.mat'],    "AngMatx")
    save([subjDir subjFileName 'rmFlag.mat'],     "rmFlag")
end

VarName = [behData.Properties.VariableNames, ...
    'Pupil_l', 'Pupil_r', 'tonicPupil_l', 'tonicPupil_r', ...
    'Rate', 'Amp', 'Vel'];
save([precDir 'VarName.mat'],     "VarName")