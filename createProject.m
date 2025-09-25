function createProject(projectDir)

if nargin < 1
    projectDir = uigetdir(pwd, 'Please choose the project folder');
end
projectSubDir = dir(projectDir);

if ismember('rawdata', {projectSubDir.name})
    rawdataDir = [projectDir '\rawdata\'];
    Status = readtable([rawdataDir, 'SubjectStatus.csv']);
else
    error('There is no rawdata.')
end

taskName = inputdlg('Task name','Input',[1 40],"Stroop");

drvDir = [projectDir '\derivatives\'];
if ~exist(drvDir, 'dir')
    mkdir(drvDir);
end

rawDir = [drvDir '\raw\'];
if ~exist(rawDir, 'dir')
    mkdir(rawDir);
end

pupFiles = dir(fullfile([rawdataDir 'Pupil'], '*.mat'));
behFiles = dir(fullfile([rawdataDir 'Behaviour'], '*.csv'));

numSubj = length(pupFiles);
for nSubj = 1:numSubj
    subjDir = [rawDir sprintf('sub-%02d',nSubj) '\'];
    mkdir(subjDir)

    subjFileName = sprintf('sub-%02d_task-%s',nSubj,taskName{1});

    pupData = load([pupFiles(nSubj).folder '\' pupFiles(nSubj).name]).F;
    save([subjDir '\\' subjFileName '_pup.mat'],"pupData")

    behData = readtable([behFiles(nSubj).folder '\' behFiles(nSubj).name]);
    save([subjDir '\\' subjFileName '_beh.mat'],"behData")

    subjInfo = struct();
    subjInfo.subjName = subjFileName;
    subjInfo.task = taskName{1};
    subjInfo.sampleRate = pupData.SAMPRATE;
    subjId = pupFiles(nSubj).name(4:end-4);
    subjStat = Status(strcmp(Status.used_name,subjId),:);
    subInfo = catstruct(subjInfo,table2struct(subjStat));
    save([subjDir subjFileName '_info.mat'],"subInfo")

    fprintf('\rProgress: %d / %d (%.1f%%)', nSubj, numSubj, nSubj/numSubj*100);
end
end

function out = catstruct(varargin)
    out = struct();
    for k = 1:nargin
        f = fieldnames(varargin{k});
        for i = 1:length(f)
            out.(f{i}) = varargin{k}.(f{i});
        end
    end
end