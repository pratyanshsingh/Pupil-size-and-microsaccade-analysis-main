function excludeOutliers(projectDir)

precDir = dir([projectDir '\derivatives\prec\sub*']);
numSubj = length(precDir);

olFlag = true(numSubj,1);
for nSubj = 1:numSubj
    %%

    subIdx = sprintf('\\sub-%02d',nSubj);
    subjDir = [precDir(nSubj).folder subIdx '\\'];

    subFile = dir([subjDir '\*beh.mat']);
    load([subFile.folder '\\' subFile.name])

    subFile = dir([subjDir '\*info.mat']);
    load([subFile.folder '\\' subFile.name])

    subFile = dir([subjDir '\*rmFlag.mat']);
    load([subFile.folder '\\' subFile.name])

    subFile = dir([subjDir '\*RateMatx.mat']);
    load([subFile.folder '\\' subFile.name])

    ACC = logical(behData.ACC);
    rmFlag.ACC = ~ACC;

    rmFlag.RT = false(length(ACC),1);
    rmFlag.RT(ACC) = isoutlier(behData.RT(ACC),'mean');

    rmTable = [rmFlag.pup_l,rmFlag.ep,rmFlag.ACC,rmFlag.RT,rmFlag.sac];
    rmFlag.overall = any(rmTable')';

    subACC = [mean(behData.ACC(behData.cond==1)), ...
        mean(behData.ACC(behData.cond==2))];

    subInfo.olFlag = false;

    if any(subACC < 0.48)
        subInfo.olFlag = true;
    end

    if sum(RateMatx,"all") == 0
        subInfo.olFlag = true;
    end

    if subInfo.olFlag
        olFlag(nSubj) = false;
    end

    subjFileName = sprintf('sub-%02d_task-%s_',nSubj,subInfo.task);
    save([subjDir subjFileName 'info.mat'],       "subInfo")
    save([subjDir subjFileName 'rmFlag.mat'],     "rmFlag")
end
save([projectDir '\derivatives\prec\olFlag.mat'],       "olFlag")
