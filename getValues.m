function op = getValues(projectDir,dv,tw,iv,ivVal)
%%
precDir = dir([projectDir '\derivatives\prec\sub*']);
numSubj = length(precDir);
op_sub = cell(numSubj, 1); % Preallocate cell array for efficiency

for nSubj = 1:numSubj
    subIdx = sprintf('\\sub-%02d',nSubj);
    subjDir = [precDir(nSubj).folder subIdx '\\'];

    subFile = dir([subjDir '\*info.mat']);
    load([subFile.folder '\\' subFile.name])

    if subInfo.olFlag
        continue
    end
    op = getValues_sj(subjDir,dv,tw,iv,ivVal);

    % if isnan(mean(op.(dv),'all'))
    %     continue
    % end
    sub = table(repmat(nSubj, height(op), 1), 'VariableNames', "Subj");
    op_sub{nSubj} = [sub, op];
end

op = vertcat(op_sub{:});
