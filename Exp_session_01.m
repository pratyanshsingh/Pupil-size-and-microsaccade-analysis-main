%%
addpath(genpath('D:\git-repos\Pupil-size-and-microsaccade-analysis'))

clearvars; clc
projectDir = 'D:\EmoStroop\\';
createProject(projectDir);

clearvars; clc
projectDir = 'D:\EmoStroop\\';
doPreprocess(projectDir)
%%
projectDir = 'D:\EmoStroop\\';
excludeOutliers(projectDir)
%%
rmpath(genpath('D:\git-repos\Pupil-size-and-microsaccade-analysis'));