function saveSession(destDir)
% saveSession  Save all workspace variables and open figures.
%
%   saveSession()            % 存到目前資料夾
%   saveSession(destDir)     % 指定資料夾

    if nargin < 1 || isempty(destDir)
        destDir = pwd;                        % 預設：目前目錄
    end
    if ~exist(destDir, 'dir'); mkdir(destDir); end

    % 時間戳
    ts = datestr(now, 'yyyymmdd_HHMMSS');

    %% ── 1. 儲存 Workspace 變數 ─────────────────────────────
    matFile = fullfile(destDir, ['workspace_' ts '.mat']);
    save(matFile);      % 預設 -v7；需大檔可用 save(matFile,'-v7.3')
    fprintf('[saveSession] Workspace saved → %s\n', matFile);

    %% ── 2. 儲存所有開啟的 Figure ──────────────────────────
    figs = findall(0, 'Type', 'figure');
    if isempty(figs)
        fprintf('[saveSession] No open figures.\n');
    else
        for k = 1:numel(figs)
            fig = figs(k);
            fname = sprintf('fig%02d_%s', k, ts);
            figFile = fullfile(destDir, [fname '.fig']);
            pngFile = fullfile(destDir, [fname '.png']);

            savefig(fig, figFile);            % 可重新開啟編輯
            saveas(fig, pngFile);             % 快速預覽用
            fprintf('[saveSession] Figure %d saved → %s / %s\n', ...
                    k, figFile, pngFile);
        end
    end
    fprintf('[saveSession] Done.\n');
end
