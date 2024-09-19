
function [meanpwr,meanLinePwr,hasNoSigLines,f_stim] = fun_getExperimentFiles(time_str)

% 5s
files = dir(['*',time_str,'*toplot.mat']);
meanpwr = [];
meanLinePwr = [];
hasNoSigLines = [];
for i = 1:length(files)
    clearvars toplot
    load(files(i).name);
    if i == 1
        f_stim = toplot.f_stim;
    end
    meanpwr = [meanpwr;toplot.meanpwr_vals(1)];
    %Get line amplitudes:
    sigAmps_tmp = toplot.sig_Amps;
    meanLinePwr = [meanLinePwr;sum(sigAmps_tmp,1,'omitnan')/size(sigAmps_tmp,1)];
    hasNoSigLines = [hasNoSigLines;sum(isnan(sigAmps_tmp),1)==size(sigAmps_tmp,1)];
end