% fmri_linespec_plot_JD.m

%Plot summary figures from analysis in fmri_linespec_JD.m

clear; clc; close all;
cd('Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW');

%% Plot Average Spectrum Power for each trial as fxn of stim freq
% Same as fxn of prox to "vaso" freq (in 50s stim trials).

% 10s
files_10s = dir('*10s*toplot.mat');
meanpwr = [];
for i = 1:length(files_10s)
    clearvars toplot
    load(files_10s(i).name);
    if i == 1
        f_stim10s = toplot.f_stim;
    end
    meanpwr10s = [meanpwr;toplot.meanpwr_vals];
end
% 15s
files_15s = dir('*15s*toplot.mat');
meanpwr = [];
for i = 1:length(files_15s)
    clearvars toplot
    load(files_15s(i).name);
    if i == 1
        f_stim15s = toplot.f_stim;
    end
    meanpwr15s = [meanpwr;toplot.meanpwr_vals];
end
% 20s
files_20s = dir('*20s*toplot.mat');
meanpwr = [];
for i = 1:length(files_20s)
    clearvars toplot
    load(files_20s(i).name);
    if i == 1
        f_stim20s = toplot.f_stim;
    end
    meanpwr20s = [meanpwr;toplot.meanpwr_vals];
end

%% Define vasomotor frequency from 50s stim 1s duration trials

% 50s 1s duration
files_50s = dir('*50sPrd1s*toplot.mat');

vasopk = zeros(length(files_50s),1);
for i = 1:length(files_50s)
    clearvars toplot
    load(files_50s(i).name);
    figure;
    plot(log10(toplot.avgpowr));
    [wind,~] = ginput(2);
    wind = fix(wind);
    [~, fpeak1] = max(toplot.avgpowr(1, wind(1):wind(2)), [], 2);
    vasopk(i) = toplot.f(fpeak1 + wind(1) - 1);
end

