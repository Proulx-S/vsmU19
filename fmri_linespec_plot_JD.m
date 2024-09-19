% fmri_linespec_plot_JD.m

%Plot summary figures from analysis in fmri_linespec_JD.m

clear; clc; close all;
cd('Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW');

%% Plot Average Spectrum Power for each trial as fxn of stim freq
% Same as fxn of prox to "vaso" freq (in 50s stim trials).

% 10s
files_10s = dir('*10s*toplot.mat');
meanpwr10s = [];
meanLinePwr10s = [];
for i = 1:length(files_10s)
    clearvars toplot
    load(files_10s(i).name);
    if i == 1
        f_stim10s = toplot.f_stim;
    end
    meanpwr10s = [meanpwr10s;toplot.meanpwr_vals(1)];
    %Get line amplitudes:
    sigAmps_tmp = toplot.sig_Amps;
    meanLinePwr10s = [meanLinePwr10s;sum(sigAmps_tmp(:,1),'omitnan')/size(sigAmps_tmp,1)];
end
% 15s
files_15s = dir('*15s*toplot.mat');
meanpwr15s = [];
meanLinePwr15s = [];
for i = 1:length(files_15s)
    clearvars toplot
    load(files_15s(i).name);
    if i == 1
        f_stim15s = toplot.f_stim;
    end
    meanpwr15s = [meanpwr15s;toplot.meanpwr_vals(1)];
    %Get line amplitudes:
    sigAmps_tmp = toplot.sig_Amps;
    meanLinePwr15s = [meanLinePwr15s;sum(sigAmps_tmp(:,1),'omitnan')/size(sigAmps_tmp,1)];
end
% 20s
files_20s = dir('*20s*toplot.mat');
meanpwr20s = [];
meanLinePwr20s = [];
for i = 1:length(files_20s)
    clearvars toplot
    load(files_20s(i).name);
    if i == 1
        f_stim20s = toplot.f_stim;
    end
    meanpwr20s = [meanpwr20s;toplot.meanpwr_vals(1)];
    %Get line amplitudes:
    sigAmps_tmp = toplot.sig_Amps;
    meanLinePwr20s = [meanLinePwr20s;sum(sigAmps_tmp(:,1),'omitnan')/size(sigAmps_tmp,1)];
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

%% Plot summary figure: Average Spectrum Power for each trial as fxn of stim freq

cd('Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW\Summary');

xdata1 = repmat(f_stim10s(1),[size(meanpwr10s,1),1]);
xdata2 = repmat(f_stim15s(1),[size(meanpwr15s,1),1]);
xdata3 = repmat(f_stim20s(1),[size(meanpwr20s,1),1]);
xdata = [xdata1;xdata2;xdata3];

ydata = [meanpwr10s;meanpwr15s;meanpwr20s];
y_mean = [mean(meanpwr10s),mean(meanpwr15s),mean(meanpwr20s)];
x_mean = [1/10,1/15,1/20];
y_meanPlusStd = y_mean + [std(meanpwr10s),std(meanpwr15s),std(meanpwr20s)];
y_meanMinusStd = y_mean - [std(meanpwr10s),std(meanpwr15s),std(meanpwr20s)];
errUP = log10(y_meanPlusStd) - log10(y_mean);
errDOWN = log10(y_mean) - log10(y_meanMinusStd);
y_meanPlusSE = y_mean + [std(meanpwr10s)/sqrt(length(meanpwr10s)),std(meanpwr15s)/sqrt(length(meanpwr15s)),std(meanpwr20s)/sqrt(length(meanpwr20s))];
y_meanMinusSE = y_mean - [std(meanpwr10s)/sqrt(length(meanpwr10s)),std(meanpwr15s)/sqrt(length(meanpwr15s)),std(meanpwr20s)/sqrt(length(meanpwr20s))];
errUPSE = log10(y_meanPlusSE) - log10(y_mean);
errDOWNSE = log10(y_mean) - log10(y_meanMinusSE);

%% PLOT AND SAVE AVERAGE POWER VS STIM FREQUENCY
figure;
scatter(xdata,log10(ydata),'filled','MarkerFaceAlpha',0.3)
xlim([0 0.15]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'log10(Power) at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Power at Stimulation Frequency, Mean $\pm$ SD',sprintf('%.0f Runs, 2 Subjects',length(xdata))},'Interpreter','latex');
hold on
scatter(x_mean,log10(y_mean),'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
erbar = errorbar(x_mean,log10(y_mean),errDOWN,errUP);
erbar.Color = 'k';
savefig('PowerAtStimFreq_MeanSD.fig');
saveas(gcf,'PowerAtStimFreq_MeanSD.png');
%.eps
figure;
scatter(xdata,log10(ydata),'filled','MarkerFaceAlpha',1)
xlim([0 0.15]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'log10(Power) at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Power at Stimulation Frequency, Mean $\pm$ SD',sprintf('%.0f Runs, 2 Subjects',length(xdata))},'Interpreter','latex');
hold on
scatter(x_mean,log10(y_mean),'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
erbar = errorbar(x_mean,log10(y_mean),errDOWN,errUP);
erbar.Color = 'k';
print(gcf,'PowerAtStimFreq_MeanSD','-depsc2','-r0')

%SE
figure;
scatter(xdata,log10(ydata),'filled','MarkerFaceAlpha',0.3)
xlim([0 0.15]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'log10(Power) at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Power at Stimulation Frequency, Mean $\pm$ SE',sprintf('%.0f Runs, 2 Subjects',length(xdata))},'Interpreter','latex');
hold on
scatter(x_mean,log10(y_mean),'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
erbar = errorbar(x_mean,log10(y_mean),errDOWNSE,errUPSE);
erbar.Color = 'k';
savefig('PowerAtStimFreq_MeanSE.fig');
saveas(gcf,'PowerAtStimFreq_MeanSE.png');
%.eps
figure;
scatter(xdata,log10(ydata),'filled','MarkerFaceAlpha',1)
xlim([0 0.15]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'log10(Power) at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Power at Stimulation Frequency, Mean $\pm$ SE',sprintf('%.0f Runs, 2 Subjects',length(xdata))},'Interpreter','latex');
hold on
scatter(x_mean,log10(y_mean),'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
erbar = errorbar(x_mean,log10(y_mean),errDOWNSE,errUPSE);
erbar.Color = 'k';
print(gcf,'PowerAtStimFreq_MeanSE','-depsc2','-r0');

%% PLOT AND SAVE AVERAGE POWER VS PROX TO 50s Peak FREQUENCY
vasofreq = mean(vasopk);
fig = figure;
scatter(xdata - vasofreq,log10(ydata),'filled','MarkerFaceAlpha',0.3)
xlim([-0.15 0.0]);
xlabel('Stimulation Frequency - Vasomotor Frequency (Hz)','Interpreter','latex');
ylabel({'log10(Power) at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({['Power at Stimulation Frequency, Mean $\pm$ SD, ',sprintf('%.0f Runs, 2 Subjects',length(xdata))],sprintf('Average Vasomotor Frequency = %.2f Hz (50s Stim Period 1s Duration Runs)',vasofreq)},'Interpreter','latex');
hold on
scatter(x_mean - vasofreq,log10(y_mean),'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
erbar = errorbar(x_mean - vasofreq,log10(y_mean),errDOWN,errUP);
erbar.Color = 'k';
fig.Position = [680 508 682 490];
savefig('PowerAtStimFreq_VasoProx_MeanSD.fig');
saveas(gcf,'PowerAtStimFreq_VasoProx_MeanSD.png');
%.eps
fig = figure;
scatter(xdata - vasofreq,log10(ydata),'filled','MarkerFaceAlpha',1)
xlim([-0.15 0.0]);
xlabel('Stimulation Frequency - Vasomotor Frequency (Hz)','Interpreter','latex');
ylabel({'log10(Power) at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({['Power at Stimulation Frequency, Mean $\pm$ SD, ',sprintf('%.0f Runs, 2 Subjects',length(xdata))],sprintf('Average Vasomotor Frequency = %.2f Hz (50s Stim Period 1s Duration Runs)',vasofreq)},'Interpreter','latex');
hold on
scatter(x_mean - vasofreq,log10(y_mean),'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
erbar = errorbar(x_mean - vasofreq,log10(y_mean),errDOWN,errUP);
erbar.Color = 'k';
fig.Position = [680 508 682 490];
print(gcf,'PowerAtStimFreq_VasoProx_MeanSD','-depsc2','-r0')

%SE
fig = figure;
scatter(xdata - vasofreq,log10(ydata),'filled','MarkerFaceAlpha',0.3)
xlim([-0.15 0.0]);
xlabel('Stimulation Frequency - Vasomotor Frequency (Hz)','Interpreter','latex');
ylabel({'log10(Power) at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({['Power at Stimulation Frequency, Mean $\pm$ SE, ',sprintf('%.0f Runs, 2 Subjects',length(xdata))],sprintf('Average Vasomotor Frequency = %.2f Hz (50s Stim Period 1s Duration Runs)',vasofreq)},'Interpreter','latex');
hold on
scatter(x_mean - vasofreq,log10(y_mean),'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
erbar = errorbar(x_mean - vasofreq,log10(y_mean),errDOWNSE,errUPSE);
erbar.Color = 'k';
fig.Position = [680 508 682 490];
savefig('PowerAtStimFreq_VasoProx_MeanSE.fig');
saveas(gcf,'PowerAtStimFreq_VasoProx_MeanSE.png');
%.eps
fig = figure;
scatter(xdata - vasofreq,log10(ydata),'filled','MarkerFaceAlpha',1)
xlim([-0.15 0.0]);
xlabel('Stimulation Frequency - Vasomotor Frequency (Hz)','Interpreter','latex');
ylabel({'log10(Power) at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({['Power at Stimulation Frequency, Mean $\pm$ SE, ',sprintf('%.0f Runs, 2 Subjects',length(xdata))],sprintf('Average Vasomotor Frequency = %.2f Hz (50s Stim Period 1s Duration Runs)',vasofreq)},'Interpreter','latex');
hold on
scatter(x_mean - vasofreq,log10(y_mean),'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
erbar = errorbar(x_mean - vasofreq,log10(y_mean),errDOWNSE,errUPSE);
erbar.Color = 'k';
fig.Position = [680 508 682 490];
print(gcf,'PowerAtStimFreq_VasoProx_MeanSE','-depsc2','-r0');

%% DO THE SAME FOR THE LINE AMPLITUDE
cd('Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW\Summary');

xdata1 = repmat(f_stim10s(1),[size(meanLinePwr10s,1),1]);
xdata2 = repmat(f_stim15s(1),[size(meanLinePwr15s,1),1]);
xdata3 = repmat(f_stim20s(1),[size(meanLinePwr20s,1),1]);
xdata = [xdata1;xdata2;xdata3];

ydata = [meanLinePwr10s;meanLinePwr15s;meanLinePwr20s];
y_mean = [mean(meanLinePwr10s),mean(meanLinePwr15s),mean(meanLinePwr20s)];
x_mean = [1/10,1/15,1/20];
y_meanPlusStd = y_mean + [std(meanLinePwr10s),std(meanLinePwr15s),std(meanLinePwr20s)];
y_meanMinusStd = y_mean - [std(meanLinePwr10s),std(meanLinePwr15s),std(meanLinePwr20s)];
% errUP = log10(y_meanPlusStd) - log10(y_mean);
% errDOWN = log10(y_mean) - log10(y_meanMinusStd);
errUP = (y_meanPlusStd) - (y_mean);
errDOWN = (y_mean) - (y_meanMinusStd);
y_meanPlusSE = y_mean + [std(meanLinePwr10s)/sqrt(length(meanLinePwr10s)),std(meanLinePwr15s)/sqrt(length(meanLinePwr15s)),std(meanLinePwr20s)/sqrt(length(meanLinePwr20s))];
y_meanMinusSE = y_mean - [std(meanLinePwr10s)/sqrt(length(meanLinePwr10s)),std(meanLinePwr15s)/sqrt(length(meanLinePwr15s)),std(meanLinePwr20s)/sqrt(length(meanLinePwr20s))];
% errUPSE = log10(y_meanPlusSE) - log10(y_mean);
% errDOWNSE = log10(y_mean) - log10(y_meanMinusSE);
errUPSE = (y_meanPlusSE) - (y_mean);
errDOWNSE = (y_mean) - (y_meanMinusSE);
%% PLOT AND SAVE AVERAGE POWER VS STIM FREQUENCY
figure;
scatter(xdata,(ydata),'filled','MarkerFaceAlpha',0.3)
xlim([0 0.15]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Line Power at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency, Mean $\pm$ SD',sprintf('%.0f Runs, 2 Subjects',length(xdata))},'Interpreter','latex');
hold on
scatter(x_mean,(y_mean),'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
erbar = errorbar(x_mean,(y_mean),errDOWN,errUP);
erbar.Color = 'k';
savefig('LinePowerAtStimFreq_MeanSD.fig');
saveas(gcf,'LinePowerAtStimFreq_MeanSD.png');
%.eps
figure;
scatter(xdata,(ydata),'filled','MarkerFaceAlpha',1)
xlim([0 0.15]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Line Power at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency, Mean $\pm$ SD',sprintf('%.0f Runs, 2 Subjects',length(xdata))},'Interpreter','latex');
hold on
scatter(x_mean,(y_mean),'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
erbar = errorbar(x_mean,(y_mean),errDOWN,errUP);
erbar.Color = 'k';
print(gcf,'LinePowerAtStimFreq_MeanSD','-depsc2','-r0')

%SE
figure;
scatter(xdata,(ydata),'filled','MarkerFaceAlpha',0.3)
xlim([0 0.15]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Line Power at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency, Mean $\pm$ SE',sprintf('%.0f Runs, 2 Subjects',length(xdata))},'Interpreter','latex');
hold on
scatter(x_mean,(y_mean),'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
erbar = errorbar(x_mean,(y_mean),errDOWNSE,errUPSE);
erbar.Color = 'k';
savefig('LinePowerAtStimFreq_MeanSE.fig');
saveas(gcf,'LinePowerAtStimFreq_MeanSE.png');
%.eps
figure;
scatter(xdata,(ydata),'filled','MarkerFaceAlpha',1)
xlim([0 0.15]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Line Power at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency, Mean $\pm$ SE',sprintf('%.0f Runs, 2 Subjects',length(xdata))},'Interpreter','latex');
hold on
scatter(x_mean,(y_mean),'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
erbar = errorbar(x_mean,(y_mean),errDOWNSE,errUPSE);
erbar.Color = 'k';
print(gcf,'LinePowerAtStimFreq_MeanSE','-depsc2','-r0');
















%% Plot all significant line amplitudes - including harmonics.
clear; clc; close all;
cd('Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2');

%% Plot Average Spectrum Power for each trial as fxn of stim freq
% Same as fxn of prox to "vaso" freq (in 50s stim trials).

% 5s
files_5s = dir('*05s*toplot.mat');
meanpwr5s = [];
meanLinePwr5s = [];
hasNoSigLines5s = [];
for i = 1:length(files_5s)
    clearvars toplot
    load(files_5s(i).name);
    if i == 1
        f_stim5s = toplot.f_stim;
    end
    meanpwr5s = [meanpwr5s;toplot.meanpwr_vals(1)];
    %Get line amplitudes:
    sigAmps_tmp = toplot.sig_Amps;
    meanLinePwr5s = [meanLinePwr5s;sum(sigAmps_tmp,1,'omitnan')/size(sigAmps_tmp,1)];
    hasNoSigLines5s = [hasNoSigLines5s;sum(isnan(sigAmps_tmp),1)==size(sigAmps_tmp,1)];
end
% 6s
files_6s = dir('*06s*toplot.mat');
meanpwr6s = [];
meanLinePwr6s = [];
hasNoSigLines6s = [];
for i = 1:length(files_6s)
    clearvars toplot
    load(files_6s(i).name);
    if i == 1
        f_stim6s = toplot.f_stim;
    end
    meanpwr6s = [meanpwr6s;toplot.meanpwr_vals(1)];
    %Get line amplitudes:
    sigAmps_tmp = toplot.sig_Amps;
    meanLinePwr6s = [meanLinePwr6s;sum(sigAmps_tmp,1,'omitnan')/size(sigAmps_tmp,1)];
    hasNoSigLines6s = [hasNoSigLines6s;sum(isnan(sigAmps_tmp),1)==size(sigAmps_tmp,1)];
end
% 8s
files_8s = dir('*08s*toplot.mat');
meanpwr8s = [];
meanLinePwr8s = [];
hasNoSigLines8s = [];
for i = 1:length(files_8s)
    clearvars toplot
    load(files_8s(i).name);
    if i == 1
        f_stim8s = toplot.f_stim;
    end
    meanpwr8s = [meanpwr8s;toplot.meanpwr_vals(1)];
    %Get line amplitudes:
    sigAmps_tmp = toplot.sig_Amps;
    meanLinePwr8s = [meanLinePwr8s;sum(sigAmps_tmp,1,'omitnan')/size(sigAmps_tmp,1)];
    hasNoSigLines8s = [hasNoSigLines8s;sum(isnan(sigAmps_tmp),1)==size(sigAmps_tmp,1)];
end
% 10s
files_10s = dir('*10s*toplot.mat');
meanpwr10s = [];
meanLinePwr10s = [];
hasNoSigLines10s = [];
for i = 1:length(files_10s)
    clearvars toplot
    load(files_10s(i).name);
    if i == 1
        f_stim10s = toplot.f_stim;
    end
    meanpwr10s = [meanpwr10s;toplot.meanpwr_vals(1)];
    %Get line amplitudes:
    sigAmps_tmp = toplot.sig_Amps;
    meanLinePwr10s = [meanLinePwr10s;sum(sigAmps_tmp,1,'omitnan')/size(sigAmps_tmp,1)];
    hasNoSigLines10s = [hasNoSigLines10s;sum(isnan(sigAmps_tmp),1)==size(sigAmps_tmp,1)];
end
% 15s
files_15s = dir('*15s*toplot.mat');
meanpwr15s = [];
meanLinePwr15s = [];
hasNoSigLines15s = [];
for i = 1:length(files_15s)
    clearvars toplot
    load(files_15s(i).name);
    if i == 1
        f_stim15s = toplot.f_stim;
    end
    meanpwr15s = [meanpwr15s;toplot.meanpwr_vals(1)];
    %Get line amplitudes:
    sigAmps_tmp = toplot.sig_Amps;
    meanLinePwr15s = [meanLinePwr15s;sum(sigAmps_tmp,1,'omitnan')/size(sigAmps_tmp,1)];
    hasNoSigLines15s = [hasNoSigLines15s;sum(isnan(sigAmps_tmp),1)==size(sigAmps_tmp,1)];
end
% 20s
files_20s = dir('*20s*toplot.mat');
meanpwr20s = [];
meanLinePwr20s = [];
hasNoSigLines20s = [];
for i = 1:length(files_20s)
    clearvars toplot
    load(files_20s(i).name);
    if i == 1
        f_stim20s = toplot.f_stim;
    end
    meanpwr20s = [meanpwr20s;toplot.meanpwr_vals(1)];
    %Get line amplitudes:
    sigAmps_tmp = toplot.sig_Amps;
    meanLinePwr20s = [meanLinePwr20s;sum(sigAmps_tmp,1,'omitnan')/size(sigAmps_tmp,1)];
    hasNoSigLines20s = [hasNoSigLines20s;sum(isnan(sigAmps_tmp),1)==size(sigAmps_tmp,1)];
end
%ADD 50s PERIOD TRIALS!!!
cd('Y:\DataAnalysis\MRI\Human240904\13685568\results_01HzHalfBW_50sPeriod');
% 50s P1 P2
[meanpwr50s,meanLinePwr50s,hasNoSigLines50s,f_stim50s] = fun_getExperimentFiles('50s');

%% PLOT ALL DATA INCLUDING HARMONIC LINE AMPLITUDES
cd('Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2\Summary');

xdata1 = repmat(f_stim10s,[size(meanLinePwr10s,1),1]);
xdata2 = repmat(f_stim15s,[size(meanLinePwr15s,1),1]);
xdata3 = repmat(f_stim20s,[size(meanLinePwr20s,1),1]);
xdata4 = repmat(f_stim5s,[size(meanLinePwr5s,1),1]);
xdata5 = repmat(f_stim6s,[size(meanLinePwr6s,1),1]);
xdata6 = repmat(f_stim8s,[size(meanLinePwr8s,1),1]);
% xdata7 = repmat(f_stim50s,[size(meanLinePwr50s,1),1]);
%Delete entries where there were no significant lines detected.
todel = find(hasNoSigLines10s);
xdata1(todel') = NaN; meanLinePwr10s(todel) = NaN;
todel = find(hasNoSigLines15s);
xdata2(todel') = NaN; meanLinePwr15s(todel) = NaN;
todel = find(hasNoSigLines20s);
xdata3(todel') = NaN; meanLinePwr20s(todel) = NaN;
todel = find(hasNoSigLines5s);
xdata4(todel') = NaN; meanLinePwr5s(todel) = NaN;
todel = find(hasNoSigLines6s);
xdata5(todel') = NaN; meanLinePwr6s(todel) = NaN;
todel = find(hasNoSigLines8s);
xdata6(todel') = NaN; meanLinePwr8s(todel) = NaN;
% todel = find(hasNoSigLines50s);
% xdata7(todel') = NaN; meanLinePwr50s(todel) = NaN;


figure;
scatter(xdata1,meanLinePwr10s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b'); hold on;
scatter(xdata2,meanLinePwr15s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b'); 
scatter(xdata3,meanLinePwr20s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
scatter(xdata4,meanLinePwr5s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
scatter(xdata5,meanLinePwr6s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
scatter(xdata6,meanLinePwr8s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
% scatter(xdata7,meanLinePwr50s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
numtrials = size(meanLinePwr20s,1)+size(meanLinePwr15s,1)+size(meanLinePwr10s,1)+size(meanLinePwr8s,1)+size(meanLinePwr6s,1)+size(meanLinePwr5s,1);
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Significant Line Power','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials)},'Interpreter','latex');
savefig('LinePowerAtStimFreq_Harmonics.fig');
saveas(gcf,'LinePowerAtStimFreq_Harmonics.png');
%.eps
figure;
scatter(xdata1,meanLinePwr10s,'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b'); hold on;
scatter(xdata2,meanLinePwr15s,'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b'); 
scatter(xdata3,meanLinePwr20s,'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b');
scatter(xdata4,meanLinePwr5s,'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b');
scatter(xdata5,meanLinePwr6s,'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b');
scatter(xdata6,meanLinePwr8s,'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b');
% scatter(xdata7,meanLinePwr50s,'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b');
numtrials = size(meanLinePwr20s,1)+size(meanLinePwr15s,1)+size(meanLinePwr10s,1)+size(meanLinePwr8s,1)+size(meanLinePwr6s,1)+size(meanLinePwr5s,1);
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Significant Line Power','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials)},'Interpreter','latex');
print(gcf,'LinePowerAtStimFreq_Harmonics','-depsc2','-r0')
%%
%Plot on log scale
figure;
scatter(xdata1,log10(meanLinePwr10s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b'); hold on;
scatter(xdata2,log10(meanLinePwr15s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b'); 
scatter(xdata3,log10(meanLinePwr20s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
scatter(xdata4,log10(meanLinePwr5s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
scatter(xdata5,log10(meanLinePwr6s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
scatter(xdata6,log10(meanLinePwr8s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
% scatter(xdata7,log10(meanLinePwr50s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
numtrials = size(meanLinePwr20s,1)+size(meanLinePwr15s,1)+size(meanLinePwr10s,1)+size(meanLinePwr8s,1)+size(meanLinePwr6s,1)+size(meanLinePwr5s,1);
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Log10 Significant Line Power','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials)},'Interpreter','latex');
savefig('LinePowerAtStimFreq_Harmonics_Log10.fig');
saveas(gcf,'LinePowerAtStimFreq_Harmonics_Log10.png');
%.eps
figure;
scatter(xdata1,log10(meanLinePwr10s),'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b'); hold on;
scatter(xdata2,log10(meanLinePwr15s),'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b'); 
scatter(xdata3,log10(meanLinePwr20s),'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b');
scatter(xdata4,log10(meanLinePwr5s),'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b');
scatter(xdata5,log10(meanLinePwr6s),'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b');
scatter(xdata6,log10(meanLinePwr8s),'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b');
% scatter(xdata7,log10(meanLinePwr50s),'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b');
numtrials = size(meanLinePwr20s,1)+size(meanLinePwr15s,1)+size(meanLinePwr10s,1)+size(meanLinePwr8s,1)+size(meanLinePwr6s,1)+size(meanLinePwr5s,1);
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Log10 Significant Line Power','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials)},'Interpreter','latex');
savefig('LinePowerAtStimFreq_Harmonics_Log10.fig');
%% Spline fit!
xtmp = [xdata1(~isnan(xdata1));xdata2(~isnan(xdata2));xdata3(~isnan(xdata3));...
    xdata4(~isnan(xdata4));xdata5(~isnan(xdata5));xdata6(~isnan(xdata6))];
ytmp = [meanLinePwr10s(~isnan(meanLinePwr10s));meanLinePwr15s(~isnan(meanLinePwr15s));meanLinePwr20s(~isnan(meanLinePwr20s));...
    meanLinePwr5s(~isnan(meanLinePwr5s));meanLinePwr6s(~isnan(meanLinePwr6s));meanLinePwr8s(~isnan(meanLinePwr8s))];
% Spline fit
mtd = 'median';
movAvgNum = 15;
[xSort,~,~,yPct,~,yOutFull] = fun_MovingAvgSD(xtmp,ytmp,movAvgNum,mtd,[25,75]);
%Plot!
figure; scatter(xtmp,ytmp,'filled','MarkerFaceAlpha',0.2)
hold on
plot(xSort,yOutFull,'k')
plot(xSort,yPct,'k')
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Significant Line Power','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials),'15 Point Median 25th 75th percentile'},'Interpreter','latex');
savefig('LinePowerAtStimFreq_Harmonics_15ptMovMedian.fig');
saveas(gcf,'LinePowerAtStimFreq_Harmonics_15ptMovMedian.png');
%.eps
figure; scatter(xtmp,ytmp,'filled','MarkerFaceAlpha',1)
hold on
plot(xSort,yOutFull,'k')
plot(xSort,yPct,'k')
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Significant Line Power','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials),'15 Point Median 25th 75th percentile'},'Interpreter','latex');
print(gcf,'LinePowerAtStimFreq_Harmonics_15ptMovMedian','-depsc2','-r0')
%% LOG SCALE
figure; scatter(xtmp,log10(ytmp),'filled','MarkerFaceAlpha',0.2)
hold on
plot(xSort,log10(yOutFull),'k')
plot(xSort,log10(yPct),'k')
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Log10 Significant Line Power','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials),'15 Point Median 25th 75th percentile'},'Interpreter','latex');
savefig('LinePowerAtStimFreq_Harmonics_log10_15ptMovMedian.fig');
saveas(gcf,'LinePowerAtStimFreq_Harmonics_log10_15ptMovMedian.png');
%.eps
figure; scatter(xtmp,log10(ytmp),'filled','MarkerFaceAlpha',1)
hold on
plot(xSort,log10(yOutFull),'k')
plot(xSort,log10(yPct),'k')
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Log10 Significant Line Power','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials),'15 Point Median 25th 75th percentile'},'Interpreter','latex');
print(gcf,'LinePowerAtStimFreq_Harmonics_log10_15ptMovMedian','-depsc2','-r0')

%% SPLINE fits
mtd = 'median';
movAvgNum = 15;
[xSort,~,~,yPct,~,yOutFull] = fun_MovingAvgSD(xtmp,ytmp,movAvgNum,mtd,[25,75]);
% FIT SPLINES!
[xFit, yFit] = prepareCurveData( xSort, yOutFull );
% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
opts.SmoothingParam = 0.99;
% Fit model to data.
[fitresultmed, gof] = fit( xFit, yFit, ft, opts );
fdata = feval(fitresultmed,xFit);
figure;
scatter(xtmp,ytmp,'filled','MarkerFaceAlpha',0.3,'MarkerFaceColor',[0,0,1]);
alpha = 0.25;
hold on
plot(xFit,fdata,'Color',[0,0,0,1],'LineWidth',1);
f1 = gcf;

[xFit, yFit] = prepareCurveData( xSort, yPct(:,2) );
% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
smparamvals = [0.99,0.999,0.9999,0.99999];
for smparam = smparamvals
    opts.SmoothingParam = smparam;
    [fitresult, ~] = fit( xFit, yFit, ft, opts );
    figure;
    plot(fitresult,xFit,yFit)
    splinein = input('Sufficent Spline? 1 yes ');
    if splinein == 1
        useparam = smparam;
    end
end
opts.SmoothingParam = useparam;
% Fit model to data.
[fitresultUP, ~] = fit( xFit, yFit, ft, opts );
xFitPlot = xFit(1):0.001:xFit(end);
fdataUP = feval(fitresultUP,xFitPlot);
figure(f1)
hold on
% plot(xFitPlot,fdataUP,'Color',[0,0,0,alpha]);

[xFit, yFit] = prepareCurveData( xSort, yPct(:,1) );
% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
smparamvals = [0.99,0.999,0.9999,0.99999];
for smparam = smparamvals
    opts.SmoothingParam = smparam;
    [fitresult, ~] = fit( xFit, yFit, ft, opts );
    figure;
    plot(fitresult,xFit,yFit)
    splinein = input('Sufficent Spline? 1 yes ');
    if splinein == 1
        useparam = smparam;
    end
end
opts.SmoothingParam = useparam;
% Fit model to data.
[fitresultDOWN, gof] = fit( xFit, yFit, ft, opts );
xFitPlot = xFit(1):0.001:xFit(end);
fdataDOWN = feval(fitresultDOWN,xFitPlot);
%%

%Plot!
figure; scatter(xtmp,ytmp,'filled','MarkerFaceAlpha',0.2)
hold on;
x2 = [xFitPlot, fliplr(xFitPlot)];
inBetween = [fdataDOWN', fliplr(fdataUP')];
fill(x2, inBetween,'k','FaceAlpha',0.25,'EdgeAlpha',0);
plot(xFit,fdata,'Color',[0,0,0,1],'LineWidth',1);

xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Significant Line Power','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials),'15 Point Median 25th 75th percentile'},'Interpreter','latex');
savefig('LinePowerAtStimFreq_Harmonics_15ptMovMedian_Spline.fig');
saveas(gcf,'LinePowerAtStimFreq_Harmonics_15ptMovMedian_Spline.png');
%.eps
figure; scatter(xtmp,ytmp,'filled','MarkerFaceAlpha',1)
hold on;
plot(xFit,fdata,'Color','k');
plot(xFitPlot,fdataDOWN,'Color','k');
plot(xFitPlot,fdataUP,'Color','k');
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Significant Line Power','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials),'15 Point Median 25th 75th percentile'},'Interpreter','latex');
print(gcf,'LinePowerAtStimFreq_Harmonics_15ptMovMedian_Spline','-depsc2','-r0')
%LOG SCALE
figure; scatter(xtmp,log10(ytmp),'filled','MarkerFaceAlpha',0.2)
hold on;
x2 = [xFitPlot, fliplr(xFitPlot)];
inBetween = [fdataDOWN', fliplr(fdataUP')];
fill(x2, log10(inBetween),'k','FaceAlpha',0.25,'EdgeAlpha',0);
plot(xFit,log10(fdata),'Color',[0,0,0,1],'LineWidth',1);
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'log10 Significant Line Power','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials),'15 Point Median 25th 75th percentile'},'Interpreter','latex');
savefig('LinePowerAtStimFreq_Harmonics_log10_15ptMovMedian_Spline.fig');
saveas(gcf,'LinePowerAtStimFreq_Harmonics_log10_15ptMovMedian_Spline.png');
%.eps
figure; scatter(xtmp,log10(ytmp),'filled','MarkerFaceAlpha',1)
hold on;
plot(xFit,log10(fdata),'Color','k');
plot(xFitPlot,log10(fdataDOWN),'Color','k');
plot(xFitPlot,log10(fdataUP),'Color','k');
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'log10 Significant Line Power','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials),'15 Point Median 25th 75th percentile'},'Interpreter','latex');
print(gcf,'LinePowerAtStimFreq_Harmonics_log10_15ptMovMedian_Spline','-depsc2','-r0')


%% Example spectra
clear; clc; close all;
% cd('Y:\DataAnalysis\MRI\Human240904\13685568\SpectraExamples')
load("Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2\P1_task_eyeOpenRest_run1_toplot.mat");
figure
plot(toplot.f,log10(toplot.avgpowr),'k');
xlabel('Frequency (Hz)','Interpreter','latex');
ylabel('log10 Power','Interpreter','latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;
title({'Resting-State Spectrum Half-BW = 0.025','Average over voxels'},'Interpreter','latex');
print(gcf,'SubjectP1_RestingStateSpectrum_BW025','-depsc2','-r0')
saveas(gcf,'SubjectP1_RestingStateSpectrum_BW025.png');

%Stim Trial Examples
clear; clc; close all;
load("Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2\P1_task_10sPrd1sDur_run3_toplot.mat")
figure
plot(toplot.f,log10(toplot.avgpowr),'k');
hold on;
plot(toplot.resid_f,log10(toplot.avgResid),'m--');
xlabel('Frequency (Hz)','Interpreter','latex');
ylabel('log10 Power','Interpreter','latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;
legend({'Average Spectrum','Residual Spectrum'})
title({'0.1Hz Visual Stimulation Spectrum Half-BW = 0.025','Average over voxels'},'Interpreter','latex');
print(gcf,'SubjectP1_10sPeriodStim_run3_Spectrum','-depsc2','-r0')
saveas(gcf,'SubjectP1_10sPeriodStim_run3_Spectrum.png');

clear; clc; close all;
load("Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2\P1_task_15sPrd1sDur_run1_toplot.mat");
figure
plot(toplot.f,log10(toplot.avgpowr),'k');
hold on;
plot(toplot.resid_f,log10(toplot.avgResid),'m--');
xlabel('Frequency (Hz)','Interpreter','latex');
ylabel('log10 Power','Interpreter','latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;
legend({'Average Spectrum','Residual Spectrum'})
title({'0.07Hz Visual Stimulation Spectrum Half-BW = 0.025','Average over voxels'},'Interpreter','latex');
print(gcf,'SubjectP1_15sPeriodStim_run1_Spectrum','-depsc2','-r0')
saveas(gcf,'SubjectP1_15sPeriodStim_run1_Spectrum.png');

clear; clc; close all;
load("Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2\P1_task_05sPrd1sDur_run2_toplot.mat");
figure
plot(toplot.f,log10(toplot.avgpowr),'k');
hold on;
plot(toplot.resid_f,log10(toplot.avgResid),'m--');
xlabel('Frequency (Hz)','Interpreter','latex');
ylabel('log10 Power','Interpreter','latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;
legend({'Average Spectrum','Residual Spectrum'})
title({'0.2Hz Visual Stimulation Spectrum Half-BW = 0.025','Average over voxels'},'Interpreter','latex');
print(gcf,'SubjectP1_5sPeriodStim_run2_Spectrum','-depsc2','-r0')
saveas(gcf,'SubjectP1_5sPeriodStim_run2_Spectrum.png');
%%
%Single-ovxel residual spectra
clear; clc; close all;
% load("Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2\P1_task_10sPrd1sDur_run3_toplot.mat")
load("Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2\P1_task_eyeOpenRest_run1_toplot.mat");
% figure
% for i = 1:size(toplot.Sresid,2)
%     plot(toplot.resid_f,log10(toplot.Stot(:,i)),'k');
%     hold on;
%     plot(toplot.resid_f,log10(toplot.Sresid(:,i)),'m--');
%     xlabel('Frequency (Hz)','Interpreter','latex');
%     ylabel('log10 Power','Interpreter','latex');
%     title(sprintf('%.0f',i));
%     ylim([3.4 4.6])
%     pause()
%     hold off
% end
% vox = 19;
% vox = 21;
% vox = 30;
vox = 34;
% vox = 16;
figure
plot(toplot.resid_f,log10(toplot.Stot(:,vox)),'k');
hold on;
plot(toplot.resid_f,log10(toplot.Sresid(:,vox)),'m--');
xlabel('Frequency (Hz)','Interpreter','latex');
ylabel('log10 Power','Interpreter','latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;
legend({'Average Spectrum','Residual Spectrum'})
title({'Rest Spectrum Half-BW = 0.025'},'Interpreter','latex');
% title({'0.1Hz Visual Stimulation Spectrum Half-BW = 0.025'},'Interpreter','latex');
% print(gcf,['SubjectP1_10sPeriodStim_run3_Voxel',num2str(vox),'_NOLine'],'-depsc2','-r0')
% saveas(gcf,['SubjectP1_10sPeriodStim_run3_Voxel',num2str(vox),'_NOLine.png']);
print(gcf,['SubjectP1_REST_Voxel',num2str(vox),'_Spectrum'],'-depsc2','-r0')
saveas(gcf,['SubjectP1_REST_Voxel',num2str(vox),'_Spectrum.png']);

%%
%Example time series
clear; clc; close all;
cd('Y:\DataAnalysis\MRI\Human240904\13685568');
dataDir = fullfile(pwd,'datav2');
dataFile = 'vsmDrivenP1.mat';
% dataFile = 'vsmDrivenP2.mat';
disp(['processing ',dataFile]);
subj = extractBetween(dataFile,'Driven','.mat');
subj = subj{1};
load(fullfile(dataDir,dataFile))
trials = fields(vfMRI);
trial = 9
toplot = struct(); %Structure to save results
toplot.trial = trial;
SVD_Q = 1; %Perform space-time SVD before line-spec?
trialName = trials{trial};
disp(['processing ',trialName])
vfMRI_tmp = vfMRI.(trialName);
toplot.SVD_Q = SVD_Q;
toplot.sub = vfMRI_tmp.sub;
toplot.label = vfMRI_tmp.label;
toplot.ses = vfMRI_tmp.ses;
run = 1
disp(['processing ',trialName,' run ',num2str(run),' of ',num2str(length(vfMRI_tmp.volTs))])
data = vfMRI_tmp.volTs(run).mri.vec;
load("Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2\P1_task_eyeOpenRest_run1_toplot.mat");
%Plot
figure; 
plot(toplot.Tvec,mean(data,2),'k');
xlabel('Time (seconds)','Interpreter','latex');
ylabel({'Signal','Average over all voxels'},'Interpreter','latex')
title('Rest (no task)','Interpreter','latex')
print(gcf,['SubjectP1_Rest_Timeseries'],'-depsc2','-r0')
saveas(gcf,['SubjectP1_Rest_Timeseries.png']);


clear; clc; close all;
cd('Y:\DataAnalysis\MRI\Human240904\13685568');
dataDir = fullfile(pwd,'datav2');
dataFile = 'vsmDrivenP1.mat';
% dataFile = 'vsmDrivenP2.mat';
disp(['processing ',dataFile]);
subj = extractBetween(dataFile,'Driven','.mat');
subj = subj{1};
load(fullfile(dataDir,dataFile))
trials = fields(vfMRI);
trial = 4
toplot = struct(); %Structure to save results
toplot.trial = trial;
SVD_Q = 1; %Perform space-time SVD before line-spec?
trialName = trials{trial};
disp(['processing ',trialName])
vfMRI_tmp = vfMRI.(trialName);
toplot.SVD_Q = SVD_Q;
toplot.sub = vfMRI_tmp.sub;
toplot.label = vfMRI_tmp.label;
toplot.ses = vfMRI_tmp.ses;
run = 3
disp(['processing ',trialName,' run ',num2str(run),' of ',num2str(length(vfMRI_tmp.volTs))])
data = vfMRI_tmp.volTs(run).mri.vec;
load("Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2\P1_task_10sPrd1sDur_run3_toplot.mat")
%Plot
figure; 
plot(toplot.Tvec,mean(data,2),'k');
xline(toplot.Stimvec)
xlabel('Time (seconds)','Interpreter','latex');
ylabel({'Signal','Average over all voxels'},'Interpreter','latex')
title('0.1Hz Stimulation','Interpreter','latex')
print(gcf,['SubjectP1_10sPeriodStim_Run3'],'-depsc2','-r0')
saveas(gcf,['SubjectP1_10sPeriodStim_Run3.png']);




% To do: Spline fit for large summary plot 
% Not all stims plot, EXCLUDE 50s STIM TRIALS
% example spectra for some trials, rest etc.
%% Subject-specific plot
%Example P1_task_05sPrd1sDur_run1_toplot.mat
clear; clc; close all;
cd('Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2');

% 5s P1 only
[meanpwr5s,meanLinePwr5s,hasNoSigLines5s,f_stim5s] = fun_getExperimentFiles('P1_task_05s');
% 6s P1 only
[meanpwr6s,meanLinePwr6s,hasNoSigLines6s,f_stim6s] = fun_getExperimentFiles('P1_task_06s');
% 8s P1 only
[meanpwr8s,meanLinePwr8s,hasNoSigLines8s,f_stim8s] = fun_getExperimentFiles('P1_task_08s');
% 10s P1 and P2
[meanpwr10sP1,meanLinePwr10sP1,hasNoSigLines10sP1,f_stim10sP1] = fun_getExperimentFiles('P1_task_10s');
[meanpwr10sP2,meanLinePwr10sP2,hasNoSigLines10sP2,f_stim10sP2] = fun_getExperimentFiles('P2_task_10s');
% 15s P1 and P2
[meanpwr15sP1,meanLinePwr15sP1,hasNoSigLines15sP1,f_stim15sP1] = fun_getExperimentFiles('P1_task_15s');
[meanpwr15sP2,meanLinePwr15sP2,hasNoSigLines15sP2,f_stim15sP2] = fun_getExperimentFiles('P2_task_15s');
% 20s P1 and P2
[meanpwr20sP1,meanLinePwr20sP1,hasNoSigLines20sP1,f_stim20sP1] = fun_getExperimentFiles('P1_task_15s');
[meanpwr20sP2,meanLinePwr20sP2,hasNoSigLines20sP2,f_stim20sP2] = fun_getExperimentFiles('P2_task_15s');

%Organize for plotting.
cd('Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2\Summary');

xdata1 = repmat(f_stim10s,[size(meanLinePwr10s,1),1]);
xdata2 = repmat(f_stim15s,[size(meanLinePwr15s,1),1]);
xdata3 = repmat(f_stim20s,[size(meanLinePwr20s,1),1]);
xdata4 = repmat(f_stim5s,[size(meanLinePwr5s,1),1]);
xdata5 = repmat(f_stim6s,[size(meanLinePwr6s,1),1]);
xdata6 = repmat(f_stim8s,[size(meanLinePwr8s,1),1]);
%Delete entries where there were no significant lines detected.
todel = find(hasNoSigLines10s);
xdata1(todel') = NaN; meanLinePwr10s(todel) = NaN;
todel = find(hasNoSigLines15s);
xdata2(todel') = NaN; meanLinePwr15s(todel) = NaN;
todel = find(hasNoSigLines20s);
xdata3(todel') = NaN; meanLinePwr20s(todel) = NaN;
todel = find(hasNoSigLines5s);
xdata4(todel') = NaN; meanLinePwr5s(todel) = NaN;
todel = find(hasNoSigLines6s);
xdata5(todel') = NaN; meanLinePwr6s(todel) = NaN;
todel = find(hasNoSigLines8s);
xdata6(todel') = NaN; meanLinePwr8s(todel) = NaN;


figure;
scatter(xdata1,meanLinePwr10s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b'); hold on;
scatter(xdata2,meanLinePwr15s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b'); 
scatter(xdata3,meanLinePwr20s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
scatter(xdata4,meanLinePwr5s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
scatter(xdata5,meanLinePwr6s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
scatter(xdata6,meanLinePwr8s,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
numtrials = size(meanLinePwr20s,1)+size(meanLinePwr15s,1)+size(meanLinePwr10s,1)+size(meanLinePwr8s,1)+size(meanLinePwr6s,1)+size(meanLinePwr6s,1);
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Line Power at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials)},'Interpreter','latex');
savefig('LinePowerAtStimFreq_Harmonics.fig');
saveas(gcf,'LinePowerAtStimFreq_Harmonics.png');








%%

% ydata = [meanLinePwr10s;meanLinePwr15s;meanLinePwr20s];
y_mean = [mean(meanLinePwr10s),mean(meanLinePwr15s),mean(meanLinePwr20s)];
x_mean = [1/10,1/15,1/20];
y_meanPlusStd = y_mean + [std(meanLinePwr10s),std(meanLinePwr15s),std(meanLinePwr20s)];
y_meanMinusStd = y_mean - [std(meanLinePwr10s),std(meanLinePwr15s),std(meanLinePwr20s)];
% errUP = log10(y_meanPlusStd) - log10(y_mean);
% errDOWN = log10(y_mean) - log10(y_meanMinusStd);
errUP = (y_meanPlusStd) - (y_mean);
errDOWN = (y_mean) - (y_meanMinusStd);
y_meanPlusSE = y_mean + [std(meanLinePwr10s)/sqrt(length(meanLinePwr10s)),std(meanLinePwr15s)/sqrt(length(meanLinePwr15s)),std(meanLinePwr20s)/sqrt(length(meanLinePwr20s))];
y_meanMinusSE = y_mean - [std(meanLinePwr10s)/sqrt(length(meanLinePwr10s)),std(meanLinePwr15s)/sqrt(length(meanLinePwr15s)),std(meanLinePwr20s)/sqrt(length(meanLinePwr20s))];
% errUPSE = log10(y_meanPlusSE) - log10(y_mean);
% errDOWNSE = log10(y_mean) - log10(y_meanMinusSE);
errUPSE = (y_meanPlusSE) - (y_mean);
errDOWNSE = (y_mean) - (y_meanMinusSE);






%%

% y_std = [std(meanpwr10s),std(meanpwr15s),std(meanpwr20s)];
% figure;
% scatter(xdata,(ydata),'filled','MarkerFaceAlpha',0.3)
% xlim([0 0.15]);
% xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
% ylabel({'Power at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% ax.FontSize = 13;
% title({'Power at Stimulation Frequency',sprintf('%.0f Trials, 2 Subjects',length(xdata))},'Interpreter','latex');
% hold on
% scatter(x_mean,y_mean,'filled','Marker','_','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
% errorbar(x_mean,y_mean,y_std)
% set(ax,'YScale','log');
% ylim([0 2*10^4])

