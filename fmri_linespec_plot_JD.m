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

%% PLOT ALL DATA INCLUDING HARMONIC LINE AMPLITUDES
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

%Plot on log scale
figure;
scatter(xdata1,log10(meanLinePwr10s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b'); hold on;
scatter(xdata2,log10(meanLinePwr15s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b'); 
scatter(xdata3,log10(meanLinePwr20s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
scatter(xdata4,log10(meanLinePwr5s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
scatter(xdata5,log10(meanLinePwr6s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
scatter(xdata6,log10(meanLinePwr8s),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','b');
numtrials = size(meanLinePwr20s,1)+size(meanLinePwr15s,1)+size(meanLinePwr10s,1)+size(meanLinePwr8s,1)+size(meanLinePwr6s,1)+size(meanLinePwr6s,1);
xlim([0 0.6]);
xlabel('Stimulation Frequency (Hz)','Interpreter','latex');
ylabel({'Log10 Line Power at Stimulation Frequency','Average Over All Vessel Voxels'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
title({'Line Power at Stimulation Frequency and Harmonics',sprintf('%.0f Runs, 2 Subjects',numtrials)},'Interpreter','latex');
savefig('LinePowerAtStimFreq_Harmonics_Log10.fig');
saveas(gcf,'LinePowerAtStimFreq_Harmonics_Log10.png');



















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

