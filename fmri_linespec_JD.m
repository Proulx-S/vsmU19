% fmri_linespec_JD.m

% Analysis code for line spectral analysis on human (single vessel)
% stimulated fmri datasets. 

%% Load data
clear; clc; close all;
cd('Y:\DataAnalysis\MRI\Human240904\13685568');
dataDir = fullfile(pwd,'data');
dataFile = 'vsmDrivenP1.mat';
disp(['processing ',dataFile]);
subj = extractBetween(dataFile,'Driven','.mat');
subj = subj{1};
load(fullfile(dataDir,dataFile))
%% 





%%



