clear all
close all

%% Download data from zenodo and put in dataDir
dataDir = fullfile(pwd,'data','bloodfest2024human');
if ~exist(dataDir,'dir'); mkdir(dataDir); end

dir(dataDir)


%% Load data (one file per subject)
load(fullfile(dataDir,'vsmDrivenP1.mat'))

fields(vfMRI)                       % substructures are different stimulus conditions

vfMRI.task_10sPrd1sDur.sub          % subject label
vfMRI.task_10sPrd1sDur.label        % stimulus condition label
vfMRI.task_10sPrd1sDur.ses          % session label (1 per acquisition run)

vfMRI.task_10sPrd1sDur.volTs(:).mri % mri data structure (1 per acquisition run)


%% THE data
run = 1;
vfMRI.task_10sPrd1sDur.volTs(run).mri.vec              % timeseries data [time x voxel] (vessel voxels only)
vfMRI.task_10sPrd1sDur.volTs(run).mri.vecInfo
figure('WindowStyle','docked');
imagesc(vfMRI.task_10sPrd1sDur.volTs(run).mri.imMean)  % timeaveraged data (whole slice)
imagesc(vfMRI.task_10sPrd1sDur.volTs(run).mri.vol2vec) % vessel mask (vol2vec) used to convert whole-slice data (vol) to vectorize data (vec)


%% Extra info
vfMRI.task_10sPrd1sDur.volTs(run).mri.tr            % 1 / sampling rate (ms)
vfMRI.task_10sPrd1sDur.volTs(run).mri.t             % acquisition time in seconds relative to the first time point of an acquisition run (note that the first few timepoints are already excluded)
vfMRI.task_10sPrd1sDur.acqTime                      % acquisition time of each run relative to the first run
vfMRI.task_10sPrd1sDur.dsgn                         % stimulus information
vfMRI.task_10sPrd1sDur.dsgn.onsetList               % stimulus onset times in seconds (same timeframe as vfMRI.task_10sPrd1sDur.volTs(run).mri.t)
vfMRI.task_10sPrd1sDur.dsgn.ondurList               % stimulus durations in seconds
1/mean(diff(vfMRI.task_10sPrd1sDur.dsgn.onsetList)) % stimulus driving frequency
vfMRI.task_10sPrd1sDur.dsgn.nullTrial               % one stimulus presentation (trial) was omitted in the middle of each run. This was meant to reproduce Das's experiment but was ignored so far (data were analyzed as if there were no trial omissions)
vfMRI.task_10sPrd1sDur.dsgn.dt                      % time frame in seconds of stimulus presentation (not exactly the same as data acquisition, but close enough)
vfMRI.task_10sPrd1sDur.dsgn.dt - vfMRI.task_10sPrd1sDur.volTs(run).mri.tr/1000




