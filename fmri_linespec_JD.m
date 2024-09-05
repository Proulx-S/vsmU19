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
%% Analyze
trials = fields(vfMRI);
% for trial = 1:length(trials)
for trial = 1:1
    output = struct(); %Structure to save results
    SVD_Q = 1; %Perform space-time SVD before line-spec?

    trialName = trials{trial};
    disp(['processing ',trialName])

    vfMRI_tmp = vfMRI.(trialName);
    output.SVD_Q = SVD_Q;
    output.sub = vfMRI_tmp.sub;
    output.label = vfMRI_tmp.label;
    output.ses = vfMRI_tmp.ses;
    % for run = 1:length(vfMRI_tmp.volTs)
    for run = 1:1
        disp(['processing ',trialName,' run ',num2str(run),' of ',num2str(length(vfMRI_tmp.volTs))])
        data = vfMRI_tmp.volTs(run).mri.vec;
        if size(data,1) > size(data,2)
            data = data';
        end
        disp(['Data matrix is ',num2str(size(data,1)),' by ',num2str(size(data,2))]) %Make sure this is space x time for following analysis.
        %DEFINE TIME FOR THE TRIAL
        Tms = vfMRI_tmp.volTs(run).mri.tr; % 1 / sampling rate (ms)
        Fs = 1/Tms*1000; % Sampling rate (Hz)
        Tvec = vfMRI_tmp.volTs(run).mri.t; % acquisition time in seconds relative to the first time point of an acquisition run (note that the first few timepoints are already excluded)
        Toffset = Tvec(1);
        Tvec = Tvec - Toffset;
        Stimvec = vfMRI_tmp.dsgn.onsetList';
        stim_omit = vfMRI_tmp.dsgn.nullTrial;
        Stimvec(stim_omit) = [];

        output.Fs = Fs;
        output.Tvec = Tvec;
        output.Stimvec = Stimvec;
        figure('WindowStyle','docked');
        plot(diff(Tvec)); ylim([0 2]);
        xlabel('Time (s)','Interpreter','latex');
        ylabel('Difference between acquisition frames','Interpreter','latex');

        %SUBTRACT MEAN FROM TIMESERIES
        [data_mean] = fun_MeanSubtract(data);

        if SVD_Q == 1
            [U,S,V]=svd(data_mean,0);
            figure('WindowStyle','docked');
            plot(log10(diag(S)));
            xlabel('Mode','Interpreter','latex');
            ylabel('Singular Value','Interpreter','latex');

            sig_modes = round(size(data,1)/2); %Start with only half the modes as test.
            output.sig_modes = sig_modes;
            Un=single(U(:,1:sig_modes));
            Sn=single(S(1:sig_modes,1:sig_modes));
            Vn=single(V(:,1:sig_modes));

            clearvars U S V
            figure('WindowStyle','docked');
            plot(Tvec,Vn(:,1),'k'); hold on;
            plot(Tvec,Vn(:,2),'b');
            plot(Tvec,Vn(:,3),'g');
            xline(Stimvec,'Color','k','Alpha',0.2);
            xlabel('Time (s)','Interpreter','latex');
            ylabel('Temporal mode value','Interpreter','latex');
            legend({'Mode 1','Mode 2','Mode 3'})

            %CALCULATE AVEAGE SPECTRUM
            num_pixel = size(Un,1);
            num_frame= size(Vn,1);
            Delta_f = 0.02; %Half-Bandwidth. 0.02-> ~11 tapers | 0.03 -> ~17 tapers | 0.04-> ~22 tapers 
            padding_ratio = 2; 
            num_frame_pad = (2 ^ ceil(log2(num_frame))) * padding_ratio;
            % toplot.pad = num_frame_pad;
            p = round(num_frame / Fs * Delta_f); % Time BW product
            % Delta_f = p * toplot.rate / num_frame;
            disp(['Bandwidth = ', num2str(Delta_f), ' Hz'])
            num_tapers = 2 * p - 1;
            [slep,~] = dpss(num_frame, p, num_tapers);
            %% Perform Spectral FFT
            % Update half-bandwidth according to the rounded p value
            Delta_f = p * Fs / num_frame;
            addpath(genpath('C:\chronux_2_12'))
            ntapers = [(num_tapers+1)/2,num_tapers];
            % Fs = toplot.rate;
            nfft = num_frame_pad;
            tapers = dpsschk(ntapers,size(Vn,1),Fs);
            [f,findx] = getfgrid(Fs,nfft,[0,Fs/2]);

            %FFT
            taperedFFT = complex(zeros(length(f),ntapers(2),size(Vn,2)));
            for i = 1:size(Vn,2) % Iterate over modes
                J = mtfftc(Vn(:,i),tapers,nfft,Fs);
                J = J(findx,:,:);
                taperedFFT(:,:,i) = J;
            end
            scores = Un*Sn;
            %AVERAGE ACROSS TAPERS
            S_tot = zeros(size(scores,1),size(taperedFFT,1),'single'); %S_tot is space x freq
            for k = 1:ntapers(2) %Iterate over tapers
                z = scores*squeeze(taperedFFT(:,k,:))'; %Reconstruct each tapered data's FT
                S_tot = S_tot + conj(z).*z; %Average Spectrum
            end
            S_tot = S_tot./ntapers(2); %Spectrum (Averaged across tapers) for each pixel
            % clear Fs ntapers nfft tapers findx
            % [num_pixel] = size(toplot.mask_ind,1);
            num_frame_pad=size(S_tot,2); %Size of frequency axis
            toplot.mpowr=[];
            toplot.mpowr=mean(S_tot,1); %Average spectrum across all pixels
            toplot.f = [];
            toplot.f = f;
            figure('WindowStyle','docked');
            plot(f, log10(mean(S_tot, 1)));
        end



    end
end





%%



