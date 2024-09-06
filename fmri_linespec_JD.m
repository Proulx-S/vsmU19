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
    for run = 1:4
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
        %Stim frequency
        stimT = str2double(extractBefore(vfMRI_tmp.dsgn.label,'s'));
        stimFreq = 1/stimT; % Hz

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
            plot(log10(diag(S).^2));
            xlabel('Mode','Interpreter','latex');
            ylabel('Log10 Eigenvalue $\sigma^2$','Interpreter','latex');

            sig_modes = round(size(data,1)/2); %Start with only half the modes as test.
            % sig_modes = 10;
            % sig_modes = round(size(data,1));
            % sig_modes = 8;
            output.sig_modes = sig_modes;
            Un=single(U(:,1:sig_modes));
            Sn=single(S(1:sig_modes,1:sig_modes));
            Vn=single(V(:,1:sig_modes));

            tmp_mask = vfMRI_tmp.volTs(run).mri.vol2vec;

            % for ii = 1:7
            %     find_mask = find(tmp_mask);
            %     mode1 = NaN(size(tmp_mask));
            %     mode1(find_mask) = U(:,ii);
            %     figure;
            %     subplot(2,1,1);
            %     imagesc(mode1);
            %     xlim([130 240])
            %     ylim([150 260])
            %     daspect([1,1,1]);
            %     colorbar
            %     title(['Mode ',num2str(ii)]);
            %     subplot(2,1,2);
            %     plot(Tvec,Vn(:,ii));
            %     xlabel('Time (s)')
            %     ylabel('Temporal mode value')
            % end


            clearvars U S V
            figure('WindowStyle','docked');
            plot(Tvec,Vn(:,1),'k'); hold on;
            plot(Tvec,Vn(:,2),'b');
            % plot(Tvec,Vn(:,3),'g');
            xline(Stimvec,'Color','k','Alpha',0.2);
            xlabel('Time (s)','Interpreter','latex');
            ylabel('Temporal mode value','Interpreter','latex');
            legend({'Mode 1','Mode 2','Stim Start'})
            %TO DO: PLOT SPATIAL MODES BACK TO MASK

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
            % [slep,~] = dpss(num_frame, p, num_tapers);
            %% Perform Spectral FFT
            % Update half-bandwidth according to the rounded p value
            Delta_f = p * Fs / num_frame;
            addpath(genpath('C:\chronux_2_12'))
            ntapers = [(num_tapers+1)/2,num_tapers];
            % Fs = toplot.rate;
            nfft = num_frame_pad;
            tapers = dpsschk(ntapers,size(Vn,1),Fs);
            %TEMPORARY
            % tapers = tapers/sqrt(Fs);

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
            summaryfig = figure('WindowStyle','docked');
            plot(f, log10(mean(S_tot, 1)),'k');
            xlabel('Frequency (Hz)','Interpreter','latex');
            ylabel('log10(Power)','Interpreter','latex');
            title({'Average Spectrum across all pixels',['Half BW = ',num2str(round(Delta_f,4))]},'Interpreter','latex');


            %% CALCULATE GLOBAL COHERENCE
            % Use same f vector from spectrum calculation
            plot_num_f_points = length(f);
            toplot.coherence = zeros(plot_num_f_points,1);
            tic
            coherence = zeros(length(f),1);
            for i = 1:length(toplot.f)
                m = scores*squeeze(taperedFFT(i,:,:))'; %m is space x tapers
                s = svd(m,0);
                coherence(i) = squeeze(s(1))^2./sum(s.^2); %"Global coherence" in Observed brain dynamics p210
            end
            toplot.coherence = coherence;
            toc;
            yyaxis right
            plot(f,coherence,'b')
            ylim([0 1]);
            ylabel('Spatial Coherence','Interpreter','latex','Color','b')
            set(gca,'YColor',[0 0 1]);

            %% CALCULATE LINE SPECTRA
            tapers_FT = fft(tapers,nfft)/Fs;
            tapers_FT = tapers_FT(1,:);
            t_norm = sum(tapers_FT.^2);

            % Amplitudes at all frequencies f
            A = zeros(length(f),size(taperedFFT,3),'single');
            A = complex(A);
            for k = 1:ntapers(2)
                A = A + tapers_FT(k).*squeeze(taperedFFT(:,k,:));
            end
            A = A./t_norm;
            mu = scores*A';
            toplot.amp = mu;

            % Amplitudes for the stim frequencies and harmonics
            f_stim = stimFreq*[1:1:50];
            f_stim(f_stim > Fs/2) = []; %Stim frequencies and harmonics.
            % Find f_stim values in f vector (and therefore tapered FFT).
            % Just do a loop for now.
            f_stim_loc = zeros(length(f_stim),1);
            for i = 1:length(f_stim)
                f_diff = abs(f - f_stim(i));
                [~,minloc] = min(f_diff);
                f_stim_loc(i) = minloc;
            end

            A = zeros(length(f_stim),size(taperedFFT,3),'single');
            A = complex(A);
            for k = 1:ntapers(2) 
                A = A + tapers_FT(k).*squeeze(taperedFFT(f_stim_loc,k,:));
            end
            A = A./t_norm;
            toplot.amp_fstim = scores*A';

            % CALCULATE F-STATISTICS
            N = length(Tvec);
            p = 0.05; % Significance level for f-test
            p= p/length(Tvec); 
            sig=finv(1-p,2,2*ntapers(2)-2);
            disp(['Calculating F-Statistics at significance level ',num2str(p),' F_sig = ',num2str(sig)])
            % All frequencies.
            Fde = zeros(size(mu),'single');
            for k = 1:ntapers(2)
                z = scores*squeeze(taperedFFT(:,k,:))' - mu*tapers_FT(k); 
                Fde = Fde + conj(z).*z;
            end
            F = (size(taperedFFT,2)-1)*(conj(mu).*mu)*t_norm./Fde; % Space x freq
            % [Fval,A,f,sig] = ftestc(data,params,p,plt); % from Chronux

            % Stim frequencies.
            F_stim = F(:,f_stim_loc);

            % CALCULATE RESIDUAL SPECTRUM
            % Already have F-values and Amplitudes
            toplot.amp=toplot.amp*Fs; %Check this.
            toplot.amp_fstim = toplot.amp_fstim*Fs;

            %From fitlinesc.m
            fmax = findpeaks(F',sig); %This is the chronux findpeaks
            C = size(F,1); %Channels
            freqs=cell(1,C); 
            Amps=cell(1,C);
            data_reconstruct = Un*Sn*Vn';
            datafit = data_reconstruct';
            %%
            clearvars i
            A = toplot.amp';
            for ch=1:C
                fsig=f(fmax(ch).loc);
                freqs{ch}=fsig;
                Amps{ch}=toplot.amp(ch,fmax(ch).loc);
                Nf=length(fsig);
                datafit(:,ch)=exp(i*2*pi*(0:N-1)'*fsig/Fs)*A(fmax(ch).loc,ch)+exp(-i*2*pi*(0:N-1)'*fsig/Fs)*conj(A(fmax(ch).loc,ch));
                %These are the sinewaves. Only subtract if there were significant sine components.
            end
            data_nolines = data_reconstruct' - datafit;


            params.tapers = ntapers;
            params.Fs = Fs;
            params.pad = 2;
            [Stot,f]=mtspectrumc(data_mean',params);
            [Sresid,f]=mtspectrumc(data_nolines,params);

            % hold on
            % plot(f,log10(mean(Stot,2)))
            figure(summaryfig);
            yyaxis left
            hold on
            plot(f,log10(mean(Sresid,2)),'m--');

            %%



        end



    end
end





%%



