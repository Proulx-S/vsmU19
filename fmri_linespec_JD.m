% fmri_linespec_JD.m

% Analysis code for line spectral analysis on human (single vessel)
% stimulated fmri datasets. 

%% Load data
clear; clc; close all;
cd('Y:\DataAnalysis\MRI\Human240904\13685568');
dataDir = fullfile(pwd,'datav2');
dataFile = 'vsmDrivenP1.mat';
% dataFile = 'vsmDrivenP2.mat';
disp(['processing ',dataFile]);
subj = extractBetween(dataFile,'Driven','.mat');
subj = subj{1};
load(fullfile(dataDir,dataFile))
%% Analyze
trials = fields(vfMRI);
% for trial = 9:length(trials)
for trial = 4:4
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


    % for run = 1:length(vfMRI_tmp.volTs)
    for run = 3:3
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
        try
            stim_omit = vfMRI_tmp.dsgn.nullTrial;
            Stimvec(stim_omit) = []; %There are omitted stims outside of onsetlist size
        catch ME
        end
        %Stim frequency
        stimT = str2double(extractBefore(vfMRI_tmp.dsgn.label,'s'));
        stimFreq = 1/stimT; % Hz

        toplot.Fs = Fs;
        toplot.Tvec = Tvec;
        toplot.Stimvec = Stimvec;
        % h.mainfig = figure('WindowStyle','docked');
        % h.mainfig = figure();
        % h.tabgroup = uitabgroup;
        % h.tab(1) = uitab(h.tabgroup,'Title','Frame Time Diff');
        % h.axes(1) = axes('Parent',h.tab(1));
        % plot(diff(Tvec)); ylim([0 2]);
        % xlabel('Time (s)','Interpreter','latex');
        % ylabel('Difference between acquisition frames','Interpreter','latex');

        %SUBTRACT MEAN FROM TIMESERIES
        [data_mean] = fun_MeanSubtract(data);

        if SVD_Q == 1
            [U,S,V]=svd(data_mean,0);
            % figure('WindowStyle','docked');
            h.mainfig = figure();
            h.tabgroup = uitabgroup;
            h.tab(1) = uitab(h.tabgroup,'Title','S-T SVD Sing. Values');
            h.axes(1) = axes('Parent',h.tab(1));
            plot(log10(diag(S).^2));
            xlabel('Mode','Interpreter','latex');
            ylabel('Log10 Eigenvalue $\sigma^2$','Interpreter','latex');
            title('Space-Time SVD Eigenspectrum','Interpreter','latex');

            sig_modes = round(size(data,1)/2); %Start with only half the modes as test.
            % sig_modes = 10;
            % sig_modes = round(size(data,1));
            % sig_modes = 8;
            toplot.sig_modes = sig_modes;
            Un=single(U(:,1:sig_modes));
            Sn=single(S(1:sig_modes,1:sig_modes));
            Vn=single(V(:,1:sig_modes));

            tmp_mask = vfMRI_tmp.volTs(run).mri.vol2vec;

            clearvars S V
            % figure('WindowStyle','docked');
            h.tab(2) = uitab(h.tabgroup,'Title','Temp.Modes');
            h.axes(2) = axes('Parent',h.tab(2))
            plot(Tvec,Vn(:,1),'k'); hold on;
            plot(Tvec,Vn(:,2),'b');
            % plot(Tvec,Vn(:,3),'g');
            try %RESTING STATE DOESN'T HAVE ANY STIM VALUES
            xline(Stimvec,'Color','k','Alpha',0.2);
            catch ME
            end
            xlabel('Time (s)','Interpreter','latex');
            ylabel('Temporal mode value','Interpreter','latex');
            legend({'Mode 1','Mode 2','Stim Start'})
            %PLOT SPATIAL MODES BACK TO MASK
            for ii = 1:7
                find_mask = find(tmp_mask);
                mode1 = NaN(size(tmp_mask));
                mode1(find_mask) = U(:,ii);
                % figure;
                h.tab(ii+2) = uitab(h.tabgroup,'Title',sprintf('M%.0f',ii));
                h.axes(ii+2) = axes('Parent',h.tab(ii+2));
                subplot(2,1,1);
                imagesc(mode1);
                xlim([130 240])
                ylim([150 260])
                daspect([1,1,1]);
                colorbar
                title(['Mode ',num2str(ii)]);
                subplot(2,1,2);
                plot(Tvec,Vn(:,ii));
                xlabel('Time (s)')
                ylabel('Temporal mode value')
            end

            %CALCULATE AVEAGE SPECTRUM
            num_pixel = size(Un,1);
            num_frame= size(Vn,1);
            % ############################# % Choose bandwidth. More
            % tapers->better for line-spectrum regression
            % Might have to do Delta_f = 0.025 to avoid BW crossover at period = 20s
            
            % if contains(trialName,'Rest')
            %     Delta_f = 0.025;
            % else
                Delta_f = 0.025; %Half-Bandwidth. 0.02-> ~11 tapers | 0.03 -> ~17 tapers | 0.04-> ~22 tapers
                % Delta_f = 0.01; %For 50s period only, for summary figure (extracted amplitude)
            % end
            toplot.Delta_f = Delta_f;
            % ############################# %
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
            toplot.avgpowr=[];
            toplot.avgpowr=mean(S_tot,1); %Average spectrum across all pixels
            toplot.f = [];
            toplot.f = f;
            %summaryfig = figure('WindowStyle','docked');
            h.tab(10) = uitab(h.tabgroup,'Title','AvgSpec');
            h.axes(10) = axes('Parent',h.tab(10));
            plot(f, log10(mean(S_tot, 1)),'k');
            xlabel('Frequency (Hz)','Interpreter','latex');
            ylabel('log10(Power)','Interpreter','latex');
            title({'Average and Residual Spectrum across all pixels',['Half BW = ',num2str(round(Delta_f,4))]},'Interpreter','latex');
            


            %% CALCULATE GLOBAL COHERENCE
            % Use same f vector from spectrum calculation
            plot_num_f_points = length(f);
            tic
            coherence = zeros(length(f),1);
            for i = 1:length(toplot.f)
                m = scores*squeeze(taperedFFT(i,:,:))'; %m is space x tapers
                s = svd(m,0);
                coherence(i) = squeeze(s(1))^2./sum(s.^2); %"Global coherence" in Observed brain dynamics p210
            end
            toplot.coherence = coherence'; %Spatial coherence
            toc;
            % yyaxis right
            % plot(f,coherence,'b')
            % ylim([0 1]);
            % ylabel('Spatial Coherence','Interpreter','latex','Color','b')
            % set(gca,'YColor',[0 0 1]);

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
            toplot.f_stim = f_stim;
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
            
            % Correct for normalization (double check this)
            toplot.amp=toplot.amp*Fs;
            toplot.amp_fstim = toplot.amp_fstim*Fs;

            % SAVE EXT AMPs and AVG POWER VALUES IN TOPLOT STRUCTURE
            meanpowr_vals = toplot.avgpowr(f_stim_loc);
            toplot.meanpwr_vals = meanpowr_vals; %Average spectrum power at stim freq and harmonics.
            Extpowr_vals = toplot.amp_fstim;
            Extpowr_vals = abs(Extpowr_vals).^2;
            Extpowr_vals = Extpowr_vals./(2*Delta_f); %Divide by full bandwidth to get line height
            toplot.Avg_Extpowr_stimfreq = mean(Extpowr_vals,1); %These are all amplitudes, irrespective of F-test result.
            % Same for all frequencies
            Extpowr_vals = abs(toplot.amp).^2;
            Extpowr_vals = Extpowr_vals./(2*Delta_f); %Divide by full bandwidth to get line height
            toplot.Avg_Extpowr_allfreq = mean(Extpowr_vals,1); %These are all amplitudes, irrespective of F-test result.

            % CALCULATE F-STATISTICS
            N = length(Tvec);
            p = 0.05; % Significance level for f-test
            p= p/length(Tvec); 
            sig=finv(1-p,2,2*ntapers(2)-2);
            toplot.F_stat_sig = sig;
            disp(['Calculating F-Statistics at significance level ',num2str(p),' F_sig = ',num2str(sig)])
            % All frequencies.
            Fde = zeros(size(mu),'single');
            for k = 1:ntapers(2)
                z = scores*squeeze(taperedFFT(:,k,:))' - mu*tapers_FT(k); 
                Fde = Fde + conj(z).*z;
            end
            F = (size(taperedFFT,2)-1)*(conj(mu).*mu)*t_norm./Fde; % Space x freq. Plot this in multi-tab figure!
            toplot.F_stat = F;
            % [Fval,A,f,sig] = ftestc(data,params,p,plt); % from Chronux

            % Stim frequencies.
            F_stim = F(:,f_stim_loc);

            % CALCULATE RESIDUAL SPECTRUM
            % Already have F-values and Amplitudes
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
                freqs{ch}=fsig; %List of significant frequencies
                Amps{ch}=toplot.amp(ch,fmax(ch).loc); %Pick out significant amp from list of all amps
                Nf=length(fsig);
                datafit(:,ch)=exp(i*2*pi*(0:N-1)'*fsig/Fs)*A(fmax(ch).loc,ch)+exp(-i*2*pi*(0:N-1)'*fsig/Fs)*conj(A(fmax(ch).loc,ch));
                %These are the sinewaves. Only subtract if there were significant sine components.
            end
            data_nolines = data_reconstruct' - datafit;

            %Save significant extracted amplitudes.
            %Iterate over freqs variable - save significantly extracted
            %amplitudes at each stim frequency. Just do a loop for now.
            sig_Amps = NaN(size(toplot.amp,1),length(f_stim));
            Stim_lowBW = f_stim - Delta_f;
            Stim_highBW = f_stim + Delta_f;
            for ii = 1:length(freqs)
                if ~isempty(freqs{ii}) %If there were significant line components detected for this pixel.
                    freqs_tmp = freqs{ii};
                    for jj = 1:length(freqs_tmp) %Iterate over frequency values
                        freq_val_tmp = freqs_tmp(jj);
                        Amp_tmp = Amps{ii};
                        Amp_tmp = Amp_tmp(jj);
                        check_low = freq_val_tmp > Stim_lowBW;
                        check_high = freq_val_tmp < Stim_highBW;
                        isfreq = find(and(check_low,check_high));
                        if ~isempty(isfreq)
                            sig_Amps(ii,isfreq) = (abs(Amp_tmp).^2)./(2*Delta_f);
                        end
                    end
                end
            end
            toplot.sig_Amps = sig_Amps; %Average over ~NaN values here at each frequency to plot summary figure.
            %Could also count how many pixels have significant line components -> some area-locking metric.

            % params.tapers = ntapers;
            params.tapers = [0.025 toplot.Tvec(end) 1];
            params.Fs = Fs;
            params.pad = 2;
            [Stot,f]=mtspectrumc(data_reconstruct',params);
            % [Stot,f]=mtspectrumc(data_mean',params);
            [Sresid,f]=mtspectrumc(data_nolines,params);
            toplot.resid_f = f;
            toplot.avgResid = mean(Sresid,2);
            toplot.Sresid = Sresid;
            toplot.Stot = Stot;

            % hold on
            % plot(f,log10(mean(Stot,2)))
            % figure(summaryfig);
            % yyaxis left
            hold on
            plot(f,log10(mean(Sresid,2)),'m--');

            %Plot f-statistic summary
            F95 = prctile(toplot.F_stat,95,1);
            F5 = prctile(toplot.F_stat,5,1);
            Favg = mean(toplot.F_stat,1);
            h.tab(11) = uitab(h.tabgroup,'Title','F-Stat');
            h.axes(11) = axes('Parent',h.tab(11));
            hp = plot(toplot.f, Favg); hp.Color = [0 0 0 1]; hold on;
            hp2 = plot(toplot.f,F95); hp2.Color = [0 0 0 0.2];
            hp3 = plot(toplot.f,F5); hp3.Color = [0 0 0 0.2];
            yline(sig);
            ylim([0 20])
            xlabel('Frequency (Hz)','Interpreter','latex');
            ylabel('F-statistic value','Interpreter','latex');
            title({'F-statistic', 'Average over pixels, 95th 5th percentile'},'Interpreter','latex');

            %%

        end

        % Save this trial's results. Just do this by trial for now.
        cd('Y:\DataAnalysis\MRI\Human240904\13685568\results_025HzHalfBW_v2');
        figname = [subj,'_',trialName,'_run',num2str(run),'.fig'];
        savefig(gcf,figname);
        close all
        %Save toplot .mat file for summary figures.
        matname = [subj,'_',trialName,'_run',num2str(run),'_toplot.mat'];
        save(matname,'toplot');
    end
end





%%



