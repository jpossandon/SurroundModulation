%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E354 surround supression spectral analysis per particpant
% it requires preproceesed files (E354_Sxxxx.set for hamburg and EEG_EXP_SS_xxxx.set)
% whic are obtained with E354_preprocessing and saved in
% <paths.preproc>/Sxxxx
% This script computes the spectral decompoistion for each experimental
% condition, both the psd of the trial average ('PxxAvg') and the average
% of the psd of each trial ('meanPxx'), and respective SNR. The result is
% saved in one file per participant:
% <paths.Analysis>/SSVEPindividual/Sxxxx_SSVEP.mat
% and in file for all participants:
% <paths.Analysis>/SSVEPgroup/allParticipants_SSVEP.mat
% In addition two figure are created for each participant, showing spectra and
% topoplots per conditions, saved in
% <paths.Analysis>/SSVEPindividual/figures/Sxxxx_SSVEP_<'meanPxx' or 'PxxAvg'>
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% where does the data comes from?
% SITE = 'HYDERABAD';
SITE = 'HAMBURGv2';
E354_paths_globalParams      % script with all paths, this file need to be changed accordingly

DO_INDIVIDUAL_PLOTS = 0;
%%
% Getting individual preprocessd data, epoching, perodogram and individual plots
for subj = 1:length(subjIDs)
    if strcmp(SITE,'HYDERABAD') % participant filename
        filename        = ['EEG_EXP_SS_' subjIDs{subj}(2:end)]; 
    elseif strcmp(SITE,'HAMBURG')
         filename        = ['E354_' subjIDs{subj}]; 
    elseif strcmp(SITE,'HAMBURGv2')
         filename        = ['E354v2_' subjIDs{subj}];      
    end
    EEG = pop_loadset('filename',[filename,'.set'],'filepath',fullfile(paths.preproc,subjIDs{subj}));  % loading preproc data

    load(fullfile(paths.Raw,subjIDs{subj},[filename '.mat']))                                       % loading file with experiment information
    triggers = unique([TRIAL.trigger]);     
    triggers = triggers(triggers<200);                                                              % triggers below 200 are non-target trials
    clear Pxx
    snrSize  = 10;      % how manu adyacent frequency bins totake for the SNR calculation
    snrKernel = ones(1,snrSize+3)./snrSize;     % the is an averaging 'kernel' taking the average of snrSize bin adjacent to any sample to compute the SNR with convolution
    snrKernel(snrSize/2+1:snrSize/2+3) = 0;     % the two most adyacent bins are bnot used in the calculation
    
    % loop trough triggers
    for trg = 1:length(triggers)
        thisTrigger     = triggers(trg);
        ixfirsttrig     = find([TRIAL.trigger]==thisTrigger,1);
        thisTriggerstr  = ['S' num2str(thisTrigger)];
        thisEpoch       = pop_epoch( EEG, {  thisTriggerstr }, [0  11], 'newname', [thisTriggerstr 'epochs'], 'epochinfo', 'yes');  % trials are 11 second duration

        if strcmp(filename,'E354v2_S025') & strcmp(thisTriggerstr,'S116') % In this participant, the 'test'/demonstration trial was included in the experimenal recording so here I am removing it
            thisEpoch.data(:,:,1) = [];
            thisEpoch.icaact(:,:,1) = [];
        end
        thisEEGdata     = thisEpoch.data(:,EEG.srate*.25+1:EEG.srate*10.25,:); % but we only take between 250 and 10250 ms after the start (10 second exact to the frequency resolution is 0.1 Hz and integer frequencies are coverd, e.g, 15,30,45)
        % fourier transform of the complete 10 second section with hamming window
        for trl = 1:size(thisEEGdata,3) % for each trial
            [Pxx(:,:,trl),F]         = periodogram(squeeze(thisEEGdata(:,:,trl))',hamming(size(thisEEGdata,2)),size(thisEEGdata,2),EEG.srate);
        end
        % and for the mean of the trials, this is better when the
        % stimulation has always the same phase (starts at the same moment with exactly the same stimulsu for a given trigger)
        [PxxAvg,F]                  = periodogram(mean(thisEEGdata,3)',hamming(size(thisEEGdata,2)),size(thisEEGdata,2),EEG.srate);
        % fill a result structure with all the relevant information
        result(trg).meanPxx         = mean(Pxx(F>1&F<70,:,:),3);  % result taking the mean of spectra of each trials
        result(trg).PxxAvg          = PxxAvg(F>1&F<70,:);         % result spectra of the trials time series mean
        avg_adj_meanPxxbins         = conv2(result(trg).meanPxx,snrKernel','same');
        avg_adj_PxxAvgbins          = conv2(result(trg).PxxAvg,snrKernel','same');
        std_adj_meanPxxbins         = nan(size(avg_adj_meanPxxbins));
        std_adj_PxxAvgbins          = nan(size(avg_adj_PxxAvgbins));
        for smpl = snrSize/2+2:size(result(trg).meanPxx,1)-(snrSize/2+1)
            std_adj_meanPxxbins(smpl,:)         = std(result(trg).meanPxx([smpl-snrSize/2-1:smpl-2,smpl+2:smpl+snrSize/2+1],:));
            std_adj_PxxAvgbins(smpl,:)          = std(result(trg).PxxAvg([smpl-snrSize/2-1:smpl-2,smpl+2:smpl+snrSize/2+1],:));
        end
        result(trg).SNR_meanPxx     = result(trg).meanPxx./avg_adj_meanPxxbins; % respective SNRs (signal divided by the average 0f 10 adyacen frequency bi without taking the two inmediatly adyacent)
        result(trg).SNR_PxxAvg      = result(trg).PxxAvg./avg_adj_PxxAvgbins;   % respective SNRs
        result(trg).Zscore_meanPxx  = (result(trg).meanPxx-avg_adj_meanPxxbins)./std_adj_meanPxxbins; % respective Zscores (signal minus by the average of 10 adyacent frequency bins divided by the std of the 10 adyacent bins,  without taking the two inmediatly adyacent)
        result(trg).Zscore_PxxAvg   = (result(trg).PxxAvg-avg_adj_PxxAvgbins)./std_adj_PxxAvgbins;
        result(trg).F               = F(F>1&F<70);
        result(trg).trigger         = thisTrigger;
        result(trg).grt_contrast    = TRIAL(ixfirsttrig).grt_contrast;
        result(trg).bkg_size        = TRIAL(ixfirsttrig).bkg_size;
        result(trg).bkg_orient      = TRIAL(ixfirsttrig).bkg_orient;
        result(trg).bkg_contrast    = TRIAL(ixfirsttrig).bkg_contrast;
        result(trg).subject         = subjIDs{subj};
    end

    %
    % summary figure per participant, this work only for experiment 1 so
    % far
    if DO_INDIVIDUAL_PLOTS
        close all
        
        % all this is for generating the figure layou
        spHsp       = 1/11;
        spVsp       = 1/5;
        clear subplotPos
        subplotPos  = [.5.*spHsp 3.25.*spVsp spHsp spVsp;...   %50 none none power
                        2.*spHsp 3.8.*spVsp spHsp spVsp;...      %50 parall 2 power
                        3.1.*spHsp 3.8.*spVsp spHsp spVsp;...      %50 parall 4 power
                        4.2.*spHsp 3.8.*spVsp spHsp spVsp;...      %50 parall full power
                          2.*spHsp 1.3.*spVsp spHsp spVsp;...      %50 orto 2 power
                        3.1.*spHsp 1.3.*spVsp spHsp spVsp;...      %50 orto 4 power
                        4.2.*spHsp 1.3.*spVsp spHsp spVsp];      %50 orto full power
       plotOrder = table([repmat({'50'},7,1);repmat({'100'},7,1)],...
                         repmat({'none' 'parallel' 'parallel' 'parallel' 'orthogonal' 'orthogonal' 'orthogonal'}',2,1),...  
                         repmat({'none' '2' '4' 'full' '2' '4' 'full'}',2,1),...
                         'VariableNames',{'grt_contrast','bkg_orient','bkg_size'});  
        subplotPos = [subplotPos;subplotPos];
        subplotPos(8:end,1) = subplotPos(8:end,1)+5.5.*spHsp;
    
        channelsToAvg = {'O1','Oz','O2','PO3','POz','PO4'};           % which channels to average for the spectra plot
        % channelsToAvg = {'Fp1','Fpz','Fp2'};
        chindx = ismember({EEG.chanlocs.labels},channelsToAvg);
        freq        = 15;   % which freuqneyc to plot the spectra
        axLim       = [0 65 0 4];
        axLimSNR    = [0 65 0 25];
        axLimZscore = [0 65 0 25];
        SSVEPtype   = {'PxxAvg','meanPxx'};
    
        for sT = 1%
            % :2
            fh          = figure;
            fh.Units    = 'centimeters';
            figSiz      = [17.6 17.6/2];
            fh.Position = [5,5, figSiz*2];
            for trg = 1:length(triggers)
              
                indxSp = find(ismember(plotOrder.grt_contrast,result(trg).grt_contrast)&ismember(plotOrder.bkg_size,result(trg).bkg_size)&ismember(plotOrder.bkg_orient,result(trg).bkg_orient));
                % power
                axes('Position',subplotPos(indxSp,:))    
                plot(result(trg).F,mean(result(trg).(SSVEPtype{sT})(:,chindx),2));
                box off
                set(gca,'XTick',[])
                axis(axLim),vline([freq freq*2])
                if ismember(result(trg).bkg_size,{'4','full'}),set(gca,'YTick',[]),else,ylabel('Power'),end
                title(sprintf('%s %s %s',result(trg).grt_contrast,result(trg).bkg_size,result(trg).bkg_orient))
                axes('Position',subplotPos(indxSp,:)+[spHsp/3 spVsp/3 -spHsp/4 -spVsp/4])
                topoplot(result(trg).(SSVEPtype{sT})(result(trg).F==freq,:),EEG.chanlocs,'maplimits',[-axLim(4) axLim(4)])
                fixtopotlines
                
                % SNR
                % axes('Position',subplotPos(indxSp,:)-[0 spVsp+.1*spVsp 0 0])
                %  plot(result(trg).F,mean(result(trg).(['SNR_' SSVEPtype{sT}])(:,chindx),2));
                %  axis(axLimSNR),vline([freq freq*2 freq*3 freq*4])
                %   box off
                %  if ismember(result(trg).bkg_size,{'4','full'}),set(gca,'YTick',[]),else,ylabel('SNR'),end
                %  axes('Position',subplotPos(indxSp,:)-[0 spVsp+.1*spVsp 0 0]+[spHsp/3 spVsp/3 -spHsp/4 -spVsp/4])
                % topoplot(result(trg).(['SNR_' SSVEPtype{sT}])(result(trg).F==freq,:),EEG.chanlocs,'maplimits',[-axLimSNR(4) axLimSNR(4)]);
                % fixtopotlines
                % 
                % axes('Position',[5.5.*spHsp,.3*spVsp,.15,1.*spVsp])
                % axis off
                % if strcmp(SSVEPtype{sT},'PxxAvg')
                %     text(-.6,.9,sprintf('PSD of trials mean\nAvg. electrodes:'))
                % else
                %     text(-.6,.9,sprintf('Mean PSD\nAvg. electrodes:\n%s',strjoin(channelsToAvg)))
                % end
                % topo_markCh(EEG.chanlocs,find(chindx))


                % Zscore
                axes('Position',subplotPos(indxSp,:)-[0 spVsp+.1*spVsp 0 0])
                 plot(result(trg).F,mean(result(trg).(['Zscore_' SSVEPtype{sT}])(:,chindx),2));
                 axis(axLimZscore),vline([freq freq*2 freq*3 freq*4])
                  box off
                  set(gca,'XTick',[0:15:60])
                 if ismember(result(trg).bkg_size,{'4','full'}),set(gca,'YTick',[]),else,ylabel('Zscore'),end
                 axes('Position',subplotPos(indxSp,:)-[0 spVsp+.1*spVsp 0 0]+[spHsp/3 spVsp/3 -spHsp/4 -spVsp/4])
                topoplot(result(trg).(['Zscore_' SSVEPtype{sT}])(result(trg).F==freq,:),EEG.chanlocs,'maplimits',[-axLimZscore(4) axLimZscore(4)]);
                fixtopotlines
        
                axes('Position',[5.5.*spHsp,.3*spVsp,.15,1.*spVsp])
                axis off
                if strcmp(SSVEPtype{sT},'PxxAvg')
                    text(-.6,.9,sprintf('PSD of trials mean\nAvg. electrodes:'))
                else
                    text(-.6,.9,sprintf('Mean PSD\nAvg. electrodes:\n%s',strjoin(channelsToAvg)))
                end
                topo_markCh(EEG.chanlocs,find(chindx))
                 % axis(axLim),vline([freq freq*2])
        
                % topoplot(result(trg).PxxAvg(result(trg).F==freq,:),EEG.chanlocs);
                % text(-.5,-.6,sprintf('%d Hz PSD of trials mean',freq))
                % subplot(2,2,2)
                 
            end
            doimage(fh,fullfile(paths.Analysis,'SSVEPindividual','figures'),'pdf',sprintf('%s_SSVEP_%s',subjIDs{subj},SSVEPtype{sT}),'300','painters',figSiz*2,1)
        end
    end
    save(fullfile(paths.Analysis,'SSVEPindividual',sprintf('%s_SSVEP',subjIDs{subj})),'result')
end

%%
% all subjects result in one structure
for subj = 1:length(subjIDs)
    load(fullfile(paths.Analysis,'SSVEPindividual',sprintf('%s_SSVEP',subjIDs{subj})),'result')
    allResult(subj,:) = result;
end
save(fullfile(paths.Analysis,'SSVEPgroup','allParticipants_SSVEP'),'allResult')

chanlocs = EEG.chanlocs;
save(fullfile(paths.Analysis,'E354_chanlocs.mat'),'chanlocs');

  
    