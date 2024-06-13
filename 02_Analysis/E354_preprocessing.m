%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E354 surround supression experiment preprocessing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% where does the data comes from?
% SITE = 'HYDERABAD';
SITE = 'HAMBURGv2';
E354_paths_globalParams      % script with all paths and participants, this file need to be changed accordingly



%%
for subj = 1:length(subjIDs)
    mkdir(fullfile(paths.preproc,subjIDs{subj}))                % create participant analysis subdirectory
    
    if strcmp(SITE,'HYDERABAD') % participant filename
        filename        = ['EEG_EXP_SS_' subjIDs{subj}(2:end)]; 
    elseif strcmp(SITE,'HAMBURG')
         filename        = ['E354_' subjIDs{subj}]; 
    elseif strcmp(SITE,'HAMBURGv2')
         filename        = ['E354v2_' subjIDs{subj}]; 
    end
   load(fullfile(paths.Raw,subjIDs{subj},[filename '.mat']))   % open experiment mat file contianing PARAMS and TRIAL info
   
    % create mat file containing experiment and preprocessing info
    preprocDataFile     = fullfile(paths.preproc,subjIDs{subj},[filename '.mat']); 
    preprocSt.subj      = subjIDs{subj};
    preprocSt.created   = datetime('now');
    save(preprocDataFile ,'PARAMS','TRIAL','preprocSt')

    % open the eeg file
    EEG         = pop_fileio(fullfile(paths.Raw,subjIDs{subj},[filename '.vhdr']), 'dataformat','auto'); % opens the EEG file
    % check triggers
    % E354_checkTriggers
    
    % removing file end for participants in which the cap was
    % removed/unpluged before ending the recording
    
    EEG         = pop_chanedit(EEG, 'lookup',paths.standardBesa);                            % get channel location according to standard naming nomenclature
    EEG         = pop_resample( EEG, preproc.rsfreq );                                                  % downsample
    EEG         = eeg_checkset( EEG );                                                        
    EEG         = pop_zapline_plus(EEG, 'noisefreqs','line','coarseFreqDetectPowerDiff',...
                    5,'chunkLength',0,'adaptiveNremove',1,'fixedNremove',1,'plotResults',0);        % remove line noise
    EEG         = pop_eegfiltnew(EEG, 'locutoff',0.1,'plotfreqz',0);                                % remove DC (~data mean) component
    EEG         = eeg_checkset( EEG );
    
    channelsTEMP    = {EEG.chanlocs.labels};
    origChanlocs    = EEG.chanlocs;               % this is used later for interpolation


     [~,dead_segments_perc,noise_segments_perc] = thresholdData_range_std(EEG.data,EEG.srate,preproc.badCh);
    removeChannel   = find(dead_segments_perc>preproc.badCh.removeAbove | noise_segments_perc>preproc.badCh.removeAbove);
    if ~isempty(removeChannel) 
        EEG             = pop_select( EEG, 'rmchannel',channelsTEMP(removeChannel));
    end
    leftchannels    = {EEG.chanlocs.labels};
    [preprocSt.preproc.badChan, preprocSt.preproc.badChanindex] = setdiff(channelsTEMP, leftchannels); % whih channels were removed
    clear channelsTEMP removeChannel leftchannels

    % cleannig data for ICA
    [preprocSt.bad_segments_forICA]  = thresholdData_range_std(EEG.data,EEG.srate,preproc.badChforICA);
    EEGforICA          = pop_select( EEG, 'rmpoint',preprocSt.bad_segments_forICA);
    
    EEG             = pop_reref(EEG, []);               % average reference (both for ICA and nonICA EEG structure)
    EEGforICA       = pop_reref(EEGforICA, []);
    EEGICA          = pop_runica(EEGforICA, 'icatype', 'runica', 'extended',1,'interrupt','on');     % ICA
    
    EEG             = pop_editset(EEG,'icaweights',EEGICA.icaweights,'icasphere',EEGICA.icasphere);  % add the ICA matrix to the originald ataset without removed data
    EEG             = eeg_checkset( EEG );
    clear EEGforICA EEGICA
    preprocSt.icaweights = EEG.icaweights;          % save orignal ICA matric in the preprocessing file
    preprocSt.icasphere  = EEG.icasphere;
  

    % removal or muscle, heart and eye component
    iclabeltresh = .9;
    EEG = iclabel(EEG); 
    EEG = pop_icflag(EEG, [NaN NaN;iclabeltresh 1;iclabeltresh 1;iclabeltresh 1;NaN NaN;iclabeltresh 1;NaN NaN]); % The 7 categories are (in order) Brain, Muscle,Eye, Heart, Line Noise, Channel Noise, Other.

    % figure for checking ICA
    pop_topoplot(EEG, 0, [1:32] ,'',[6 6] ,0,'electrodes','off');
    annotation('textbox',[.1 .015 .5 .075],'String',sprintf('Muscle components: %s\nEye components: %s\nECG components: %s\nChnoise components: %s',...
        num2str(sum(EEG.etc.ic_classification.ICLabel.classifications(:,2) >= iclabeltresh)),num2str(sum(EEG.etc.ic_classification.ICLabel.classifications(:,3) >= iclabeltresh)),...
        num2str(sum(EEG.etc.ic_classification.ICLabel.classifications(:,4) >= iclabeltresh)),num2str(sum(EEG.etc.ic_classification.ICLabel.classifications(:,6) >= iclabeltresh))))
    doimage(gcf,fullfile(paths.preproc,subjIDs{subj}),'png',[filename '_ICA'],'100','painters',[],0)
    
    preprocSt.comprejectIClabel = EEG.reject.gcompreject;
    save(preprocDataFile,'preprocSt')
    EEG = pop_subcomp( EEG, [], 0);  % actual ICA component removal
    EEG = eeg_checkset( EEG );
    
    % interpolate rejetedchannels
    if ~isempty(preprocSt.preproc.badChanindex)
        EEG = pop_interp(EEG, origChanlocs, 'spherical');
    end
    save(preprocDataFile,'preprocSt')
    EEG = pop_saveset( EEG, 'filename',[filename,'.set'],'filepath',fullfile(paths.preproc,subjIDs{subj}));
end


