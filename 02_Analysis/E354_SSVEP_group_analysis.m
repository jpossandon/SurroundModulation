%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E354 surround supression spectral analysis for all participants
% it requires the spectral result file for all participants generated by
% E354_subjProc.m:
% <paths.Analysis>/SSVEPgroup/allParticipants_SSVEP.mat
%%%%%%%%%%%%%%%%%%%%%%

%%
clear all

% where does the data comes from?
% SITE = 'HYDERABAD';
SITE = 'HAMBURGv2';
E354_paths_globalParams      % script with all paths, this file need to be changed accordingly

load(fullfile(paths.Analysis,'E354_chanlocs.mat'))

RECOMPUTE_SPECTRAL_AVERAGES = 0;
%%

 % spectral average for conditions without backgroud
if RECOMPUTE_SPECTRAL_AVERAGES 
    load(fullfile(paths.Analysis,'SSVEPgroup','allParticipants_SSVEP'),'allResult')
    
    measure                 = 'PxxAvg';
    allNoneNone.psd         = [];
    allNoneNone.SNR         = [];
    allNoneNone.Zscore      = [];
    for subj = 1:size(allResult,1) % allResult(participant,condition)
    
       
        ixcond          = find(ismember({allResult(subj,:).bkg_size},'none') & ismember({allResult(subj,:).bkg_orient},'none')); % which conditions to use (it is redundant to find these indexes per participant and both for size and orient but doing it like this cannot fail)
        thisParticipant = [];       thisParticipantSNR = [];     thisParticipantZscore = [];
        for cc = 1:length(ixcond)
            thisParticipant      = cat(3,thisParticipant,allResult(subj,ixcond(cc)).(measure));
            thisParticipantSNR   = cat(3,thisParticipant,allResult(subj,ixcond(cc)).(['SNR_' measure]));
            thisParticipantZscore   = cat(3,thisParticipant,allResult(subj,ixcond(cc)).(['Zscore_' measure]));
        end
        allNoneNone.psd(:,:,subj)    = mean(thisParticipant,3);
        allNoneNone.SNR(:,:,subj) = mean(thisParticipantSNR,3);
        allNoneNone.Zscore(:,:,subj) = mean(thisParticipantZscore,3);
    end
    FF = allResult(subj,ixcond(cc)).F;
    clear thisParticipant thisParticipantSNR thisParticipantZscore ixcond subj cc
    
    % spectral average per participant for all conditions, diffrent measures
    freqs = [15 30 45 60];
    grpChSel    = find(median(allNoneNone.Zscore(FF==15,:,:),3)>anlys.zscoreThresh);
    for freq = 1:length(freqs)
        subjs = {}; bkgSize = {}; bkgOrient = {};foreContrast = {};  bkgContrast = {};
        psd_indvChSel = []; SNR_indvChSel = []; zscore_indvChSel = [];
        psd_grpChSel = []; SNR_grpChSel = []; zscore_grpChSel = [];
        for subj = 1:size(allResult,1)
            for cond = 1:size(allResult,2)
                thisSt           = allResult(subj,cond);
                subjs            = [subjs;thisSt.subject];
                bkgSize          = [bkgSize;thisSt.bkg_size];
                bkgOrient        = [bkgOrient;thisSt.bkg_orient];
                bkgContrast      = [bkgContrast;thisSt.bkg_contrast];
                foreContrast     = [foreContrast;thisSt.grt_contrast];
                psd_grpChSel     = [psd_grpChSel;mean(thisSt.(measure)(thisSt.F==freqs(freq),grpChSel),2)];
                SNR_grpChSel     = [SNR_grpChSel;mean(thisSt.(['SNR_' measure])(thisSt.F==freqs(freq),grpChSel),2)];
                zscore_grpChSel  = [zscore_grpChSel;mean(thisSt.(['Zscore_' measure])(thisSt.F==freqs(freq),grpChSel),2)];
                indvChSel        = find(allNoneNone.Zscore(FF==15,:,subj)>anlys.zscoreThresh);
                psd_indvChSel    = [psd_indvChSel;mean(thisSt.(measure)(thisSt.F==freqs(freq),indvChSel),2)];
                SNR_indvChSel    = [SNR_indvChSel;mean(thisSt.(['SNR_' measure])(thisSt.F==freqs(freq),indvChSel),2)];
                zscore_indvChSel = [zscore_indvChSel;mean(thisSt.(['Zscore_' measure])(thisSt.F==freqs(freq),indvChSel),2)];
            end
        end
        resTables.([sprintf('table_%d',freqs(freq))]) = table(categorical(subjs),categorical(foreContrast),categorical(bkgSize),categorical(bkgOrient),categorical(bkgContrast),psd_grpChSel,SNR_grpChSel,zscore_grpChSel,psd_indvChSel,SNR_indvChSel,zscore_indvChSel);
        resTables.([sprintf('table_%d',freqs(freq))]).Properties.VariableNames(1:5) = {'subject','contrast','bkgSize','bkgOrient','bkgContrast'};
    end
    clear freqs freq subjs bkgSize bkgOrient bkgContrast foreContrast psd_indvChSel SNR_indvChSel zscore_indvChSel psd_grpChSel SNR_grpChSel zscore_grpChSel cond subj thisSt
    
    save(fullfile(paths.Analysis,'SSVEPgroup','SSVEPtables'),'allNoneNone','resTables')
    
end
%%

% summary figures for SSVEP when there is not background, it creates
% figures for different frequencies (15,30,45,60) and measures (uV,log(uV),SNR,Zscore) with a 
% topoplot for each participant and a figure with  topoplots with the
% average across participants (for the different frequencies and measure)

 E354_summaryfigs_noBkg


 E354_SSVEPresult_scattters_v2


 % E354_SSVEP_psych_corr


      
