
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Background
% Contrast: 0 (gray)/100%
% Size: Small/Full
% Orientation: parallel/orthogonal
%
% Gratings
% Contrast: 25/50/100

%PsychDefaultSetup(0);
%PsychDebugWindowConfiguration(0)  
beep OFF
clear all
clc

KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');

EXP_LOCATION = 3; % 1 - LVPEI 2 - Hamburg 3 - JPO laptop

if EXP_LOCATION == 1
    resultPATH = 'C:\Users\BPN\Documents\EEG_Experiments\Surround_supression\RawData\psychometric'; % LVPEI
elseif EXP_LOCATION == 2
    resultPATH = 'C:\Users\bpn\Documents\Jose\Surround_supression\RawData\psychometric'; % Babylab EEG
elseif EXP_LOCATION == 3
    resultPATH = '/Users/jossando/trabajo/SSVEP_surroundsuppression/06_RawData/psychometric'; % jose MAC
end

if ~isdir(resultPATH)
    mkdir(resultPATH)
end
testType = questdlg(sprintf('Select task type:\n''Low vision'' and ''Test Trials'' start at high contrast. ''Test trials'' does not save an ouput. Choose ''Normal Vision'' for control participants without filter and ''Low Vision'' for control participants using blurring filters'), ...
	'Surround supression 2AFC', ...
	'Low_Vision','Normal_Vision','Test_Trials','Normal_Vision');

if strcmp(testType,'Test_Trials')
    answer      = {'ERASEME';'jose'};
    fileMAT     = 'ERASEME.mat';
else
    answer=inputdlg({'Participant ID:','Experimenter:'},'SS 2AFC');
    fileMAT             = sprintf('SSpsych_v2_S%03d_%s.mat',str2num(answer{1}),testType);
    if exist(fullfile(resultPATH,fileMAT))
        selection = questdlg(sprintf('%s already exists. Do you want to:',fileMAT),'SS 2AFC',...
                   'Cancel','Overwrite','Cancel');
        switch selection
           case 'Cancel'
               error(sprintf('Experiment cancelled by experimenter, %s already exists',fileMAT))
           case 'Overwrite'
               selection2 = questdlg(sprintf('Are you sure you want to overwrite %s?',fileMAT),'SS 2AFC',...
                        'Cancel','Overwrite','Cancel');
                 switch selection2
                     case 'Cancel'
                         error(sprintf('Experiment cancelled by experimenter, %s already exists',fileMAT))
                     case 'Overwrite'
                         sprintf('%s will be overwritten. Previous file is save in ~/overwritten as CS_%sbackup.mat',fileMAT,answer{1})
                         mkdir(fullfile(resultPATH,'overwritten'))
                         copyfile(fullfile(resultPATH,fileMAT),fullfile(resultPATH,'overwritten',['CS_' answer{1} 'backup.mat']))
                 end
        end
    end
end

timelist = {'.200' '.400' '.800' '1.6'};
stimdur = listdlg('PromptString','Select a stimulus duration (seconds):',...
                      'SelectionMode','single',...
                      'ListString',timelist);
stimdur = str2double(timelist(stimdur));
if isempty(stimdur)
    error('You must select a stimulus duration')
end
PARAMS.subjectID    = answer{1};
PARAMS.experimenter = answer{2};
PARAMS.date         = datestr(now);
clear answer selection selection2 selection3 
Screen('Preference', 'SkipSyncTests', 1);
Screen('Preference', 'ConserveVRAM',  4096)
ListenChar(2)
 % Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;


screens         = Screen('Screens');
screenNumber    = max(screens);
bkg_gray        = 127;

% open full window
[win, windowRect] = PsychImaging('OpenWindow', screenNumber,bkg_gray);
[PARAMS.setup.scr_wdth, PARAMS.setup.scr_hgt] = Screen('WindowSize', win);
PARAMS.setup.screen_ratio = PARAMS.setup.scr_wdth/PARAMS.setup.scr_hgt;
% Make sure the GLSL shading language is supported:
AssertGLSL;
Screen('TextSize', win, 30);
HideCursor(win);
% REFRESH PARAMETERS

PARAMS.setup.ifi                = Screen('GetFlipInterval', win);               % time that the screen needs for the flip command                                  % get the size of the screen in cm


Screen('FillRect',win,127); 
Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);
Screen('Flip', win);

PARAMS.setup.distance_to_screen = 60; % in cms
PARAMS.setup.screen_width       = 53; % in cms (viewPixx EEG 24-inch display size (diagonal),53.5 cm recchek)
PARAMS.setup.pixels_per_degree  = PARAMS.setup.scr_wdth/(2*atan(PARAMS.setup.screen_width/2/PARAMS.setup.distance_to_screen)*180/pi); % TODO acutal calcualtion

% GRATING AND BACKGROUND PARAMETERS
PARAMS.spat_freq_deg            = 1; % in cycles per degree (Sourav, 2018 used 2 cycles per degree)
PARAMS.spat_freq_pix            = PARAMS.spat_freq_deg./PARAMS.setup.pixels_per_degree; % in cycles per pixel
PARAMS.centerdotsize            = PARAMS.setup.pixels_per_degree; %in pixels
PARAMS.ovalPenWidth             = 2;
PARAMS.ovalColor                = [50 50 50];

% GRATING PARAMETERS
PARAMS.GRT.radius               = 2/2*PARAMS.setup.pixels_per_degree; % diameter 2.5 visual degrees
PARAMS.GRT.width                = round(PARAMS.GRT.radius*2); 
PARAMS.GRT.height               = round(PARAMS.GRT.radius*2);
PARAMS.GRT.angle                = 0;
%PARAMS.GRT.contrast    = [0.25 0.5]; % 0.5 is for fullcontrast. From DriftDemo: Amplitude of the grating in units of absolute display intensity range: A setting of 0.5 means that the grating will extend over a range from -0.5 up to 0.5, i.e., it will cover a total range of 1.0 == 100% of the total displayable range. As we select a background color and offset for the grating of 0.5 (== 50% nominal intensity == a nice neutral gray), this will extend the sinewaves values from 0 = total black in the minima of the sine wave up to 1 = maximum white in the maxima. Amplitudes of more than 0.5 don't make sense, as parts of the grating would lie outside the displayable range for your computers displays:
PARAMS.GRT.eccenticity          = 4; % in visual degrees
PARAMS.GRT.angularpos           = [155 25 225 315]; % in degress for TL TR BL BR  

[X,Y] = pol2cart(PARAMS.GRT.angularpos*2*pi/360,PARAMS.GRT.eccenticity*PARAMS.setup.pixels_per_degree);
PARAMS.GRT.positions            = round([X;-Y]);% x,y pixel posicion center of grating from screen center for TL TR BL BR, Y is negative because it would be used as upward form center and y pizels start from zero at top of the screen
PARAMS.GRT.phase                = 0;
clear X Y

% TRIAL SPECIFICATION

PARAMS.TRIAL.trial_per_block     = 25;%25;
PARAMS.TRIAL.trial_per_condition = 50; %50 makeit a multiple of trial_per_block
%PARAMS.TRIAL.trial_duration     = 11; % in seconds
PARAMS.TRIAL.stim_duration      = stimdur; % in seconds
%PARAMS.TRIAL.colortarget_trials = 16; % as a percentage of trials
%PARAMS.TRIAL.target_duration    = .5; % in seconds
PARAMS.TRIAL.grt_positions      = {'peripheral'}; % peripheral,(foveal?)
%PARAMS.TRIAL.grt_contrasts      = {'50' '100'}; % 50,100
PARAMS.TRIAL.bkg_contrasts      = {'0' '10' '30'};%
PARAMS.TRIAL.bkg_sizes          = {'.25','1'};% size of bakground. full is the complete screen, otherwise number n specify a radial background thath extend (n-1)*GRT.radius beyond the grating outermos edge 
PARAMS.TRIAL.bkg_orients        = {'parallel','orthogonal'}; % parallel or orthogonal
if strcmp(testType,'Normal vision')
    PARAMS.TRIAL.PsychCurve.tGuess      = -1;   %low threshold prior 10% contrat
else
     PARAMS.TRIAL.PsychCurve.tGuess      = -.1; %high threhsold prior 79% contrast
end
if strcmp(testType,'Test Trials')
    PARAMS.TRIAL.PsychCurve.tGuessSD    = 1;
else
    PARAMS.TRIAL.PsychCurve.tGuessSD    = 3;
end
% for test trials
% PARAMS.TRIAL.PsychCurve.tGuess      = -.5; 
% PARAMS.TRIAL.PsychCurve.tGuessSD    = 2;
PARAMS.TRIAL.PsychCurve.pThreshold  = 0.82;
PARAMS.TRIAL.PsychCurve.beta        = 3;
PARAMS.TRIAL.PsychCurve.delta       = 0.02;
PARAMS.TRIAL.PsychCurve.gamma       = 0.5;

% The q95Limit is the width of the posterior distribution used for stopping
% the testing for a condition. When the width of the posterior is below the
% limit for at least 3 trials, remaining trials for the condition are
% elimined
PARAMS.TRIAL.PsychCurve.q95Limit    = .05;  
% PARAMS.TRIAL.PsychCurve.range       = 5;

% BACKGROUND PARAMETERS
PARAMS.BKG.angle                = [0 90];
%PARAMS.BKG.contrast    = [0 0.5]; 
PARAMS.BKG.phase                = 0;
% PARAMS.BKG.widths               = [];
% PARAMS.BKG.heights              = [];
for bi = 1:length(PARAMS.TRIAL.bkg_sizes)
        PARAMS.BKG.radius(bi)  = PARAMS.GRT.radius+str2num(PARAMS.TRIAL.bkg_sizes{bi})*PARAMS.setup.pixels_per_degree; % diameter 2.5 visual degrees

end
PARAMS.BKG.widths                = round(PARAMS.BKG.radius*2); 
PARAMS.BKG.heights               = round(PARAMS.BKG.radius*2);
PARAMS.BKG.angle                = 0;
%PARAMS.GRT.contrast    = [0.25 0.5]; % 0.5 is for fullcontrast. From DriftDemo: Amplitude of the grating in units of absolute display intensity range: A setting of 0.5 means that the grating will extend over a range from -0.5 up to 0.5, i.e., it will cover a total range of 1.0 == 100% of the total displayable range. As we select a background color and offset for the grating of 0.5 (== 50% nominal intensity == a nice neutral gray), this will extend the sinewaves values from 0 = total black in the minima of the sine wave up to 1 = maximum white in the maxima. Amplitudes of more than 0.5 don't make sense, as parts of the grating would lie outside the displayable range for your computers displays:
PARAMS.BKG.eccenticity          = 4; % in visual degrees
PARAMS.BKG.angularpos           = [155 25 225 315]; % in degress for TL TR BL BR  

[X,Y] = pol2cart(PARAMS.BKG.angularpos*2*pi/360,PARAMS.BKG.eccenticity*PARAMS.setup.pixels_per_degree);
PARAMS.BKG.positions            = round([X;-Y]);% x,y pixel posicion center of grating from screen center for TL TR BL BR, Y is negative because it would be used as upward form center and y pizels start from zero at top of the screen
PARAMS.BKG.phase                = 0;
clear X Y

                        
% possible conditions
% all combination of factors
pcombs              = cell(1,4);
[pcombs{:}]         = ndgridnotfull(PARAMS.TRIAL.grt_positions,PARAMS.TRIAL.bkg_contrasts,PARAMS.TRIAL.bkg_sizes,PARAMS.TRIAL.bkg_orients);
postrials           = cat(2,pcombs{1}(:),pcombs{2}(:),pcombs{3}(:),pcombs{4}(:));
postrials(ismember(postrials(:,2),'0'),[3 4]) = {'none'}; % change bakground factors that are not possible for contrast 0 to 'none'
postrials           = table2cell(unique(cell2table(postrials)));  % remove trials wth the same factors (becaue of liena bove)

if ~exist('CONDpsych','var') && ~exist('TRIAL','var')
    
    CONDpsych          = struct('grt_position',postrials(:,1),...
        'bkg_contrast',postrials(:,2), 'bkg_size',postrials(:,3), 'bkg_orient',postrials(:,4),'RT',[],'q95',[]);
    for qi = 1:length(CONDpsych)
        % Intensity value from Quest should go from 0 to 1 (~ 0 to 100% contrast)
        CONDpsych(qi).q = QuestCreate(PARAMS.TRIAL.PsychCurve.tGuess,PARAMS.TRIAL.PsychCurve.tGuessSD,...
            PARAMS.TRIAL.PsychCurve.pThreshold,PARAMS.TRIAL.PsychCurve.beta,...
            PARAMS.TRIAL.PsychCurve.delta,PARAMS.TRIAL.PsychCurve.gamma); 
    end
    %structure used during the experimetns and saved

    TRIAL               = struct('grt_position',[],...
        'bkg_contrast',[], 'bkg_size',[], 'bkg_orient',[],... 
        'grt_contrast',[],'grt_position_present',[],'correct',[],'RT',[]);
   trialcond = repmat(1:length(CONDpsych),1,PARAMS.TRIAL.trial_per_condition);
   trialcond = trialcond(randperm(length(CONDpsych)*PARAMS.TRIAL.trial_per_condition,length(CONDpsych)*PARAMS.TRIAL.trial_per_condition));
else
    % check that conditions in CONDpsych are consistent with current
    % experiemtn configuration (same numebr of conidition only)
    if length(CONDpsych)~=size(postrials,1)
       error('Previous psychometric assesment have a different number of conditions as speficied currently in this experiment') 
    end
end
clear pcombs ixtarget targettimes tgttrials
    

[gratingid, gratingrect]        = CreateProceduralSineGrating(win, PARAMS.GRT.width, PARAMS.GRT.height, [0.5 0.5 0.5 0.0], PARAMS.GRT.radius, []);

% LOCATIONS AND PHASE OF GRATING IN SCREEN
% location is given by rect specification [xini,yini,xend,yend] in pixels
% from PARAM.GRT eccentricity and anglarpos with respect to center of
% screen
PARAMS.GRT.Rects.TL             = OffsetRect(CenterRectOnPoint(gratingrect,PARAMS.GRT.positions(1,1),PARAMS.GRT.positions(2,1)),PARAMS.setup.scr_wdth/2,PARAMS.setup.scr_hgt/2);
PARAMS.GRT.Rects.TR             = OffsetRect(CenterRectOnPoint(gratingrect,PARAMS.GRT.positions(1,2),PARAMS.GRT.positions(2,2)),PARAMS.setup.scr_wdth/2,PARAMS.setup.scr_hgt/2);
PARAMS.GRT.Rects.BL             = OffsetRect(CenterRectOnPoint(gratingrect,PARAMS.GRT.positions(1,3),PARAMS.GRT.positions(2,3)),PARAMS.setup.scr_wdth/2,PARAMS.setup.scr_hgt/2);
PARAMS.GRT.Rects.BR             = OffsetRect(CenterRectOnPoint(gratingrect,PARAMS.GRT.positions(1,4),PARAMS.GRT.positions(2,4)),PARAMS.setup.scr_wdth/2,PARAMS.setup.scr_hgt/2);
% phase corrections so PARMAS.GRT/BKG.phase are aligned
PARAMS.GRT.phase_cor.TL         = rem(PARAMS.GRT.Rects.TL(1)/(1/PARAMS.spat_freq_pix),floor(PARAMS.GRT.Rects.TL(1)/(1/PARAMS.spat_freq_pix)))*360+PARAMS.GRT.phase;
PARAMS.GRT.phase_cor.TR         = rem(PARAMS.GRT.Rects.TR(1)/(1/PARAMS.spat_freq_pix),floor(PARAMS.GRT.Rects.TR(1)/(1/PARAMS.spat_freq_pix)))*360+PARAMS.GRT.phase;
PARAMS.GRT.phase_cor.BL         = rem(PARAMS.GRT.Rects.BL(1)/(1/PARAMS.spat_freq_pix),floor(PARAMS.GRT.Rects.BL(1)/(1/PARAMS.spat_freq_pix)))*360+PARAMS.GRT.phase;
PARAMS.GRT.phase_cor.BR         = rem(PARAMS.GRT.Rects.BR(1)/(1/PARAMS.spat_freq_pix),floor(PARAMS.GRT.Rects.BR(1)/(1/PARAMS.spat_freq_pix)))*360+PARAMS.GRT.phase;

for bi = 1:length(PARAMS.TRIAL.bkg_sizes)

    [gratingidBKG(bi), gratingrectBKG]        = CreateProceduralSineGrating(win, PARAMS.BKG.widths(bi), PARAMS.BKG.heights(bi), [0.5 0.5 0.5 0.0], PARAMS.BKG.radius(bi), []);


    PARAMS.BKG.Rects(bi).TL             = OffsetRect(CenterRectOnPoint(gratingrectBKG,PARAMS.BKG.positions(1,1),PARAMS.BKG.positions(2,1)),PARAMS.setup.scr_wdth/2,PARAMS.setup.scr_hgt/2);
    PARAMS.BKG.Rects(bi).TR             = OffsetRect(CenterRectOnPoint(gratingrectBKG,PARAMS.BKG.positions(1,2),PARAMS.BKG.positions(2,2)),PARAMS.setup.scr_wdth/2,PARAMS.setup.scr_hgt/2);
    PARAMS.BKG.Rects(bi).BL             = OffsetRect(CenterRectOnPoint(gratingrectBKG,PARAMS.BKG.positions(1,3),PARAMS.BKG.positions(2,3)),PARAMS.setup.scr_wdth/2,PARAMS.setup.scr_hgt/2);
    PARAMS.BKG.Rects(bi).BR             = OffsetRect(CenterRectOnPoint(gratingrectBKG,PARAMS.BKG.positions(1,4),PARAMS.BKG.positions(2,4)),PARAMS.setup.scr_wdth/2,PARAMS.setup.scr_hgt/2);
    % phase corrections so PARMAS.BKG/BKG.phase are aligned
    PARAMS.BKG.phase_cor(bi).TL         = rem(PARAMS.BKG.Rects(bi).TL(1)/(1/PARAMS.spat_freq_pix),floor(PARAMS.BKG.Rects(bi).TL(1)/(1/PARAMS.spat_freq_pix)))*360+PARAMS.BKG.phase;
    PARAMS.BKG.phase_cor(bi).TR         = rem(PARAMS.BKG.Rects(bi).TR(1)/(1/PARAMS.spat_freq_pix),floor(PARAMS.BKG.Rects(bi).TR(1)/(1/PARAMS.spat_freq_pix)))*360+PARAMS.BKG.phase;
    PARAMS.BKG.phase_cor(bi).BL         = rem(PARAMS.BKG.Rects(bi).BL(1)/(1/PARAMS.spat_freq_pix),floor(PARAMS.BKG.Rects(bi).BL(1)/(1/PARAMS.spat_freq_pix)))*360+PARAMS.BKG.phase;
    PARAMS.BKG.phase_cor(bi).BR         = rem(PARAMS.BKG.Rects(bi).BR(1)/(1/PARAMS.spat_freq_pix),floor(PARAMS.BKG.Rects(bi).BR(1)/(1/PARAMS.spat_freq_pix)))*360+PARAMS.BKG.phase;

end

Screen('FillRect',win,127); 
DrawFormattedText(win, 'Surround supression experiment\n\n (Press any key to start)','center','center');
Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);
Screen('Flip', win);

KbWait([], 2);
WaitSecs(2)
flagscape =0;
trl = 1;
while trl <= length(trialcond) 
% for trl = 1:length(TRIAL)
    if trl>1 & rem(trl,PARAMS.TRIAL.trial_per_block)==1
        Screen('FillRect',win,127); 
        DrawFormattedText(win, sprintf('Block %d of %d finished\n\n (Press any key to start)',floor(trl/PARAMS.TRIAL.trial_per_block),ceil(length(trialcond)/PARAMS.TRIAL.trial_per_block)),'center','center');
        Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);
        Screen('Flip', win);
        WaitSecs(1)
        KbWait([], 2);
        WaitSecs(1  )
    end
    % choose randomly hich condition to do next from available ones
    thisCOND = trialcond(trl);
    TRIAL(trl).bkg_contrast = CONDpsych(thisCOND).bkg_contrast;
    % this are the trial textures and parameters to use
    % from "CreateProceduralGabor" ‘contrast’ (in 'DrawTexture')
    % For a zero degrees grating:
    % g(x,y) = modulatecolor * contrast * contrastPreMultiplicator * sin(x*2*pi*freq + phase) + Offset
    % so with a grayscale from 0 to 1, for a 0.5 background color offset (specified in CreateProce...),
    % 100% 'contrast' (a sine grating that goes from 0 to 1) the max value of contrast that should be given
    % is 0.5, higher value will go above the range and trasnform a sine
    % grating to bars
    thisBKGcontrast = str2num(TRIAL(trl).bkg_contrast)/100/2; % contrast
    
    %thisGRTcontrast = str2num(TRIAL(trl).grt_contrast)/100/2;
    tIntensity  = QuestQuantile(CONDpsych(thisCOND).q);
    if tIntensity>0
        tIntensity=-rand(1).*.1;
    end
    TRIAL(trl).grt_contrast  = (10^tIntensity)/2;
    
    thisGRTcontrast = TRIAL(trl).grt_contrast;
    % foregorund grating rects and parameters
    thisRects    = [PARAMS.GRT.Rects.TL;PARAMS.GRT.Rects.TR;PARAMS.GRT.Rects.BL;PARAMS.GRT.Rects.BR]';
    % where is the target above or belo
    upordown     =   round(rand(1));
    if upordown ==1
        % left or right
        if rand(1)>.5
            thisGRTcontrastTL =  thisGRTcontrast;
            thisGRTcontrastTR =  0;
        else
            thisGRTcontrastTL =  0;
            thisGRTcontrastTR =  thisGRTcontrast;
        end
        thisauxParam   = [PARAMS.GRT.phase_cor.TL, PARAMS.spat_freq_pix, thisGRTcontrastTL, 0;...
                    PARAMS.GRT.phase_cor.TR, PARAMS.spat_freq_pix, thisGRTcontrastTR, 0;...
                    PARAMS.GRT.phase_cor.BL, PARAMS.spat_freq_pix, 0, 0;...
                    PARAMS.GRT.phase_cor.BR, PARAMS.spat_freq_pix, 0, 0]';
        TRIAL(trl).grt_position_present = 'Up';
    elseif upordown ==0
        if rand(1)>.5
            thisGRTcontrastBL =  thisGRTcontrast;
            thisGRTcontrastBR =  0;
        else
            thisGRTcontrastBL =  0;
            thisGRTcontrastBR =  thisGRTcontrast;
        end
        thisauxParam = [PARAMS.GRT.phase_cor.TL, PARAMS.spat_freq_pix, 0, 0;...
                    PARAMS.GRT.phase_cor.TR, PARAMS.spat_freq_pix, 0, 0;...
                    PARAMS.GRT.phase_cor.BL, PARAMS.spat_freq_pix,  thisGRTcontrastBL, 0;...
                    PARAMS.GRT.phase_cor.BR, PARAMS.spat_freq_pix,  thisGRTcontrastBR, 0]';
         TRIAL(trl).grt_position_present = 'Down';
    end
  
    TRIAL(trl).bkg_orient = CONDpsych(thisCOND).bkg_orient;
    TRIAL(trl).bkg_size  = CONDpsych(thisCOND).bkg_size;
    if strcmp(TRIAL(trl).bkg_size,'none')
        BKGix = 1; % it does not matter, there is no background
    else
        BKGix = find(ismember(PARAMS.TRIAL.bkg_sizes,TRIAL(trl).bkg_size));
    end    
    thisRectsBKG    = [PARAMS.BKG.Rects(BKGix).TL;PARAMS.BKG.Rects(BKGix).TR;PARAMS.BKG.Rects(BKGix).BL;PARAMS.BKG.Rects(BKGix).BR]';
    thisauxParamBKG = [PARAMS.BKG.phase_cor(BKGix).TL, PARAMS.spat_freq_pix, thisBKGcontrast, 0;...
                    PARAMS.BKG.phase_cor(BKGix).TR, PARAMS.spat_freq_pix, thisBKGcontrast, 0;...
                    PARAMS.BKG.phase_cor(BKGix).BL, PARAMS.spat_freq_pix, thisBKGcontrast, 0;...
                    PARAMS.BKG.phase_cor(BKGix).BR, PARAMS.spat_freq_pix, thisBKGcontrast, 0]';

    if strcmp(TRIAL(trl).bkg_orient,'parallel') 
        thisbkgangle  = 0;
         % thisbkgphase = PARAMS.BKG.pase_cor(ismember(PARAMS.TRIAL.bkg_sizes,TRIAL(trl).bkg_size));
    elseif strcmp(TRIAL(trl).bkg_orient,'orthogonal')
        thisbkgangle  = 90;
         % thisbkgphase = PARAMS.BKG.pase_cor(ismember(PARAMS.TRIAL.bkg_sizes,TRIAL(trl).bkg_size));
    elseif strcmp(TRIAL(trl).bkg_orient,'none')
        thisbkgangle  = 0;
        thisbkgphase = 0;
    end
    % target trial parameters
%     if str2num(TRIAL(trl).target)/1000>0
%         tgtflag = 1; % when 1 a trigger can be presented
%         tgtTime = str2num(TRIAL(trl).target)/1000;
%     else
%         tgtflag = 0;
%         tgtTime = NaN;
%     end
    
    % fixation dot presentation
    Screen('FillRect',win,127);
     Screen('Flip', win);
    WaitSecs(1)
    Screen('FillRect',win,127);
    Screen('FrameOval', win,PARAMS.ovalColor , [PARAMS.GRT.Rects.TL',PARAMS.GRT.Rects.TR',PARAMS.GRT.Rects.BL',PARAMS.GRT.Rects.BR'],PARAMS.ovalPenWidth); % fixation dot
    
    Screen('FillOval', win,[255 0 0] , [PARAMS.setup.scr_wdth/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_wdth/2+PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2+PARAMS.centerdotsize/2]);
    %Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);   
    Screen('Flip', win);
    WaitSecs(1)
         
    % start of the trial with grating background, fixation dot and first
    % graing frame presentation, + trigger
      Screen('DrawTextures', win, gratingidBKG(BKGix),[],thisRectsBKG,thisbkgangle ,[],[],[],[],[], thisauxParamBKG); % background
  
    % Screen('DrawTexture', win, thisbkg,[],[],thisbkgangle ,[],[],[],[],[], [thisbkgphase, PARAMS.spat_freq_pix, thisBKGcontrast, 0]); % background
    Screen('FillOval', win,[255 0 0] , [PARAMS.setup.scr_wdth/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_wdth/2+PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2+PARAMS.centerdotsize/2]); % fixation dot
    
    % for testing
    % Screen('DrawText', win, sprintf('GRT contrast %1.3f, BKG contrast %s, size %s, orient %s',TRIAL(trl).grt_contrast,TRIAL(trl).bkg_contrast,TRIAL(trl).bkg_size,TRIAL(trl).bkg_orient) ,...
    %    PARAMS.setup.scr_wdth*.1,PARAMS.setup.scr_hgt*.9, [255 0 0] ); %for testing

%     Screen('FillRect', win, generateVPixxTrigger(TRIAL(trl).trigger,true) , [-5 -5 2 2]); % FOR EEG TRIGGERING

    Screen('DrawTextures', win, gratingid,[],thisRects ,[],[],[],[],[],[], thisauxParam); % the foregorund gratings
    Screen('FrameOval', win,PARAMS.ovalColor , [PARAMS.GRT.Rects.TL',PARAMS.GRT.Rects.TR',PARAMS.GRT.Rects.BL',PARAMS.GRT.Rects.BR'],PARAMS.ovalPenWidth); % fixation dot
    
    [VBLTimestamp] = Screen('Flip', win,[],1); % this is when the first trial (and trigger) flips
%     timeStamps = nan(1200,1);
%     firstFlip = VBLTimestamp;
%     Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);   
 %   Screen('Flip', win, 0,1,0);
    %flag = 1; % this is used to determine which grating frame to flip (grating in top and gray circles bottom or viceversa)
    %ii   = 1;
    
    % %         to make stimuli screenshot
   %         imageArray=Screen('GetImage', win);
 %           fh = figure, imshow(imresize(imageArray,.5));
%            doimage(fh,['/Users/jossando/trabajo/SSVEP_surroundsuppression/02_Planning'],'png',sprintf('psychGRcont%1.3f_BKGcont%s_size%s_orient%s%d',TRIAL(trl).grt_contrast,TRIAL(trl).bkg_contrast,TRIAL(trl).bkg_size,TRIAL(trl).bkg_orient),'300','painters',[],1)


    thistime = GetSecs;
     WaitSecs(PARAMS.TRIAL.stim_duration);
  %  WaitSecs(5);
    % answer screen,
    Screen('FillRect',win,127); 
    Screen('FrameOval', win,PARAMS.ovalColor , [PARAMS.GRT.Rects.TL',PARAMS.GRT.Rects.TR',PARAMS.GRT.Rects.BL',PARAMS.GRT.Rects.BR'],PARAMS.ovalPenWidth); % fixation dot
    Screen('FillOval', win,[255 0 0] , [PARAMS.setup.scr_wdth/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_wdth/2+PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2+PARAMS.centerdotsize/2]);
   
    DrawFormattedText(win, 'Up or Down?\n(Press Esc to abort the experiment)','center',round(PARAMS.setup.scr_hgt*.25),[0 0 0]);
%     Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);   
    Screen('Flip', win);
    while KbCheck; end % Wait until all keys are released.

    while 1
        % Check the state of the keyboard.
        [ keyIsDown, seconds, keyCode ] = KbCheck;

        % If the user is pressing a key, then display its code number and name.
        if keyIsDown
            if strcmp(KbName(keyCode),'UpArrow')
               TRIAL(trl).grt_position_response = 'Up';
               if  strcmp(TRIAL(trl).grt_position_present,'Up')
                   TRIAL(trl).correct =1;
               else
                   TRIAL(trl).correct =0;
               end
                break
            elseif strcmp(KbName(keyCode),'DownArrow')   
                TRIAL(trl).grt_position_response = 'Down';
                if strcmp(TRIAL(trl).grt_position_present,'Down')
                   TRIAL(trl).correct =1;
               else
                   TRIAL(trl).correct =0;
                end
               
                break
            elseif keyCode(escapeKey)
%                 sca
                flagscape = 1;
                break
                %error(sprintf('Experiment aborted by experimenter on trial %d',trl))
            end
        end
    end
    if flagscape
        break
    else
        TRIAL(trl).RT = seconds-thistime;
        CONDpsych(thisCOND).RT =[CONDpsych(thisCOND).RT seconds-thistime];
        if TRIAL(trl).correct == 1
            CONDpsych(thisCOND).q =QuestUpdate(CONDpsych(thisCOND).q,tIntensity,1); 
        elseif TRIAL(trl).correct == 0   
            CONDpsych(thisCOND).q =QuestUpdate(CONDpsych(thisCOND).q,tIntensity,0); 
        end
        CONDpsych(thisCOND).q95 = [CONDpsych(thisCOND).q95 (10.^QuestQuantile(CONDpsych(thisCOND).q,[.975])-10.^QuestQuantile(CONDpsych(thisCOND).q,[.025]))/2];
        if length(CONDpsych(thisCOND).q95)>5
            if sum(CONDpsych(thisCOND).q95(end-2:end)<PARAMS.TRIAL.PsychCurve.q95Limit)==3

                trialcond(trialcond==thisCOND & [1:length(trialcond)]>trl) = [];
            end
        end
        save(fullfile(resultPATH,fileMAT),'PARAMS','TRIAL','CONDpsych')
        trl = trl+1;
    end
end

 if flagscape
     Screen('FillRect',win,127); 
    DrawFormattedText(win, sprintf('Experiment interrupted',floor(trl/PARAMS.TRIAL.trial_per_block),length(TRIAL)/PARAMS.TRIAL.trial_per_block),'center','center');
    Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);
    Screen('Flip', win);
    WaitSecs(2)
    sca
 else
    Screen('FillRect',win,127); 
    DrawFormattedText(win, sprintf('Block %d/%d finished\nExperiment finished!',floor(trl/PARAMS.TRIAL.trial_per_block),length(TRIAL)/PARAMS.TRIAL.trial_per_block),'center','center');
    Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);
    Screen('Flip', win);
    WaitSecs(2)
    sca
 end
%%
ListenChar(0)
if trl>length(CONDpsych)*PARAMS.TRIAL.trial_per_condition/2 % only plot if at least half of the trials were done
    fh = figure;,hold on
    %col = cbrewer('qual','Set1',9);
    col = [.894 .102 .11;.215 .494 .721;.302 .686 .29;.596 .305 .639;1 .498 0;1 1 .2;.65 .33 .156; .968 .505 .749;.894 .102 .11;.215 .494 .721;.302 .686 .29;.596 .305 .639;1 .498 0;1 1 .2;.65 .33 .156; .968 .505 .749];
    for qq = 1:length(CONDpsych)
        rangeI = 0:.005:1;
        indxData = CONDpsych(qq).q.trialCount;
        [lC_unc,dev_unc,stat] = glmfit(10.^CONDpsych(qq).q.intensity(1:indxData)',CONDpsych(qq).q.response(1:indxData)','binomial','logit');
        [logitFit_unc logitFit_unc_lo logitFit_unc_hi] = glmval(lC_unc,rangeI ,'logit',stat);
    
        %TODO se
        
        % vline(0,'k--')
        %hline(.5,'k--')
        plot(10.^CONDpsych(qq).q.intensity(1:indxData)',CONDpsych(qq).q.response(1:indxData)','o','MarkerSize',6,'Color',col(qq,:))
        h(qq) = plot(rangeI ,logitFit_unc,'r','LineWidth',2,'Color',col(qq,:));
        axis([0 rangeI(end) 0 1])
        text(.7,qq/15,sprintf('%s %% bkg %s %s\nthreshold:%1.3f (sd:%1.3f)',CONDpsych(qq).bkg_contrast,CONDpsych(qq).bkg_size,CONDpsych(qq).bkg_orient,10.^QuestMean(CONDpsych(qq).q),10.^QuestSd(CONDpsych(qq).q)),'Color',col(qq,:),'FontSize',6)
    end
    xlabel('Sine grating ''contrast''')
    ylabel('Incorrect(0)/Correct(1')
    doimage(fh,resultPATH,'pdf',[fileMAT(1:end-4) '_curve'],'300','painters',[],0)
end
ListenChar(0)