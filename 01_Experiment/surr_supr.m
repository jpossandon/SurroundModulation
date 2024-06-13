
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
clear all
clc

KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');
EXP_LOCATION =3; % 1 - LVPEI 2 - Hamburg 3 - JPO laptop

if EXP_LOCATION == 1
    resultPATH = 'C:\Users\BPN\Documents\EEG_Experiments\Surround_supression\RawData'; % LVPEI
elseif EXP_LOCATION == 2
    resultPATH = 'C:\Users\bpn\Documents\Jose\Surround_supression\RawData'; % Babylab EEG
elseif EXP_LOCATION == 3
    resultPATH = '/Users/jossando/trabajo/SSVEP_surroundsuppression/06_RawData/'; % jose MAC
end

testType = questdlg('', ...
	'Surround supression EEG', ...
	'Experiment','Test Trials','Experiment');

if strcmp(testType,'Test Trials')
    answer      = {'ERASEME';'jose'};
    fileMAT     = 'ERASEME.mat';
elseif strcmp(testType,'Experiment')
    answer=inputdlg({'Participant ID:','Experimenter:'},'SS SSVEP');
    fileMAT             = ['CS_' answer{1} '.mat'];
    if exist(fullfile(resultPATH,fileMAT))
      selection = questdlg(sprintf('%s already exists. Do you want to:',fileMAT),'CS SSVEP',...
          'Cancel','Overwrite','Save as new','Cancel');
      switch selection
          case 'Cancel'
              error(sprintf('Experiment cancelled by experimenter, %s already exists',fileMAT))
          case 'Overwrite'
              selection2 = questdlg(sprintf('Are you sure you want to overwrite %s?',fileMAT),'CS SSVEP',...
          'Cancel','Overwrite','Cancel');
                switch selection2
                    case 'Cancel'
                        error(sprintf('Experiment cancelled by experimenter, %s already exists',fileMAT))
                    case 'Overwrite'
                        sprintf('%s will be overwritten. Previous file is save in ~/overwritten as CS_%sbackup.mat',fileMAT,answer{1})
                        mkdir(fullfile(resultPATH,'overwritten'))
                        copyfile(fullfile(resultPATH,fileMAT),fullfile(resultPATH,'overwritten',['CS_' answer{1} 'backup.mat']))
                end
          case 'Save as new'
              fileMAT = ['CS_' answer{1} '_' datestr(now,'HH_MM') '.mat']; 
              selection3 = questdlg(sprintf('File will be saved as %s (current time)?',fileMAT),'CS SSVEP',...
          'Cancel','Confirm','Cancel');
              if strcmp(selection3,'Cancel')
                error(sprintf('Experiment cancelled by experimenter, %s already exists',['CS_' answer{1} '.mat']))
              end
      end
    end
end
PARAMS.subjectID    = answer{1};
PARAMS.experimenter = answer{2};
PARAMS.date         = datestr(now);
clear answer selection selection2 selection3 
Screen('Preference', 'SkipSyncTests', 1);
Screen('Preference', 'ConserveVRAM',  4096)
 
 % Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;


%Screen('Preference', 'SkipSyncTests', 1);
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
slack                           = PARAMS.setup.ifi/2;
PARAMS.setup.refreshRate        = round(1/PARAMS.setup.ifi);
PARAMS.flick_rate               = 15; % in Hz
if rem(PARAMS.setup.refreshRate,PARAMS.flick_rate)==surr0
    PARAMS.frames_refresh       = PARAMS.setup.refreshRate/PARAMS.flick_rate/2;
    fprintf('\n%%%%%%%%%%%%%%%%%%%%\nMonitor refresh rate: %d Hz\nRequested stimulation rate: %d Hz\nScreen will flip every %d frames.\n',PARAMS.setup.refreshRate,PARAMS.flick_rate,PARAMS.frames_refresh)
else
    error('Not possible to stimulate at %d Hz with a monitor refresh rate of %d',PARAMS.flick_rate,PARAMS.setup.refreshRate )
end
WaitSecs(5)


Screen('FillRect',win,[127 127 127]); 
Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);
Screen('Flip', win);

PARAMS.setup.distance_to_screen = 60; % in cms
PARAMS.setup.screen_width       = 53; % in cms (viewPixx EEG 24-inch display size (diagonal),53.5 cm recchek)
PARAMS.setup.pixels_per_degree  = PARAMS.setup.scr_wdth/(2*atan(PARAMS.setup.screen_width/2/PARAMS.setup.distance_to_screen)*180/pi); % TODO acutal calcualtion

% GRATING AND BACKGROUND PARAMETERS
PARAMS.spat_freq_deg            = 1; % in cycles per degree (Sourav, 2018 used 2 cycles per degree)
PARAMS.spat_freq_pix            = PARAMS.spat_freq_deg./PARAMS.setup.pixels_per_degree; % in cycles per pixel
PARAMS.centerdotsize            = PARAMS.setup.pixels_per_degree; %in pixels

% GRATING PARAMETERS
PARAMS.GRT.radius               = 2.5/2*PARAMS.setup.pixels_per_degree; % diameter 2.5 visual degrees
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
PARAMS.TRIAL.trial_per_condition= 6;
PARAMS.TRIAL.trial_per_block    = 10;
PARAMS.TRIAL.trial_duration     = 11; % in seconds
if strcmp(testType,'Test Trials')
    PARAMS.TRIAL.colortarget_trials = 32; % as a percentage of trials
else
    PARAMS.TRIAL.colortarget_trials = 16; % as a percentage of trials
end
PARAMS.TRIAL.target_duration    = .5; % in seconds
PARAMS.TRIAL.grt_positions      = {'peripheral'}; % peripheral,(foveal?)
PARAMS.TRIAL.grt_contrasts      = {'50' '100'}; % 50,100
PARAMS.TRIAL.bkg_contrasts      = {'0' '100'}; % 0,(50?)100
PARAMS.TRIAL.bkg_sizes          = {'2','4','full'};% size of bakground. full is the complete screen, otherwise number n specify a radial background thath extend (n-1)*GRT.radius beyond the grating outermos edge 
PARAMS.TRIAL.bkg_orients        = {'parallel','orthogonal'}; % parallel or orthogonal

% BACKGROUND PARAMETERS
PARAMS.BKG.angle                = [0 90];
%PARAMS.BKG.contrast    = [0 0.5]; 
PARAMS.BKG.phase                = 0;
PARAMS.BKG.widths               = [];
PARAMS.BKG.heights              = [];
for bi = 1:length(PARAMS.TRIAL.bkg_sizes)
    if strcmp(PARAMS.TRIAL.bkg_sizes(bi),'full')
        PARAMS.BKG.widths(bi)   = PARAMS.setup.scr_wdth;
        PARAMS.BKG.heights(bi)  = PARAMS.setup.scr_hgt;
    else
        PARAMS.BKG.widths(bi)   = round(PARAMS.setup.screen_ratio*(PARAMS.GRT.eccenticity*PARAMS.setup.pixels_per_degree+PARAMS.GRT.radius*str2num(PARAMS.TRIAL.bkg_sizes{bi}))*2);
        PARAMS.BKG.heights(bi)  = round((PARAMS.GRT.eccenticity*PARAMS.setup.pixels_per_degree+PARAMS.GRT.radius*str2num(PARAMS.TRIAL.bkg_sizes{bi}))*2);
    end
end
                        
% TRIAL STRUCTURE
% all combination of factors
pcombs              = cell(1,5);
[pcombs{:}]         = ndgridnotfull(PARAMS.TRIAL.grt_positions,PARAMS.TRIAL.grt_contrasts,PARAMS.TRIAL.bkg_contrasts,PARAMS.TRIAL.bkg_sizes,PARAMS.TRIAL.bkg_orients);
postrials           = cat(2,pcombs{1}(:),pcombs{2}(:),pcombs{3}(:),pcombs{4}(:),pcombs{5}(:));
postrials(ismember(postrials(:,3),'0'),[4 5]) = {'none'}; % change bakground factors that are not possible for contrast 0 to 'none'
postrials           = table2cell(unique(cell2table(postrials)));  % remove trials wth the same factors (becaue of liena bove)
postrials           = [postrials,num2cell(100+[2:2:2*size(postrials,1)])']; % triggers
postrials           = repmat(postrials,PARAMS.TRIAL.trial_per_condition,1); % expand for conditon repetitions
%ixtarget    = randsample(size(postrials,1),ceil(size(postrials,1)*100/(100-PARAMS.TRIAL.colortarget_trials))-size(postrials,1)); % select randomly conditions to be target trials
ixtarget            = randperm(size(postrials,1));
ixtarget            = ixtarget(1:ceil(size(postrials,1)*100/(100-PARAMS.TRIAL.colortarget_trials))-size(postrials,1));
tgttrials           = postrials(ixtarget,:);
tgttrials(:,end)    = num2cell(cell2mat(tgttrials(:,end))+100);
postrials           = [postrials;tgttrials]; % add target trials
targettimes = num2cell(round(rand(length(ixtarget),1)*(PARAMS.TRIAL.trial_duration-.75)*1000));%im ms
postrials           = [postrials,[repmat({'0'},size(postrials,1)-length(ixtarget),1);cellfun(@num2str,targettimes,'UniformOutput',false)]]; % add column indicating tjhat trial is target
%postrials   = postrials(randsample(size(postrials,1),size(postrials,1)),:); % and randomize order of trials
postrials           = postrials(randperm(size(postrials,1)),:); % and randomize order of trials

%structure used during the experimetns and saved
TRIAL               = struct('grt_position',postrials(:,1),'grt_contrast',postrials(:,2), 'bkg_contrast',postrials(:,3), 'bkg_size',postrials(:,4), 'bkg_orient',postrials(:,5), 'trigger',postrials(:,6), 'target',postrials(:,7),'answerColorChange',[],'correct',[]);
clear pcombs postrials ixtarget targettimes tgttrials
    


% CREATE TEXTURES
% backgrounds
for bi = 1:length(PARAMS.TRIAL.bkg_sizes)
    [bkgidvert(bi), vertrect]   = CreateProceduralSineGrating(win, PARAMS.BKG.widths(bi), PARAMS.BKG.heights(bi),   [0.5 0.5 0.5 0.0],PARAMS.BKG.heights(bi)/2, []);
    [bkgidhorz(bi), horzrect]   = CreateProceduralSineGrating(win, PARAMS.BKG.heights(bi), PARAMS.BKG.widths(bi),  [0.5 0.5 0.5 0.0],PARAMS.BKG.heights(bi)/2, []);
    
    PARAMS.BKG.Rects(bi).vert   = OffsetRect(CenterRectOnPoint(vertrect,0,0),PARAMS.setup.scr_wdth/2,PARAMS.setup.scr_hgt/2);
    PARAMS.BKG.Rects(bi).horz   = OffsetRect(CenterRectOnPoint(horzrect,0,0),PARAMS.setup.scr_wdth/2,PARAMS.setup.scr_hgt/2);
    if strcmp(PARAMS.TRIAL.bkg_sizes(bi),'full')
        PARAMS.BKG.pase_cor(bi) = 0;
    else
        PARAMS.BKG.pase_cor(bi) = rem(PARAMS.BKG.Rects(bi).vert(1)/(1/PARAMS.spat_freq_pix),floor(PARAMS.BKG.Rects(bi).vert(1)/(1/PARAMS.spat_freq_pix)))*360+PARAMS.GRT.phase; 
    end
end
%grating
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


Screen('FillRect',win,[127 127 127]); 
DrawFormattedText(win, 'Surround supression experiment\n\n (Press any key to start)','center','center');
Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);
Screen('Flip', win);
WaitSecs(2)
KbWait([], 2);
WaitSecs(2)
    

for trl = 1:length(TRIAL)
    if trl>1 & rem(trl,PARAMS.TRIAL.trial_per_block)==1
        Screen('FillRect',win,[127 127 127]); 
        DrawFormattedText(win, sprintf('Block %d/%d finished\n\n (Press any key to start)',floor(trl/PARAMS.TRIAL.trial_per_block),floor(length(TRIAL)/PARAMS.TRIAL.trial_per_block)),'center','center');
        Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);
        Screen('Flip', win);
        WaitSecs(2)
        KbWait([], 2);
        WaitSecs(2)
    end
    % this are the trial textures and parameters to use
    thisBKGcontrast = str2num(TRIAL(trl).bkg_contrast)/100/2; % contrast
    thisGRTcontrast = str2num(TRIAL(trl).grt_contrast)/100/2;
    
    % foregorund grating rects and parameters
    thisRects    = [PARAMS.GRT.Rects.TL;PARAMS.GRT.Rects.TR;PARAMS.GRT.Rects.BL;PARAMS.GRT.Rects.BR]';
    thisauxParam1 = [PARAMS.GRT.phase_cor.TL, PARAMS.spat_freq_pix, thisGRTcontrast, 0;...
                    PARAMS.GRT.phase_cor.TR, PARAMS.spat_freq_pix, thisGRTcontrast, 0;...
                    PARAMS.GRT.phase_cor.BL, PARAMS.spat_freq_pix, 0, 0;...
                    PARAMS.GRT.phase_cor.BR, PARAMS.spat_freq_pix, 0, 0]';
    thisauxParam2 = [PARAMS.GRT.phase_cor.TL, PARAMS.spat_freq_pix, 0, 0;...
                    PARAMS.GRT.phase_cor.TR, PARAMS.spat_freq_pix, 0, 0;...
                    PARAMS.GRT.phase_cor.BL, PARAMS.spat_freq_pix,  thisGRTcontrast, 0;...
                    PARAMS.GRT.phase_cor.BR, PARAMS.spat_freq_pix,  thisGRTcontrast, 0]';
    % background grating rects and parameters
    if strcmp(TRIAL(trl).bkg_orient,'parallel') 
        thisbkg       = bkgidvert(ismember(PARAMS.TRIAL.bkg_sizes,TRIAL(trl).bkg_size));
        thisbkgangle  = 0;
         thisbkgphase = PARAMS.BKG.pase_cor(ismember(PARAMS.TRIAL.bkg_sizes,TRIAL(trl).bkg_size));
    elseif strcmp(TRIAL(trl).bkg_orient,'orthogonal')
        thisbkg       = bkgidhorz(ismember(PARAMS.TRIAL.bkg_sizes,TRIAL(trl).bkg_size));
        thisbkgangle  = 90;
         thisbkgphase = PARAMS.BKG.pase_cor(ismember(PARAMS.TRIAL.bkg_sizes,TRIAL(trl).bkg_size));
    elseif strcmp(TRIAL(trl).bkg_orient,'none')
        thisbkg       = bkgidhorz(1); % doesn matter because contrast is 0
        thisbkgangle  = 0;
        thisbkgphase = 0;
    end
    
    % target trial parameters
    if str2num(TRIAL(trl).target)/1000>0
        tgtflag = 1; % when 1 a trigger can be presented
        tgtTime = str2num(TRIAL(trl).target)/1000;
    else
        tgtflag = 0;
        tgtTime = NaN;
    end
    
    % fixation dot presentation
    Screen('FillRect',win,[127 127 127]);
    Screen('FillOval', win,[255 0 0] , [PARAMS.setup.scr_wdth/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_wdth/2+PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2+PARAMS.centerdotsize/2]);
    Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);   
    Screen('Flip', win);
    WaitSecs(2)
         
    % start of the trial with grating background, fixation dot and first
    % graing frame presentation, + trigger
    Screen('DrawTexture', win, thisbkg,[],[],thisbkgangle ,[],[],[],[],[], [thisbkgphase, PARAMS.spat_freq_pix, thisBKGcontrast, 0]); % background
    Screen('FillOval', win,[255 0 0] , [PARAMS.setup.scr_wdth/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_wdth/2+PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2+PARAMS.centerdotsize/2]); % fixation dot
    
    % for testing
    % Screen('DrawText', win, sprintf('GRT contrast %s, BKG contrast %s, size %s, orient %s, target %s, trigger %d',TRIAL(trl).grt_contrast,TRIAL(trl).bkg_contrast,TRIAL(trl).bkg_size,TRIAL(trl).bkg_orient,TRIAL(trl).target,TRIAL(trl).trigger) ,...
       % PARAMS.setup.scr_wdth*.1,PARAMS.setup.scr_hgt*.9, [255 0 0] ); %for testing

    Screen('FillRect', win, generateVPixxTrigger(TRIAL(trl).trigger,true) , [-5 -5 2 2]); % FOR EEG TRIGGERING
%    sprintf('GRT contrast %s, BKG contrast %s, size %s, orient %s, target %s, trigger %d',TRIAL(trl).grt_contrast,TRIAL(trl).bkg_contrast,TRIAL(trl).bkg_size,TRIAL(trl).bkg_orient,TRIAL(trl).target,TRIAL(trl).trigger) % for the experimenter command window
    Screen('DrawTextures', win, gratingid,[],thisRects ,[],[],[],[],[],[], thisauxParam1); % the foregorund gratings
    [VBLTimestamp] = Screen('Flip', win,[],1); % this is when the first trial (and trigger) flips
%     timeStamps = nan(1200,1);
    firstFlip = VBLTimestamp;
    Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);   
    Screen('Flip', win, 0,1,0);
    flag = 1; % this is used to determine which grating frame to flip (grating in top and gray circles bottom or viceversa)
    ii   = 1;
    thistime = GetSecs;tic
    while GetSecs<thistime+PARAMS.TRIAL.trial_duration 
        % this is for target trials
        if tgtflag==1 && GetSecs>firstFlip+tgtTime 
            Screen('FillOval', win,[0 0 255] , [PARAMS.setup.scr_wdth/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_wdth/2+PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2+PARAMS.centerdotsize/2]);
            tgtflag = 2;
            
        end
        if tgtflag==2 && GetSecs>firstFlip+tgtTime+PARAMS.TRIAL.target_duration
            Screen('DrawTexture', win, thisbkg,[],[],thisbkgangle ,[],[],[],[],[], [thisbkgphase, PARAMS.spat_freq_pix, thisBKGcontrast, 0]);
            Screen('FillOval', win,[255 0 0] , [PARAMS.setup.scr_wdth/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_wdth/2+PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2+PARAMS.centerdotsize/2]);
            tgtflag = 0;
        
        end
        % which grating to present next
        if flag==1
            Screen('DrawTextures', win, gratingid,[],thisRects ,[],[],[],[],[],[], thisauxParam2);
            Screen('FillRect', win, generateVPixxTrigger(30,true) , [-5 -5 2 2]); % FOR EEG TRIGGERING
            flag = 2;
        elseif flag==2   
            Screen('DrawTextures', win, gratingid,[],thisRects ,[],[],[],[],[],[], thisauxParam1);
            Screen('FillRect', win, generateVPixxTrigger(80,true) , [-5 -5 2 2]); % FOR EEG TRIGGERING
            flag = 1;
        end
        Screen('DrawingFinished', win,1);
       % toc
        %WaitSecs((PARAMS.frames_refresh-2).*PARAMS.setup.ifi); Hamburg
        WaitSecs((PARAMS.frames_refresh-2).*PARAMS.setup.ifi); %LVPEI
        [VBLTimestamp] = Screen('Flip', win, VBLTimestamp+PARAMS.frames_refresh.*PARAMS.setup.ifi-slack,1,0); % here is when the change of texture is ordered, to occur  PARAMS.frames_refresh frames after last one (slack is used to order this to occur half of a frame before is expeted so in the next monitor refresh occurs)
    %[VBLTimestamp] = Screen('Flip', win, 0,1,0); % here is when the change of texture is ordered, to occur  PARAMS.frames_refresh frames after last one (slack is used to order this to occur half of a frame before is expeted so in the next monitor refresh occurs)
   % [VBLTimestamp] = Screen('AsyncFlipBegin', win, VBLTimestamp+0.03,1,0); % here is when the change of texture is ordered, to occur  PARAMS.frames_refresh frames after last one (slack is used to order this to occur half of a frame before is expeted so in the next monitor refresh occurs)
   %  Screen('AsyncFlipEnd', win);
        Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);   
        Screen('Flip', win, 0,1,0);
          
% %         to make stimuli screenshot
%         if ii<3 
%           imageArray=Screen('GetImage', win);
%           fh = figure, imshow(imresize(imageArray,.5))
%           doimage(fh,['/Users/jossando/trabajo/SSVEP_surroundsuppression/02_Planning'],'png',sprintf('GRcont%s_BKGcont%s_size%s_orient%s%d',TRIAL(trl).grt_contrast,TRIAL(trl).bkg_contrast,TRIAL(trl).bkg_size,TRIAL(trl).bkg_orient,ii),'300','painters',[],1)
%         end
%         ii = ii+1;
   % befoer Flip
    end
    % anser screen
    Screen('FillRect',win,[127 127 127]); 
    DrawFormattedText(win, 'Color change?(y/n)\n(Press Esc to abort the experiment)','center','center',[0 0 0]);
    Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);   
    Screen('Flip', win);
    while KbCheck; end % Wait until all keys are released.

    while 1
        % Check the state of the keyboard.
        [ keyIsDown, seconds, keyCode ] = KbCheck;

        % If the user is pressing a key, then display its code number and name.
        if keyIsDown
            if strcmp(KbName(keyCode),'y')
               TRIAL(trl).answerColorChange = 1;
               if str2num(TRIAL(trl).target)>0
                   TRIAL(trl).correct =1;
               else
                   TRIAL(trl).correct =0;
               end
               break
            elseif strcmp(KbName(keyCode),'n')   
                TRIAL(trl).answerColorChange = 0;
                if str2num(TRIAL(trl).target)==0
                   TRIAL(trl).correct =1;
               else
                   TRIAL(trl).correct =0;
               end
                break
            elseif keyCode(escapeKey)
                sca
                error(sprintf('Experiment aborted by experimenter on trial %d',trl))
            end
        end
    end
    save(fullfile(resultPATH,fileMAT),'PARAMS','TRIAL')
end

Screen('FillRect',win,[127 128 127]); 
DrawFormattedText(win, sprintf('Block %d/%d finished\nExperiment finished!',floor(trl/PARAMS.TRIAL.trial_per_block),length(TRIAL)/PARAMS.TRIAL.trial_per_block),'center','center');
Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);
Screen('Flip', win);
WaitSecs(10)
sca
