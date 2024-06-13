
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

EXP_LOCATION = 2; % 1 - LVPEI 2 - Hamburg 3 - JPO laptop

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

    answer=inputdlg({'Participant ID:','Experimenter:'},'SS Visibility');
    fileMAT             = ['SSvis_' answer{1} '.mat'];
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
                         copyfile(fullfile(resultPATH,fileMAT),fullfile(resultPATH,'overwritten',['SSvis_' answer{1} 'backup.mat']))
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
Screen('Preference', 'SkipSyncTests', 0); % zero during expriment
Screen('Preference', 'ConserveVRAM',  4096)
ListenChar(2) 
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
%slack                           = PARAMS.setup.ifi/2;
%PARAMS.setup.refreshRate        = round(1/PARAMS.setup.ifi);
%PARAMS.flick_rate               = 15; % in Hz
% if rem(PARAMS.setup.refreshRate,PARAMS.flick_rate)==0
%     PARAMS.frames_refresh       = PARAMS.setup.refreshRate/PARAMS.flick_rate/2;
%     fprintf('\n%%%%%%%%%%%%%%%%%%%%\nMonitor refresh rate: %d Hz\nRequested stimulation rate: %d Hz\nScreen will flip every %d frames.\n',PARAMS.setup.refreshRate,PARAMS.flick_rate,PARAMS.frames_refresh)
% else
%     error('Not possible to stimulate at %d Hz with a monitor refresh rate of %d',PARAMS.flick_rate,PARAMS.setup.refreshRate )
% end
% WaitSecs(5)


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

PARAMS.TRIAL.trial_per_block     = 10;
PARAMS.TRIAL.trial_per_condition = 10; % makeit a multiple of trial_per_block
%PARAMS.TRIAL.trial_duration     = 11; % in seconds
PARAMS.TRIAL.stim_duration      = stimdur; % in seconds
%PARAMS.TRIAL.colortarget_trials = 16; % as a percentage of trials
%PARAMS.TRIAL.target_duration    = .5; % in seconds
PARAMS.TRIAL.grt_positions      = {'peripheral'}; % peripheral,(foveal?)
PARAMS.TRIAL.grt_contrasts      = {'100'}; % 50,100
PARAMS.TRIAL.grt_orients        = {'parallel','orthogonal'}; % parallel or orthogonal
                        
% possible conditions
% all combination of factors
pcombs              = cell(1,2);
%[pcombs{:}]         = ndgridnotfull(PARAMS.TRIAL.grt_positions,PARAMS.TRIAL.grt_orients);
pcombs = {{'peripheral','peripheral'},{'parallel','orthogonal'}};
postrials           = cat(2,pcombs{1}(:),pcombs{2}(:));
postrials           = repmat(postrials,PARAMS.TRIAL.trial_per_condition,1); % expand for conditon repetitions
postrials           = postrials(randperm(size(postrials,1)),:); % and randomize order of trials

%structure used during the experimetns and saved
TRIAL               = struct('grt_position',postrials(:,1), 'grt_orient',postrials(:,2),'grt_position_present',[],'answerOrientation',[],'correct',[]);
clear pcombs postrials ixtarget targettimes tgttrials
clear pcombs ixtarget targettimes tgttrials
    

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
thisGRTcontrast = str2num(PARAMS.TRIAL.grt_contrasts{1})/100/2;

Screen('FillRect',win,127); 
DrawFormattedText(win, 'Surround supression visibility test\n\n (Press any key to start)','center','center');
Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);
Screen('Flip', win);

KbWait([], 2);
WaitSecs(2)
flagscape =0;  
 for trl = 1:length(TRIAL)
    if trl>1 & rem(trl,PARAMS.TRIAL.trial_per_block)==1
        Screen('FillRect',win,127); 
        DrawFormattedText(win, sprintf('Block %d/%d finished\n\n (Press any key to start)',floor(trl/PARAMS.TRIAL.trial_per_block),length(TRIAL)/PARAMS.TRIAL.trial_per_block),'center','center');
        Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);
        Screen('Flip', win);
        WaitSecs(1)
        KbWait([], 2);
        WaitSecs(1  )
    end
    
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
  
    % background grating rects and parameters
    if strcmp(TRIAL(trl).grt_orient,'parallel') 
        thisgrtangle  = 0;
    elseif strcmp(TRIAL(trl).grt_orient,'orthogonal')
        thisgrtangle  = 90;
    end
    

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
    
    % for testing
    % Screen('DrawText', win, sprintf('GRT contrast %1.3f, BKG contrast %s, size %s, orient %s',TRIAL(trl).grt_contrast,TRIAL(trl).bkg_contrast,TRIAL(trl).bkg_size,TRIAL(trl).bkg_orient) ,...
    %    PARAMS.setup.scr_wdth*.1,PARAMS.setup.scr_hgt*.9, [255 0 0] ); %for testing

%     Screen('FillRect', win, generateVPixxTrigger(TRIAL(trl).trigger,true) , [-5 -5 2 2]); % FOR EEG TRIGGERING
    Screen('FillOval', win,[255 0 0] , [PARAMS.setup.scr_wdth/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_wdth/2+PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2+PARAMS.centerdotsize/2]);
    Screen('DrawTextures', win, gratingid,[],thisRects ,thisgrtangle,[],[],[],[],[], thisauxParam); % the foregorund gratings
    Screen('FrameOval', win,PARAMS.ovalColor , [PARAMS.GRT.Rects.TL',PARAMS.GRT.Rects.TR',PARAMS.GRT.Rects.BL',PARAMS.GRT.Rects.BR'],PARAMS.ovalPenWidth); % fixation dot
    
    [VBLTimestamp] = Screen('Flip', win,[],1); % this is when the first trial (and trigger) flips

    thistime = GetSecs;
     WaitSecs(PARAMS.TRIAL.stim_duration);
  %  WaitSecs(5);
    % answer screen,
    Screen('FillRect',win,127); 
%     Screen('FrameOval', win,PARAMS.ovalColor , [PARAMS.GRT.Rects.TL',PARAMS.GRT.Rects.TR',PARAMS.GRT.Rects.BL',PARAMS.GRT.Rects.BR'],PARAMS.ovalPenWidth); % fixation dot
%     Screen('FillOval', win,[255 0 0] , [PARAMS.setup.scr_wdth/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2-PARAMS.centerdotsize/2,PARAMS.setup.scr_wdth/2+PARAMS.centerdotsize/2,PARAMS.setup.scr_hgt/2+PARAMS.centerdotsize/2]);
   
    DrawFormattedText(win, 'Vertical (UpArrow) or horizontal (RightArrow)?\n(Press Esc to abort the experiment)','center',round(PARAMS.setup.scr_hgt*.25),[0 0 0]);
%     Screen('FillRect', win, generateVPixxTrigger(0,true) , [-5 -5 2 2]);   
    Screen('Flip', win);
    while KbCheck; end % Wait until all keys are released.

    while 1
        % Check the state of the keyboard.
        [ keyIsDown, seconds, keyCode ] = KbCheck;

        % If the user is pressing a key, then display its code number and name.
        if keyIsDown
            if strcmp(KbName(keyCode),'UpArrow')
               TRIAL(trl).grt_position_response = 'vertical';
               if  strcmp(TRIAL(trl).grt_orient,'parallel')
                   TRIAL(trl).correct =1;
               else
                   TRIAL(trl).correct =0;
               end
                break
            elseif strcmp(KbName(keyCode),'RightArrow')   
                TRIAL(trl).grt_position_response = 'horizontal';
                if strcmp(TRIAL(trl).grt_orient,'orthogonal')
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
        
        save(fullfile(resultPATH,fileMAT),'PARAMS','TRIAL')
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
fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%d CORRECT OF %d TRIALS\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',sum([TRIAL.correct]),sum([TRIAL.correct])+sum(~[TRIAL.correct]))

