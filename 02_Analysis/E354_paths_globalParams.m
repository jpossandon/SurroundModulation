paths.Parent        = '/Users/jossando/trabajo/SSVEP_surroundsuppression/';
paths.standardBesa  = '/Users/jossando/trabajo/matlab/eeglab2023.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp';
paths.Raw           = fullfile(paths.Parent,['06_RawData/' SITE]);
paths.Analysis      = fullfile(paths.Parent,['07_Analysis/' SITE]);
paths.preproc       = fullfile(paths.Analysis,'preproc');
% paths.pathbadchannel= fullfile(paths.Parent,'06_RawData/HAMBURG/');

% analysis parameters

% how many std zscore for a channel to be considered to have a reliable SSVEP response
% using 4 for now --> X = norminv(.0001,0,1) = -3.719, .05/75=0.0006 (75 is numbers of channels)
anlys.zscoreThresh = 4;   

% preprocessing
preproc.rsfreq            = 500;

preproc.badCh.segments_size_sec   = 2;
preproc.badCh.dead_thresh_var     = 1;
preproc.badCh.variance_thresh     = 250;
preproc.badCh.range_thresh        = 500;
preproc.badCh.removeAbove         = .15; %amount of dead/noise allowd in a channel

preproc.badChforICA.segments_size_sec   = 2;
preproc.badChforICA.dead_thresh_var     = 1;
preproc.badChforICA.variance_thresh     = 150;
preproc.badChforICA.range_thresh        = 300;


% figures parameters
figs.txtSiz         = 6;
figs.lblSiz         = 8;
figs.figSizBig      = [17.6 17.6];
figs.cmapSet1       = cbrewer('qual','Set1',9);
figs.cmapTopo1      = colormap('turbo');
close(gcf)

figs.cmapPaired       = cbrewer('qual','Paired',12);
figs.cmapTopo1      = colormap('turbo');
close(gcf)

figs.cmapGray       = [.2 .2 .2;1 .6 .6];

figs.scatter.mkSiz   = 20;
figs.scatter.MEC     = [.5 .5 .5];
figs.scatter.MFA     = .5;
figs.scatter.MEA     = .5;
figs.scatter.LW      = .25;
figs.scatter.options = {'LineWidth',figs.scatter.LW ,...
                        'MarkerEdgeColor',figs.scatter.MEC,...
                        'MarkerFaceAlpha',figs.scatter.MFA ,...
                        'MarkerEdgeAlpha',figs.scatter.MEA};

% ADD HERE NEW SITE e..g HAMBURG_children
if strcmp(SITE,'HYDERABAD')
    subjIDs  = {'S1608','S1615','S1618','S1628'}; 
elseif strcmp(SITE,'HYDERABAD_controls')
    subjIDs  = {'S1640','S999999','S9999'}; %HERE ADD NEW PARTICIPANTS
elseif strcmp(SITE,'HYDERABAD_controls_blur')
    subjIDs  = {'S1640','S999999','S9999'}; %HERE ADD NEW PARTICIPANTS
elseif strcmp(SITE,'HAMBURG')
     subjIDs  = {'S0001','S0002','S0003','S0004','S0005','S0006','S0007','S0008','S0009','S0010',...
    'S0011','S0012','S0013','S0014','S0015','S0016','S0017','S0018','S0020','S0021','S0022','S0023',...
    'S0024','S0025'};
elseif strcmp(SITE,'HAMBURGv2')
     subjIDs  = {'S001','S002','S003','S005','S006','S007','S008','S010',...
    'S011','S012','S013','S014','S015','S016','S017','S018','S019','S020','S021','S022','S023','S024','S025','S026','S027','S028','S029','S030'};
end 
