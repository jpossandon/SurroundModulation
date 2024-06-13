    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E354 surround supression behavioral analysis
% participants behavioral file ('SSpsych_xxxx) must be in each participant
% individual folder (~/06_RawData/<SITE>/Sxxxx/)(don't forget the 'S')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% where does the data comes from?
%SITE = 'HYDERABAD_controls';
 SITE = 'HAMBURGv2';

E354_paths_globalParams    % script with all paths and parameters for analysis and figures, this file need to be changed accordingly

MAKE_RESULTFILE_INDIVIDUALFIGURES = 0;
% list of subjects to analyze
%%
if MAKE_RESULTFILE_INDIVIDUALFIGURES

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
    
    % table variables
    subjId = {}; noBkg = []; smallOrtho = []; fullOrtho = []; smallParallel = []; fullParallel = [];
    for subj = 1:length(subjIDs)
        % GET INDIVIDUAL FILES
        if strcmp(SITE,'HYDERABAD')
            load(fullfile(paths.Raw,subjIDs{subj},['SSpsych_' subjIDs{subj}(2:end)]))
        elseif strcmp(SITE,'HYDERABAD_controls')
            load(fullfile(paths.Raw,subjIDs{subj},['SSpsych_v2_' subjIDs{subj} '_Normal_Vision']))
        elseif strcmp(SITE,'HYDERABAD_controls_blur')
            load(fullfile(paths.Raw,subjIDs{subj},['SSpsych_v2_' subjIDs{subj} '_Low_Vision']))
        elseif strcmp(SITE,'HAMBURG')
             load(fullfile(paths.Raw,subjIDs{subj},['SSpsych_' subjIDs{subj}]))
        elseif strcmp(SITE,'HAMBURGv2')
             load(fullfile(paths.Raw,subjIDs{subj},['E354v2_SSpsych_' subjIDs{subj}]))
        end     
        
        % TODO check stimulus duration, for now just print it to the command
        % window
        fprintf('\nParticipant %s stimulus duration was %d ms',subjIDs{subj},PARAMS.TRIAL.stim_duration*1000)
        subjId{subj} = subjIDs{subj};

        figure
        for qq = 1:length(CONDpsych)            %loop thorough conditions
                %INDIVIDUAL PARTICIPANT FIGURE
                tcount = CONDpsych(qq).q.trialCount;                            % trials in a given condition
                xT = 1:tcount;                                                  % to plot them sequentially
                intensities = CONDpsych(qq).q.intensity(1:tcount);              % the trials contrast
                responses = CONDpsych(qq).q.response(1:tcount);
                subplot(3,3,qq), hold on    
                colororder({'k','b'})
                 yyaxis left
                plot(xT(responses==1),100*10.^intensities(responses==1),'.k')   % correct trials
                plot(xT(responses==0),100*10.^intensities(responses==0),'.r')   % incorrect trials
               ylim([0 100])
               if qq==1 || qq==4 || qq==7,ylabel('Contrast'),end
                 yyaxis right
                  plot(xT,CONDpsych(qq).q95(1:tcount),'.-b')   % correct trials
                  if qq==3 || qq==6 || qq==9,ylabel('q95'),end
                ylim([0 10])
                

                title(sprintf('BKG:siz %s %s cont %s\n82%% thesh: %3.1f q95:%1.3f',...
                    CONDpsych(qq).bkg_size,CONDpsych(qq).bkg_orient,CONDpsych(qq).bkg_contrast,100*10.^QuestMode(CONDpsych(qq).q),...
                    (10.^QuestQuantile(CONDpsych(qq).q,[.975])-10.^QuestQuantile(CONDpsych(qq).q,[.025]))/2))
                

                % this is to gather each participant/condition data
                % together
                if ismember(SITE,{'HYDERABAD','HAMBURG'})
                    if strcmp(CONDpsych(qq).bkg_size,'none') && strcmp(CONDpsych(qq).bkg_orient,'none')
                        noBkg(subj) = 100*10.^QuestMean(CONDpsych(qq).q);
                    elseif strcmp(CONDpsych(qq).bkg_size,'2') && strcmp(CONDpsych(qq).bkg_orient,'parallel')
                        smallParallel(subj) = 100*10.^QuestMean(CONDpsych(qq).q);    
                    elseif strcmp(CONDpsych(qq).bkg_size,'2') && strcmp(CONDpsych(qq).bkg_orient,'orthogonal')
                        smallOrtho(subj) = 100*10.^QuestMean(CONDpsych(qq).q);    
                    elseif strcmp(CONDpsych(qq).bkg_size,'full') && strcmp(CONDpsych(qq).bkg_orient,'parallel')
                        fullParallel(subj) = 100*10.^QuestMean(CONDpsych(qq).q);    
                    elseif strcmp(CONDpsych(qq).bkg_size,'full') && strcmp(CONDpsych(qq).bkg_orient,'orthogonal')
                        fullOrtho(subj) = 100*10.^QuestMean(CONDpsych(qq).q); 
                    end
                end
                 if ismember(SITE,{'HAMBURGv2','HYDERABAD_controls','HYDERABAD_controls_blur'})
                    if strcmp(CONDpsych(qq).bkg_size,'none') && strcmp(CONDpsych(qq).bkg_orient,'none')
                        noBkg(subj) = 100*10.^QuestMean(CONDpsych(qq).q);
                    elseif strcmp(CONDpsych(qq).bkg_size,'.25') && strcmp(CONDpsych(qq).bkg_orient,'parallel') && strcmp(CONDpsych(qq).bkg_contrast,'10')
                        size025ParallelCont10(subj) = 100*10.^QuestMean(CONDpsych(qq).q);    
                    elseif strcmp(CONDpsych(qq).bkg_size,'.25') && strcmp(CONDpsych(qq).bkg_orient,'parallel') && strcmp(CONDpsych(qq).bkg_contrast,'30')
                        size025ParallelCont30(subj) = 100*10.^QuestMean(CONDpsych(qq).q);    
                    elseif strcmp(CONDpsych(qq).bkg_size,'1') && strcmp(CONDpsych(qq).bkg_orient,'parallel') && strcmp(CONDpsych(qq).bkg_contrast,'10')
                        size1ParallelCont10(subj) = 100*10.^QuestMean(CONDpsych(qq).q);    
                    elseif strcmp(CONDpsych(qq).bkg_size,'1') && strcmp(CONDpsych(qq).bkg_orient,'parallel') && strcmp(CONDpsych(qq).bkg_contrast,'30')
                        size1ParallelCont30(subj) = 100*10.^QuestMean(CONDpsych(qq).q);    
                    elseif strcmp(CONDpsych(qq).bkg_size,'.25') && strcmp(CONDpsych(qq).bkg_orient,'orthogonal') && strcmp(CONDpsych(qq).bkg_contrast,'10')
                        size025OrthoCont10(subj) = 100*10.^QuestMean(CONDpsych(qq).q);    
                    elseif strcmp(CONDpsych(qq).bkg_size,'.25') && strcmp(CONDpsych(qq).bkg_orient,'orthogonal') && strcmp(CONDpsych(qq).bkg_contrast,'30')
                        size025OrthoCont30(subj) = 100*10.^QuestMean(CONDpsych(qq).q);    
                    elseif strcmp(CONDpsych(qq).bkg_size,'1') && strcmp(CONDpsych(qq).bkg_orient,'orthogonal') && strcmp(CONDpsych(qq).bkg_contrast,'10')
                        size1OrthoCont10(subj) = 100*10.^QuestMean(CONDpsych(qq).q);    
                    elseif strcmp(CONDpsych(qq).bkg_size,'1') && strcmp(CONDpsych(qq).bkg_orient,'orthogonal') && strcmp(CONDpsych(qq).bkg_contrast,'30')
                        size1OrthoCont30(subj) = 100*10.^QuestMean(CONDpsych(qq).q);    
                      
                    end
                end
        end
        doimage(gcf,fullfile(paths.Analysis,'psychometric'),'png',[subjIDs{subj},'psych_Result'],'300','painters',[],1) % this print the image as a pdf, might not work in a windows computer
    end

    % here we put participans data per condition in one table and save it
    if ismember(SITE,{'HYDERABAD','HAMBURG'})
        % removing imposible contrast values
        noBkg(noBkg>100)=100;           smallOrtho(smallOrtho>100)=100;
        fullOrtho(fullOrtho>100)=100;   smallParallel(smallParallel>100)=100;
        fullParallel(fullParallel>100)=100;
        contThresh = table(subjId',noBkg',smallParallel',fullParallel',smallOrtho',fullOrtho','VariableNames',{'ID','no Bkg','small ||','full ||','small #','full #'});
    elseif ismember(SITE,{'HAMBURGv2','HYDERABAD_controls','HYDERABAD_controls_blur'})
        noBkg(noBkg>100)=100;           size025ParallelCont10(size025ParallelCont10>100)=100;
        size025ParallelCont30(size025ParallelCont30>100)=100;   size1ParallelCont10(size1ParallelCont10>100)=100;
        size1ParallelCont30(size1ParallelCont30>100)=100;
        size025OrthoCont10(size025OrthoCont10>100)=100;
        size025OrthoCont30(size025OrthoCont30>100)=100;   size1OrthoCont10(size1OrthoCont10>100)=100;
        size1OrthoCont30(size1OrthoCont30>100)=100;
        contThresh = table(subjId',noBkg',size025ParallelCont10',size025ParallelCont30',size1ParallelCont10',size1ParallelCont30',size025OrthoCont10',size025OrthoCont30',size1OrthoCont10',size1OrthoCont30',...
            'VariableNames',{'ID','no Bkg','|| 0.25deg 10cont','|| 0.25deg 30cont','|| 1deg 10cont','|| 1deg 30cont',' # 0.25deg 10cont',' # 0.25deg 30cont',' # 1deg 10cont',' # 1deg 30cont'});
    end
    save(fullfile(paths.Analysis,'psychometric',['contThresh_Results_' SITE]),'contThresh')
end
%%
% Group analysis
% FIGURES FIRST VERSION EXPERIMENT
if ismember(SITE,{'HYDERABAD','HAMBURG'})
    % figures (many parameters for plotting are set in E354_paths_globalParams figs structure)
    load(fullfile(paths.Analysis,'psychometric',['contThresh_Results_' SITE]),'contThresh')
    colorAxis   = [6,2,4,1,3];
    xjitter     = (rand(size(contThresh,1),1)-.5)/2;
    figure;
    line(repmat([1:5]',1,size(contThresh,1))+xjitter',table2array(contThresh(:,2:end))','Color',[.9 .9 .9]),hold on
    for conds = 1:size(contThresh,2)-1
        scatter(conds+xjitter,table2array(contThresh(:,conds+1)),figs.scatter.mkSiz,'o','MarkerFaceColor',figs.cmapPaired(colorAxis(conds),:),...
            figs.scatter.options{:}),hold on
    end
    errorbar(1:size(contThresh,2)-1,mean(table2array(contThresh(:,2:end))),std(table2array(contThresh(:,2:end)))./sqrt(size(contThresh,1)),'Color',[0 0 0])
    axis([.5 5.5 -5 120])
    % singrank test for all possibl contrast, controlled for multiple
    % comparison
    allContrast = nchoosek(1:size(contThresh,2)-1,2);
    alphathr = .05/size(allContrast,1);
    a=0;
    for comps = 1:size(allContrast,1)
        p = signtest(contThresh.(allContrast(comps,1)+1),contThresh.(allContrast(comps,2)+1));
        if p<alphathr
             if p<.0001,tx = '***';bkg = 3;elseif p<.001,tx = '**';bkg = 2;elseif p<.01,tx = '*';bkg = 1;end
             line([allContrast(comps,1) mean(allContrast(comps,:))-.07*bkg; mean(allContrast(comps,:))+.07*bkg allContrast(comps,2)]',[120-a 120-a;120-a 120-a],'Color',[0 0 0],'LineWidth',.2)
           
            % tx = '*';
             text(mean(allContrast(comps,:)),120-a-1.5,tx,'HorizontalAlignment','Center','FontSize',5,'VerticalAlignment','middle')
        
            a = a+2;
        end
    end
    
    ylabel('contrast threshold','FontSize',figs.lblSiz)
    set(gca,'XTick',1:size(contThresh,2)-1,'XTickLabels',contThresh.Properties.VariableNames(2:end),'YTick',0:25:100,'FontSize',figs.txtSiz)
    
    doimage(gcf,fullfile(paths.Analysis,'psychometric'),'pdf',['contThreshAll_psych_Result'],'300','painters',figs.figSizBig/4,0)
    
    %
    %
    % figures log (many parameters for plotting are set in E354_paths_globalParams figs structure)
    xjitter     = (rand(size(contThresh,1),1)-.5)/2;
    figure;
    values = log10(table2array(contThresh(:,3:end))./repmat(table2array(contThresh(:,2)),1,4));
    line(repmat([1:4]',1,size(contThresh,1))+xjitter',values','Color',[.9 .9 .9]),hold on
    hline(0)
    for conds = 1:size(contThresh,2)-2
        scatter(conds+xjitter,values(:,conds),figs.scatter.mkSiz,'o','MarkerFaceColor',figs.cmapPaired(colorAxis(conds+1),:),...
            figs.scatter.options{:}),hold on
    end
    errorbar(1:size(contThresh,2)-2,mean(values),std(values)./sqrt(size(contThresh,1)),'Color',[0 0 0])
    axis([.5 4.5 -1 2.5])
    % singrank test for all possibl contrast, controlled for multiple
    % comparison
    allContrast = nchoosek(3:size(contThresh,2),2);
    alphathr = .05/size(allContrast,1);
    a=0;
    for comps = 1:size(allContrast,1)
        [h,p] = ttest(log10(contThresh.(allContrast(comps,1))./contThresh.(2)),log10(contThresh.(allContrast(comps,2))./contThresh.(2)))
        if p<alphathr
             if p<.0001,tx = '***';bkg = 3;elseif p<.001,tx = '**';bkg = 2;elseif p<.01,tx = '*';bkg = 1;end
             line([allContrast(comps,1) mean(allContrast(comps,:))-.07*bkg; mean(allContrast(comps,:))+.07*bkg allContrast(comps,2)]'-2,[2.5-a 2.5-a;2.5-a 2.5-a],'Color',[0 0 0],'LineWidth',.2)
    
            % tx = '*';
             text(mean(allContrast(comps,:))-2,2.5-a-.05,tx,'HorizontalAlignment','Center','FontSize',5,'VerticalAlignment','middle')
    
            a = a+.1;
        end
    end
    
    ylabel('thresh bkg./no bkg','FontSize',figs.lblSiz)
    tickvals = [.1 .5 1 5 10 50 100];
    set(gca,'XTick',1:size(contThresh,2)-2,'XTickLabels',contThresh.Properties.VariableNames(3:end),'YTick',log10(tickvals),'YTickLabels',tickvals,'FontSize',figs.txtSiz)
    
    doimage(gcf,fullfile(paths.Analysis,'psychometric'),'pdf',['contThreshAll_psych_ResultLog'],'300','painters',figs.figSizBig/4,0)
end

%%
% FIGURES SECOND VERSION EXPERIMENT
if ismember(SITE,{'HAMBURGv2','HYDERABAD_controls','HYDERABAD_controls_blur'})
    % figures (many parameters for plotting are set in E354_paths_globalParams figs structure)
    load(fullfile(paths.Analysis,'psychometric',['contThresh_Results_' SITE]),'contThresh')
    colorAxis   = [6,2,4,8,10,1,3,7,9];
    xjitter     = (rand(size(contThresh,1),1)-.5)/2;
    figure;
    line(repmat([1:9]',1,size(contThresh,1))+xjitter',table2array(contThresh(:,2:end))','Color',[.9 .9 .9]),hold on
    for conds = 1:size(contThresh,2)-1
        scatter(conds+xjitter,table2array(contThresh(:,conds+1)),figs.scatter.mkSiz,'o','MarkerFaceColor',figs.cmapPaired(colorAxis(conds),:),...
            figs.scatter.options{:}),hold on
    end
    errorbar(1:size(contThresh,2)-1,mean(table2array(contThresh(:,2:end))),std(table2array(contThresh(:,2:end)))./sqrt(size(contThresh,1)),'Color',[0 0 0])
    axis([.5 9.5 -5 120])
    % singrank test for all possibl contrast, controlled for multiple
    % comparison
    % allContrast = nchoosek(1:size(contThresh,2)-1,2);
    % alphathr = .05/size(allContrast,1);
    % a=0;
    % for comps = 1:size(allContrast,1)
    %     p = signtest(contThresh.(allContrast(comps,1)+1),contThresh.(allContrast(comps,2)+1));
    %     if p<alphathr
    %          if p<.0001,tx = '***';bkg = 3;elseif p<.001,tx = '**';bkg = 2;elseif p<.01,tx = '*';bkg = 1;end
    %          line([allContrast(comps,1) mean(allContrast(comps,:))-.07*bkg; mean(allContrast(comps,:))+.07*bkg allContrast(comps,2)]',[120-a 120-a;120-a 120-a],'Color',[0 0 0],'LineWidth',.2)
    % 
    %         % tx = '*';
    %          text(mean(allContrast(comps,:)),120-a-1.5,tx,'HorizontalAlignment','Center','FontSize',5,'VerticalAlignment','middle')
    % 
    %         a = a+2;
    %     end
    % end
    
    ylabel('contrast threshold','FontSize',figs.lblSiz)
    set(gca,'XTick',1:size(contThresh,2)-1,'XTickLabels',contThresh.Properties.VariableNames(2:end),'YTick',0:25:100,'FontSize',figs.txtSiz)
    
     doimage(gcf,fullfile(paths.Analysis,'psychometric'),'pdf',['contThreshAll_psych_Result'],'300','painters',figs.figSizBig./[2 4],1)
    
    %
    %
    % figures log (many parameters for plotting are set in E354_paths_globalParams figs structure)
    xjitter     = (rand(size(contThresh,1),1)-.5)/2;
    figure;
    values = log10(table2array(contThresh(:,3:end))./repmat(table2array(contThresh(:,2)),1,8));
    % for stats
    conthThreshStat = stack([contThresh(:,1) array2table(table2array(contThresh(:,3:end))./repmat(table2array(contThresh(:,2)),1,8),...
        'VariableNames',contThresh.Properties.VariableNames(3:end))],2:9,'NewDataVariableName','value');
conthThreshStat = [conthThreshStat splitvars(table(split(string(conthThreshStat.value_Indicator))),'Var1','NewVariableNames',{'Orientation','Size','Contrast'})];
lme = fitlme(conthThreshStat ,'value~Orientation*Size*Contrast+(1|ID)','DummyVarCoding','effects');

    line(repmat([1:8]',1,size(contThresh,1))+xjitter',values','Color',[.9 .9 .9]),hold on
    hline(0)
    for conds = 1:size(contThresh,2)-2
        scatter(conds+xjitter,values(:,conds),figs.scatter.mkSiz,'o','MarkerFaceColor',figs.cmapPaired(colorAxis(conds+1),:),...
            figs.scatter.options{:}),hold on
    end
    errorbar(1:size(contThresh,2)-2,mean(values),std(values)./sqrt(size(contThresh,1)),'Color',[0 0 0])
    axis([.5 8.5 -1 2.5])
    % % singrank test for all possibl contrast, controlled for multiple
    % % comparison
    % allContrast = nchoosek(3:size(contThresh,2),2);
    % alphathr = .05/size(allContrast,1);
    % a=0;
    % for comps = 1:size(allContrast,1)
    %     [h,p] = ttest(log10(contThresh.(allContrast(comps,1))./contThresh.(2)),log10(contThresh.(allContrast(comps,2))./contThresh.(2)))
    %     if p<alphathr
    %          if p<.0001,tx = '***';bkg = 3;elseif p<.001,tx = '**';bkg = 2;elseif p<.01,tx = '*';bkg = 1;end
    %          line([allContrast(comps,1) mean(allContrast(comps,:))-.07*bkg; mean(allContrast(comps,:))+.07*bkg allContrast(comps,2)]'-2,[2.5-a 2.5-a;2.5-a 2.5-a],'Color',[0 0 0],'LineWidth',.2)
    % 
    %         % tx = '*';
    %          text(mean(allContrast(comps,:))-2,2.5-a-.05,tx,'HorizontalAlignment','Center','FontSize',5,'VerticalAlignment','middle')
    % 
    %         a = a+.1;
    %     end
    % end
    
    ylabel('thresh bkg./no bkg','FontSize',figs.lblSiz)
    tickvals = [.1 .5 1 5 10 50 100];
    set(gca,'XTick',1:size(contThresh,2)-2,'XTickLabels',contThresh.Properties.VariableNames(3:end),'YTick',log10(tickvals),'YTickLabels',tickvals,'FontSize',figs.txtSiz)
    
    doimage(gcf,fullfile(paths.Analysis,'psychometric'),'pdf',['contThreshAll_psych_ResultLog'],'300','painters',figs.figSizBig./[2 4],1)

    % scatter ortho vs parallel
    %
     figure;
     for conds = 3:6
         thisBkg = contThresh.Properties.VariableNames{conds}(4:end);
         ixconds = find(contains(contThresh.Properties.VariableNames,thisBkg));
         scatter(values(:,ixconds(1)-2),values(:,ixconds(2)-2),figs.scatter.mkSiz,'o','MarkerFaceColor',figs.cmapPaired(colorAxis(conds-2+1),:),...
            figs.scatter.options{:}),hold on
     end
     for conds = 3:6
         thisBkg = contThresh.Properties.VariableNames{conds}(4:end);
         ixconds = find(contains(contThresh.Properties.VariableNames,thisBkg));
         errorbar(mean(values(:,ixconds(1)-2)),mean(values(:,ixconds(2)-2)),...
             std(values(:,ixconds(2)-2))./sqrt(size(values,1)),std(values(:,ixconds(2)-2))./sqrt(size(values,1)),...
             std(values(:,ixconds(1)-2))./sqrt(size(values,1)),std(values(:,ixconds(1)-2))./sqrt(size(values,1)),'LineWidth',2,'Color',figs.cmapPaired(colorAxis(conds-2+1),:))
         scatter(mean(values(:,ixconds(1)-2)),mean(values(:,ixconds(2)-2)),2*figs.scatter.mkSiz,'s','MarkerFaceColor',figs.cmapPaired(colorAxis(conds-2+1),:),...
            'LineWidth',.2,'MarkerEdgeColor',[0 0 0]),hold on
         text(-.4,2.3-conds/8,thisBkg,'Color',figs.cmapPaired(colorAxis(conds-2+1),:))
     end
     axis([-.5 2 -.5 2])
     line([-1 2.5]',[-1 2.5]','LineStyle',':','Color',[.5 .5 .5])
      line([0 0]',[-.5 2.5]','LineStyle',':','Color',[.5 .5 .5])
       line([-.5 2.5]',[0 0]','LineStyle',':','Color',[.5 .5 .5])
       text(1.6,0,'= no bkg.','FontSize',figs.lblSiz,'Color',[.5 .5 .5],'VerticalAlignment','top')
        text(0.1,1.3,'= no bkg.','FontSize',figs.lblSiz,'Color',[.5 .5 .5],'Rotation',-90)
     xlabel('parallel thresh.','FontSize',figs.lblSiz)
     ylabel('orthogonal thresh.','FontSize',figs.lblSiz)
     tickvals = [.1 .5 1 5 10 50 100];
     set(gca,'XTick',log10(tickvals),'XTickLabels',tickvals,'YTick',log10(tickvals),'YTickLabels',tickvals,'FontSize',figs.txtSiz)
       doimage(gcf,fullfile(paths.Analysis,'psychometric'),'pdf',['scatter_contThreshAll_psych_ResultLog'],'300','painters',figs.figSizBig./[2 2],1)


end