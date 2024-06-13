close all

load(fullfile(paths.Analysis,'SSVEPgroup','SSVEPtables'),'allNoneNone','resTables')

freqsToPlot     = [15 30 45 60];
measures        = {'psd','zscore'};
channAvgType    = {'grpChSel','indvChSel'};


cAx.psd        = [0 12;0 1;0 .6;0 .6];
cAx.zscore     = [0 151;0 40;0 50;0 40];
cAx.psdAvg     = [0 4;0 .3;0 .2;0 .1];
cAx.zscoreAvg  = [0 70;0 12;0 30;0 6];
cAx.psdR       = [-20 10;-20 10;-20 10;-20 10];
cAx.zscoreR    = [-5 5;-20 10;-20 10;-20 10];

if strcmp(SITE,'HAMBURG')
xpos        = [1 1.8 2.2 2.8 3.2 3.8 4.2 6 6.8 7.2 7.8 8.2 8.8 9.2];
xposR       = [.8 1.2 1.8 2.2 2.8 3.2 4.8 5.2 5.8 6.2 6.8 7.2];
contrasts   = [repmat({'50'},1,7) repmat({'100'},1,7)];
sizes       = repmat({'none','2','2','4','4','full','full'},1,2);
orients     = repmat({'none','parallel','orthogonal','parallel','orthogonal','parallel','orthogonal'},1,2);
colorAxis   = repmat([6,2,1,8,7,4,3],1,2);
colorAxisR  = repmat([2,1,8,7,4,3],1,2);
end

for freq =1%:length(freqsToPlot)
    mkdir(fullfile(paths.Analysis,'SSVEPgroup',['figures_SSVEP_' num2str(freqsToPlot(freq))]))
    thisTbl = resTables.([sprintf('table_%d',freqsToPlot(freq))]);
    for mes = 1%:length(measures)
        for chAvTyp = 1%:length(channAvgType)
            nsubjs      = length(unique(thisTbl.subject));
            xjitter     = (rand(nsubjs,1)-.5)/4;
            options_errorbar = {'-s',"MarkerSize",2,"MarkerEdgeColor",[0 0 0],'Color',[0 0 0],"MarkerFaceColor",[0.3 0.3 0.3]};
            % plots with averages and individual data
            fh = figure;, hold on
            thisFigureData = [];
            for conds = 1:length(xpos)
                thisFigureData = [thisFigureData,thisTbl(thisTbl.contrast==contrasts(conds) &...
                    thisTbl.bkgSize==sizes(conds) & thisTbl.bkgOrient==orients(conds),:).([measures{mes} '_' channAvgType{chAvTyp}])];
            end
            line(repmat(xpos(1:7)',1,size(xjitter,1))+xjitter',thisFigureData(:,1:7)','Color',[.9 .9 .9]),hold on
            line(repmat(xpos(8:14)',1,size(xjitter,1))+xjitter',thisFigureData(:,8:14)','Color',[.9 .9 .9]),hold on
            for conds = 1:length(xpos)

                scatter(xpos(conds)+xjitter,thisFigureData(:,conds),figs.scatter.mkSiz/2,'o','MarkerFaceColor',figs.cmapPaired(colorAxis(conds),:),...
                    figs.scatter.options{:}),hold on
            end
            for eb = [0,7]
                options_errorbar{end} = figs.cmapSet1(1,:);
                errorbar(xpos(1+eb),nanmean(thisFigureData(:,1+eb)),nanstd(thisFigureData(:,1+eb))./sqrt(nsubjs),options_errorbar{:})
                options_errorbar{end} = figs.cmapGray(1,:);
                errorbar(xpos([2,4,6]+eb),nanmean(thisFigureData(:,[2,4,6]+eb)),nanstd(thisFigureData(:,[2,4,6]+eb))./sqrt(nsubjs),options_errorbar{:})
                options_errorbar{end} = figs.cmapGray(2,:);
                errorbar(xpos([2,4,6]+1+eb),nanmean(thisFigureData(:,[2,4,6]+1+eb)),nanstd(thisFigureData(:,[2,4,6]+1+eb))./sqrt(nsubjs),options_errorbar{:})
            end
            set(gca,'XTick',[1:4,6:9],'XTickLabels',{'no Bkg','small','medium','full','no Bkg','small','medium','full'},'FontSize',figs.txtSiz)
            ylim(cAx.(measures{mes})(freq,:))
            text(2.5,cAx.(measures{mes})(freq,2)*.9,'50%','FontSize',figs.lblSiz,'HorizontalAlignment','center')
            text(5,cAx.(measures{mes})(freq,2)*.8,'|| parall.','Color',[.2 .2 .2],'FontSize',figs.lblSiz,'HorizontalAlignment','center')
            text(5,cAx.(measures{mes})(freq,2)*.7,'# ortho.','Color',[.6 .6 .6],'FontSize',figs.lblSiz,'HorizontalAlignment','center')
            text(7.5,cAx.(measures{mes})(freq,2)*.9,'100%','FontSize',figs.lblSiz,'HorizontalAlignment','center')
            ylabel(measures{mes} ,'FontSize',figs.lblSiz)
            doimage(fh,fullfile(paths.Analysis,'SSVEPgroup',['figures_SSVEP_' num2str(freqsToPlot(freq))]),'pdf',...
                sprintf('SSVEP%dHz_%s_%s',freqsToPlot(freq),channAvgType{chAvTyp},measures{mes}),...
                '300','painters',[figs.figSizBig(1).*.5 figs.figSizBig(1).*.25],1)


            % plots only with averages
            fh = figure; hold on

            for eb = [0,7]
                options_errorbar{end} = figs.cmapSet1(1,:);
                errorbar(xpos(1+eb),nanmean(thisFigureData(:,1+eb)),nanstd(thisFigureData(:,1+eb))./sqrt(nsubjs),options_errorbar{:});
                options_errorbar{end} = figs.cmapGray(1,:);
                errorbar(xpos([2,4,6]+eb),nanmean(thisFigureData(:,[2,4,6]+eb)),nanstd(thisFigureData(:,[2,4,6]+eb))./sqrt(nsubjs),options_errorbar{:});
                options_errorbar{end} = figs.cmapGray(2,:);
                errorbar(xpos([2,4,6]+1+eb),nanmean(thisFigureData(:,[2,4,6]+1+eb)),nanstd(thisFigureData(:,[2,4,6]+1+eb))./sqrt(nsubjs),options_errorbar{:});
            end
            set(gca,'XTick',[1:4,6:9],'XTickLabels',{'no Bkg','small','medium','full','no Bkg','small','medium','full'},'FontSize',figs.txtSiz)
            axis([0 10 cAx.([measures{mes} 'Avg'])(freq,:)])
            text(2.5,cAx.([measures{mes} 'Avg'])(freq,2)*.9,'50%','FontSize',figs.lblSiz,'HorizontalAlignment','center')
            text(5,cAx.([measures{mes} 'Avg'])(freq,2)*.8,'|| parall.','Color',[.2 .2 .2],'FontSize',figs.lblSiz,'HorizontalAlignment','center')
            text(5,cAx.([measures{mes} 'Avg'])(freq,2)*.7,'# ortho.','Color',[.6 .6 .6],'FontSize',figs.lblSiz,'HorizontalAlignment','center')
            text(7.5,cAx.([measures{mes} 'Avg'])(freq,2)*.9,'100%','FontSize',figs.lblSiz,'HorizontalAlignment','center')
            ylabel(measures{mes} ,'FontSize',figs.lblSiz)
            % "MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90])
            doimage(fh,fullfile(paths.Analysis,'SSVEPgroup',['figures_SSVEP_' num2str(freqsToPlot(freq))]),'pdf',...
                sprintf('SSVEP%dHzAvg_%s_%s',freqsToPlot(freq),channAvgType{chAvTyp},measures{mes}),...
                '300','painters',[figs.figSizBig(1).*.5 figs.figSizBig(1).*.25],1)





            % rarioplots with averages and individual data
            options_errorbar = {'-s',"MarkerSize",2,"MarkerEdgeColor",[0 0 0],'Color',[0 0 0],"MarkerFaceColor",[0.3 0.3 0.3]};
            if strcmp(measures{mes},'psd')
                fh = figure;, hold on

                thisFigureDataRatio = 10*log10([thisFigureData(:,2:7)./repmat(thisFigureData(:,1),1,6) thisFigureData(:,9:14)./repmat(thisFigureData(:,8),1,6)]);

                % mixed model
                participants = repmat([1:size(thisFigureDataRatio,1)],size(thisFigureDataRatio,2),1)';
                cont         = repmat(contrasts([2:7,9:14]),size(thisFigureDataRatio,1),1);
                siz          = repmat(sizes([2:7,9:14]),size(thisFigureDataRatio,1),1);
                ori          = repmat(orients([2:7,9:14]),size(thisFigureDataRatio,1),1);
                ratioStatTable = table(double(thisFigureDataRatio(:)),categorical(cont(:)),categorical(siz(:)), categorical(ori(:)),categorical(participants(:)),...
                    'VariableNames',{'ratio','FGcontrast','BKGsize','BKGorientation','participant'} );
                lme = fitlme(ratioStatTable,'ratio~FGcontrast*BKGsize*BKGorientation+(1|participant)','DummyVarCoding','effects')
                anova(lme)




                % line(repmat(xpos(1:7)',1,size(xjitter,1))+xjitter',thisFigureData(:,1:7)','Color',[.9 .9 .9]),hold on
                % line(repmat(xpos(8:14)',1,size(xjitter,1))+xjitter',thisFigureData(:,8:14)','Color',[.9 .9 .9]),hold on
                % axis([0 8 cAx.([measures{mes} 'R'])(freq,:)])
                hline(0)
                for conds = 1:length(xposR)

                    scatter(xposR(conds)+xjitter,thisFigureDataRatio(:,conds),figs.scatter.mkSiz/2,'o','MarkerFaceColor',figs.cmapPaired(colorAxisR(conds),:),...
                        figs.scatter.options{:}),hold on
                end
                for eb = [0,6]
                    options_errorbar{end} = figs.cmapGray(1,:);
                    errorbar(xposR([1,3,5]+eb),nanmean(thisFigureDataRatio(:,[1,3,5]+eb)),nanstd(thisFigureDataRatio(:,[1,3,5]+eb))./sqrt(nsubjs),options_errorbar{:})
                    options_errorbar{end} = figs.cmapGray(2,:);
                    errorbar(xposR([1,3,5]+1+eb),nanmean(thisFigureDataRatio(:,[1,3,5]+1+eb)),nanstd(thisFigureDataRatio(:,[1,3,5]+1+eb))./sqrt(nsubjs),options_errorbar{:})
                end
                set(gca,'XTick',[1:3,5:7],'XTickLabels',{'small','medium','full','small','medium','full'},'FontSize',figs.txtSiz)

                text(2,cAx.([measures{mes} 'R'])(freq,1)+range(cAx.([measures{mes} 'R'])(freq,:))*.95,'50%','FontSize',figs.lblSiz,'HorizontalAlignment','center')
                text(4,cAx.([measures{mes} 'R'])(freq,1)+range(cAx.([measures{mes} 'R'])(freq,:))*.8,'|| parall.','Color',[.2 .2 .2],'FontSize',figs.lblSiz,'HorizontalAlignment','center')
                text(4,cAx.([measures{mes} 'R'])(freq,1)+range(cAx.([measures{mes} 'R'])(freq,:))*.7,'# ortho.','Color',[.6 .6 .6],'FontSize',figs.lblSiz,'HorizontalAlignment','center')
                text(6,cAx.([measures{mes} 'R'])(freq,1)+range(cAx.([measures{mes} 'R'])(freq,:))*.95,'100%','FontSize',figs.lblSiz,'HorizontalAlignment','center')
                ylabel([measures{mes} ' bkg/noBkg. (dB)'],'FontSize',figs.lblSiz)
                doimage(fh,fullfile(paths.Analysis,'SSVEPgroup',['figures_SSVEP_' num2str(freqsToPlot(freq))]),'pdf',...
                    sprintf('SSVEPratio%dHz_%s_%s',freqsToPlot(freq),channAvgType{chAvTyp},measures{mes}),...
                    '300','painters',[figs.figSizBig(1).*.5 figs.figSizBig(1).*.25],1)
            end

            if strcmp(measures{mes},'psd') & freqsToPlot(freq) == 15
                fh = figure;, hold on
                thisFigureDataRationoOrient = [mean(thisFigureDataRatio(:,1:2),2) mean(thisFigureDataRatio(:,3:4),2) mean(thisFigureDataRatio(:,5:6),2) mean(thisFigureDataRatio(:,7:8),2) mean(thisFigureDataRatio(:,9:10),2) mean(thisFigureDataRatio(:,11:12),2)];
                participants2 = participants(:,1:2:end);
                cont2 = cont(:,1:2:end);
                siz2 = siz(:,1:2:end);
                ratioStatTable = table(double(thisFigureDataRationoOrient(:)),categorical(cont2(:)),categorical(siz2(:)),categorical(participants2(:)),...
                    'VariableNames',{'ratio','FGcontrast','BKGsize','participant'} );
                lme = fitlme(ratioStatTable,'ratio~FGcontrast*BKGsize+(1|participant)','DummyVarCoding','effects')
                anova(lme)

                errorbar(1:3,nanmean(thisFigureDataRationoOrient(:,1:3)),nanstd(thisFigureDataRationoOrient(:,1:3))./sqrt(nsubjs),options_errorbar{:});
                errorbar(1:3,nanmean(thisFigureDataRationoOrient(:,4:6)),nanstd(thisFigureDataRationoOrient(:,4:6))./sqrt(nsubjs),options_errorbar{:},'LineStyle','--');

                set(gca,'XTick',[1:3],'XTickLabels',{'small','medium','full'},'FontSize',figs.txtSiz)
                axis([0.5 3.5 -10 0])
                text(2.5,-1,'50%','FontSize',figs.lblSiz)
                text(2.5,-2,'100%','FontSize',figs.lblSiz)
                ylabel([measures{mes} ' bkg/noBkg. (dB)'],'FontSize',figs.lblSiz)
                % "MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90])
                doimage(fh,fullfile(paths.Analysis,'SSVEPgroup',['figures_SSVEP_' num2str(freqsToPlot(freq))]),'pdf',...
                    sprintf('SSVEP%dHzAvg_noOrient%s_%s',freqsToPlot(freq),channAvgType{chAvTyp},measures{mes}),...
                    '300','painters',[figs.figSizBig(1).*.25 figs.figSizBig(1).*.25],1)
                
                allContrast = [1 2;1 3;2 3;4 5;4 6;5 6;1 4;2 5;3 6];
                alphathr = .05/size(allContrast,1);
                a=0;
                for comps = 1:size(allContrast,1)
                    [h,p,ci,stats] = ttest(thisFigureDataRationoOrient(:,allContrast(comps,1)),thisFigureDataRationoOrient(:,allContrast(comps,2)));
                    sprintf('cont%s_siz%s vs cont%s_siz%s p = %1.4f',cont2{1,allContrast(comps,1)},siz2{1,allContrast(comps,1)},cont2{1,allContrast(comps,2)},siz2{1,allContrast(comps,2)},p)
                end
            end
        end
    end
end