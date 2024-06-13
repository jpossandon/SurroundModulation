load(fullfile(paths.Analysis,'SSVEPgroup','SSVEPtables'),'resTables')
load(fullfile(paths.Analysis,'psychometric','contThresh_Results'))

psychCols = {'no Bkg', 'small ||',  'full ||',  'small #',  'full #'};
SSVEPRows = {{'none','none'}, {'2','parallel'}, {'full','parallel'} , {'2','orthogonal'}, {'full','orthogonal'}};
fgCont    = categorical([50,100]);

%%
  for cont = 1:2
      thisContallSSVEP=[];
      thisContallbehav=[];
for ncorrs = 2:length(psychCols)
    thisBehavValues =log10(table2array(contThresh(:,psychCols{ncorrs}))./table2array(contThresh(:,psychCols{1})));

  
        thisSSVEPTable = resTables.table_15(resTables.table_15.bkgSize==SSVEPRows{ncorrs}{1} & ...
            resTables.table_15.bkgOrient==SSVEPRows{ncorrs}{2} & resTables.table_15.contrast == fgCont(cont),1:5);
        noneSSVEPTable  = resTables.table_15(resTables.table_15.bkgSize=='none' & resTables.table_15.bkgOrient=='none' & resTables.table_15.contrast == fgCont(cont),1:5);
        
        thisSSVEPValues = [];
        if all(thisSSVEPTable.subject == noneSVEPTable.subject)
            thisSSVEPValues = 10*log10(thisSSVEPTable.psd_grpChSel./noneSSVEPTable.psd_grpChSel);
        
            if categorical(contThresh.ID)==thisSSVEPTable.subject % check participants are the same
                thisContallSSVEP=[thisContallSSVEP,thisSSVEPValues];
      thisContallbehav=[thisContallbehav,thisBehavValues]; 
                
                fh = figure;, hold on
                  [r,p] = corr(thisSSVEPValues,thisBehavValues);
                scatter(thisSSVEPValues,thisBehavValues,figs.scatter.mkSiz,'o','MarkerFaceColor',figs.cmapSet1(ncorrs,:),...
                    figs.scatter.options{:}),hold on
                text(-18,-.8,sprintf('r: %1.2f, p = %1.3f',r,p),'FontSize',figs.txtSiz)
                axis([-20 10 -1 2.5])
                tickvals = [.1 .5 1 5 10 50 100];
               
            set(gca,'YTick',log10(tickvals),'YTickLabels',tickvals,'FontSize',figs.txtSiz)

            xlabel([psychCols{ncorrs} '/no Bkg. EEG'] ,'FontSize',figs.lblSiz)
            ylabel([psychCols{ncorrs} '/no Bkg. behav'] ,'FontSize',figs.lblSiz)
            
            mkdir(fullfile(paths.Analysis,'SSVEPgroup','SSVEP_behav_corr'))
            doimage(fh,fullfile(paths.Analysis,'SSVEPgroup','SSVEP_behav_corr'),'pdf',...
                sprintf('corr_FGcont_%s_BkgOrient_%s_BkgSize_%s',fgCont(cont),SSVEPRows{ncorrs}{2},SSVEPRows{ncorrs}{1}),...
                 '300','painters',figs.figSizBig/4,1)
            end
        end
end

thisContallSSVEP = mean(thisContallSSVEP,2);
thisContallbehav =mean(thisContallbehav,2);
[r,p] = corr(thisContallSSVEP,thisContallbehav )
end