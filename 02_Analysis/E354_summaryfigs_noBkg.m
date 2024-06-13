%%E354_summaryfigs_noBkg

% plotting
close all

freqsToPlot = [15 30 45 60];
nrows   = 5;           ncols = 6;
hmargin = .05;         vmargin = .05;

cAxpsd = [0 15;0 .5;0 .5;0 .5];
cAxdb  = [-10 10;-20 0;-20 0;-20 0];
cAxSNR = [0 30;0 10;0 10;0 10];
cAxZscore = [0 30;0 10;0 10;0 10];
mkdir(fullfile(paths.Analysis,'SSVEPgroup','figures_topos_each'))
for sameColorRange = [0:1]
    for freq = 1:length(freqsToPlot)
        thisFreq = freqsToPlot(freq);

        % psd
        fh = figure;
        fh.Position = [1 1 750 750];
        pp = 1;
        for subj = 1:size(allNoneNone.psd,3)
            subplot('Position',subplotFull(pp,-1,nrows,ncols,hmargin,vmargin,0))
            topoplot(allNoneNone.psd(FF==thisFreq,:,subj),chanlocs,'maplimits','maxmin','colormap',figs.cmapTopo1);  fixtopotlines
            if sameColorRange
                clim(cAxpsd(freq,:))
            else
                hc = colorbar;
                hc.Position = [hc.Position(1)+.02 hc.Position(2)+hc.Position(4)/6 hc.Position(3) hc.Position(4)*2/3];
            end
            text(0,.6,subjIDs{subj},'HorizontalAlignment','center')
            pp = pp+1;
        end
        if sameColorRange
           hc = colorbar;   
           hc.Position = [hc.Position(1)+.1 hc.Position(2) hc.Position(3)*1.5 hc.Position(4)];
           hc.Label.String = 'psd';
        end
        doimage(fh,fullfile(paths.Analysis,'SSVEPgroup','figures_topos_each'),'pdf',sprintf('toposAll_noBkg_psd_%dHz_sameRange%d',thisFreq,sameColorRange),'300','painters',figs.figSizBig,1)

         % dB
        fh = figure;
        fh.Position = [1 1 750 750];
        pp = 1;
        for subj = 1:size(allNoneNone.psd,3)
            subplot('Position',subplotFull(pp,-1,nrows,ncols,hmargin,vmargin,0))
            topoplot(10*log10(allNoneNone.psd(FF==thisFreq,:,subj)),chanlocs,'maplimits','maxmin','colormap',figs.cmapTopo1);  fixtopotlines
            if sameColorRange
               clim(cAxdb(freq,:))
            else
                hc = colorbar;
                hc.Position = [hc.Position(1)+.02 hc.Position(2)+hc.Position(4)/6 hc.Position(3) hc.Position(4)*2/3];
            end
            text(0,.6,subjIDs{subj},'HorizontalAlignment','center')
            pp = pp+1;
        end
        if sameColorRange
           hc = colorbar;   
           hc.Position = [hc.Position(1)+.1 hc.Position(2) hc.Position(3)*1.5 hc.Position(4)];
           hc.Label.String = 'dB';
        end
        doimage(fh,fullfile(paths.Analysis,'SSVEPgroup','figures_topos_each'),'pdf',sprintf('toposAll_noBkg_dB_%dHz_sameRange%d',thisFreq,sameColorRange),'300','painters',figs.figSizBig,1)

        % SNR
        fh = figure;
        fh.Position = [1 1 750 750];
        pp = 1;
        for subj = 1:size(allNoneNone.SNR,3)
            subplot('Position',subplotFull(pp,-1,nrows,ncols,hmargin,vmargin,0))
            topoplot(allNoneNone.SNR(FF==thisFreq,:,subj),chanlocs,'maplimits','maxmin','colormap',figs.cmapTopo1);  fixtopotlines
            if sameColorRange
                clim(cAxSNR(freq,:))
            else
                hc = colorbar;
                hc.Position = [hc.Position(1)+.02 hc.Position(2)+hc.Position(4)/6 hc.Position(3) hc.Position(4)*2/3];
            end
            text(0,.6,subjIDs{subj},'HorizontalAlignment','center')
            pp = pp+1;
        end
        if sameColorRange
           hc = colorbar;   
           hc.Position = [hc.Position(1)+.1 hc.Position(2) hc.Position(3)*1.5 hc.Position(4)];
           hc.Label.String = 'SNR';
        end
        doimage(fh,fullfile(paths.Analysis,'SSVEPgroup','figures_topos_each'),'pdf',sprintf('toposAll_noBkg_SNR_%dHz_sameRange%d',thisFreq,sameColorRange),'300','painters',figs.figSizBig,1)

        % Zscore
        fh = figure;
        fh.Position = [1 1 750 750];
        pp = 1;
        for subj = 1:size(allNoneNone.Zscore,3)
            subplot('Position',subplotFull(pp,-1,nrows,ncols,hmargin,vmargin,0))
            topoplot(allNoneNone.Zscore(FF==thisFreq,:,subj),chanlocs,...
                'emarker2',{find(allNoneNone.Zscore(FF==thisFreq,:,subj)>anlys.zscoreThresh),'.','k',8,1},...  
                'maplimits','maxmin','colormap',figs.cmapTopo1);  fixtopotlines
            if sameColorRange
                clim(cAxZscore(freq,:))
            else
                hc = colorbar;
                hc.Position = [hc.Position(1)+.02 hc.Position(2)+hc.Position(4)/6 hc.Position(3) hc.Position(4)*2/3];
            end
            text(0,.6,subjIDs{subj},'HorizontalAlignment','center')
            pp = pp+1;
        end
        text(1,0.4,sprintf('* Zscore > %1.2f',anlys.zscoreThresh))
        if sameColorRange
           hc = colorbar;   
           hc.Position = [hc.Position(1)+.1 hc.Position(2) hc.Position(3)*1.5 hc.Position(4)];
           hc.Label.String = 'Zscore';
        end
        
        doimage(fh,fullfile(paths.Analysis,'SSVEPgroup','figures_topos_each'),'pdf',sprintf('toposAll_noBkg_Zscore_%dHz_sameRange%d',thisFreq,sameColorRange),'300','painters',figs.figSizBig,1)

    end
end

%%
fh = figure;
fh.Position = [68 1 600 600];
freqsToPlot = [15 30 45 60];
nrows       = length(freqsToPlot);
ncols       = 4; hmargin = .05; vmargin = .05;

for freq = 1:length(freqsToPlot)
    thisFreq = freqsToPlot(freq);
    
    % median psd
    subplot('Position',subplotFull(freq,1,nrows,ncols,hmargin,vmargin))
    % subplot(nrows,ncols,ncols*(freq-1)+1)
    topoplot(median(allNoneNone.psd(FF==thisFreq,:,:),3),chanlocs,'maplimits','maxmin','colormap',figs.cmapTopo1);  fixtopotlines%'maplimits',[-axLimSNR(4) axLimSNR(4)]
    hc = colorbar; hc.Limits = [0 hc.Limits(2)]; hc.Ticks = floor(hc.Limits*100)/100; hc.FontSize = 6;
    hc.Position = [hc.Position(1)+.03 hc.Position(2)+hc.Position(4)/6 hc.Position(3) hc.Position(4)*2/3];
    text(-.75,0,sprintf('%d\nHz',thisFreq),'FontWeight','bold','FontSize',9)
    if freq==1,text(0,.65,'median(psd)','HorizontalAlignment','center','FontSize',10),end
    
    % average log(psd)
    subplot('Position',subplotFull(freq,2,nrows,ncols,hmargin,vmargin))
    topoplot(mean(10*log10(allNoneNone.psd(FF==thisFreq,:,:)),3),chanlocs,'maplimits','maxmin','colormap',figs.cmapTopo1);  fixtopotlines%'maplimits',[-axLimSNR(4) axLimSNR(4)]
    hc = colorbar; hc.Ticks = [ceil(hc.Limits(1)) floor(hc.Limits(2))]; hc.FontSize = 6;
    hc.Position = [hc.Position(1)+.03 hc.Position(2)+hc.Position(4)/6 hc.Position(3) hc.Position(4)*2/3];
    if freq==1,text(0,.65,'avg(psd) (dB)','HorizontalAlignment','center','FontSize',10),end
   
    % median (SNR)
    subplot('Position',subplotFull(freq,3,nrows,ncols,hmargin,vmargin))
    topoplot(median(allNoneNone.SNR(FF==thisFreq,:,:),3),chanlocs,'maplimits','maxmin','colormap',figs.cmapTopo1);  fixtopotlines%'maplimits',[-axLimSNR(4) axLimSNR(4)]
    hc = colorbar; hc.Limits = [0 hc.Limits(2)]; hc.Ticks = floor(hc.Limits*10)/10; hc.FontSize = 6;
    hc.Position = [hc.Position(1)+.03 hc.Position(2)+hc.Position(4)/6 hc.Position(3) hc.Position(4)*2/3];
    if freq==1,text(0,.65,'median(SNR)','HorizontalAlignment','center','FontSize',10),end

    % median (Zscore)
    subplot('Position',subplotFull(freq,4,nrows,ncols,hmargin,vmargin))
    topoplot(median(allNoneNone.Zscore(FF==thisFreq,:,:),3),chanlocs,...
         'emarker2',{find(median(allNoneNone.Zscore(FF==thisFreq,:,:),3)>anlys.zscoreThresh),'*',[.99 .99 .99],3,1} ,...
         'maplimits','maxmin','colormap',figs.cmapTopo1);  fixtopotlines
    hc = colorbar; hc.Limits = [0 hc.Limits(2)]; hc.Ticks = floor(hc.Limits*10)/10; hc.FontSize = 6;
    hc.Position = [hc.Position(1)+.03 hc.Position(2)+hc.Position(4)/6 hc.Position(3) hc.Position(4)*2/3];
    if freq==1,
    text(0,.65,'median(Zscore)','HorizontalAlignment','center','FontSize',10),
    text(0,-.65,sprintf('* median(Zscore) > %1.2f',anlys.zscoreThresh),'HorizontalAlignment','center','FontSize',6),
    end

end
doimage(fh,fullfile(paths.Analysis,'SSVEPgroup'),'pdf',sprintf('toposAllavgs_noBkg'),'300','painters',figs.figSizBig.*.75,1)
         