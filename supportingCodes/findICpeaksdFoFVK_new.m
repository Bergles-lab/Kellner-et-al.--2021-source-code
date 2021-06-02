%Kellner et al., 2021 
% function [wholeROIinfo] = findICpeaksdFoFVK(ICsignal,plotFlag,cellAns,thresh)
function [wholeROIinfo,pkData] = findICpeaksdFoFVK_new(ICsignal,plotFlag,cellAns)

% [wholeROIinfo] = findICpeaksdFoFVK(ICsignal,plotFlag,cellAns,LICvar2,RICvar2)
% if ~exist('analysisName','var')
%     analysisName = 'stock';
% end
% if ~exist('plotFlag','var')
%     plotFlag = 1;
% end

time = [1:1:size(ICsignal,1)]';
LIC = smooth(ICsignal(:,1));
RIC = smooth(ICsignal(:,2));
ctx = smooth(ICsignal(:,3));

% if LICvar2<0.001
%     LIC=LIC-ctx;
% end
% if RICvar2<0.001
%     RIC=RIC-ctx;
% end
% LIC=LIC-ctx;
% RIC=RIC-ctx;
%parameters
% if strcmp(cellAns,'Astrocyte')
%     pkThreshold = 0.005;
%     pkMinHeight = 0.002;
% else
%     pkThreshold = .01;
%     pkMinHeight = .01;
% end
pkThreshold=0.01;
% pkMinHeight=thresh;
% pkDistance = 10; %in frame, 10 = 1s %20
[pks,locs,w] = findpeaks(LIC,'MinPeakHeight',pkThreshold,'MinPeakProminence',pkThreshold/2,...
    'WidthReference','halfheight','MinPeakWidth',5,'Annotate','extents');

% [pks,locs,w] = findpeaks(LIC,'MinPeakProminence',pkThreshold,'MinPeakHeight',pkThreshold,'MinPeakDistance',pkDistance,'WidthReference','halfheight','Annotate','extents');
[LIC_itk,LIC_rmv] = ctx_PkRemoval(ctx, LIC, locs);
pksLIC = pks(LIC_itk);
locsLIC = locs(LIC_itk);
% locsLICRmv = locs(LIC_rmv);
wLIC = w(LIC_itk);
% if median(LIC)<0
%     pks=pks-median(LIC);% add baseline to amplitude
% else
%     pks=pks+median(LIC);
% end
% LICinfo = [pks locs w];
[pksR,locsR,wR] = findpeaks(RIC,'MinPeakProminence',pkThreshold/2,'MinPeakHeight',pkThreshold,'MinPeakWidth',5,'WidthReference','halfheight');
% [pksR,locsR,wR] = findpeaks(RIC,'MinPeakProminence',pkThreshold,'MinPeakHeight',pkThreshold,'MinPeakDistance',pkDistance,'WidthReference','halfheight','Annotate','extents');
[RIC_itk, RIC_rmv] = ctx_PkRemoval(ctx, RIC, locsR);
pksRIC = pksR(RIC_itk);
locsRIC = locsR(RIC_itk);
% locsRICRmv = locsR(RIC_rmv);
wRIC = wR(RIC_itk);
% if median(RIC)<0
%     pksR=pksR-median(RIC);% add baseline to amplitude
% else
%     pksR=pksR+median(RIC);
% end
% RICinfo = [pksR locsR wR];

% % if the "cortex" response is larger than either the left or right peak -
% % remove the peak from both left and right IC
% for i=1:length(locsLICRmv)
%     if sum(ismember(locsRIC,[locsLICRmv(i)-5:locsLICRmv(i)+5]))>=1
%         ind=find(ismember(locsRIC,[locsLICRmv(i)-5:locsLICRmv(i)+5]));
%         locsRIC(ind)=[];
%         pksRIC(ind)=[];
%         wRIC(ind)=[];
%     end
% end
% for ii=1:length(locsRICRmv)
%     if sum(ismember(locsLIC,[locsRICRmv(ii)-5:locsRICRmv(ii)+5]))>=1
%         ind=find(ismember(locsLIC,[locsRICRmv(ii)-5:locsRICRmv(ii)+5]));
%         locsLIC(ind)=[];
%         pksLIC(ind)=[];
%         wLIC(ind)=[];
%     end
% end
pkData = ICcompare2(LIC, RIC, pksLIC,locsLIC,pksRIC,locsRIC); %get data about dominance
%col 1: LIC loc, 2: LIC pk, 3:RIC loc, 4:RIC pk, 5:delta, 6: peak type (1=matched, 2=LIC only, 3=RIC only) 7: which is bigger (1=LIC, 2=RIC)
LICinfo = [pksLIC locsLIC wLIC];
RICinfo = [pksRIC locsRIC wRIC];
%     if ~isempty(pksRIC)
%     pkData = ICcompare2(LIC, RIC, pksLIC,locsLIC,pksRIC,locsRIC);
%     end

%% plot the average event for left and right
if strcmp(cellAns,'Astrocyte')
    tEvent=4; %4 seconds
else
    tEvent=2;
end
sampRate=10; %in hz
% tEvent=2;
% t=1:tEvent*sampRate*2+1;
if ~isempty(LICinfo)
    LICeventMat=nan(length(LICinfo(:,2)),(tEvent*2*sampRate)+1);
    for a=1:length(LICinfo(:,2))
        if LICinfo(a,2)-(tEvent*sampRate)>0 && (LICinfo(a,2)+(tEvent*sampRate))<=length(LIC)
            LICeventMat(a,:)=LIC(LICinfo(a,2)-(tEvent*sampRate):LICinfo(a,2)+(tEvent*sampRate)); %taking 40 frames as arbitrary
        elseif LICinfo(a,2)-(tEvent*sampRate)<=0 && (LICinfo(a,2)+(tEvent*sampRate))<=length(LIC)
            LICeventMat(a,:)=[nan(abs(LICinfo(a,2)-(tEvent*sampRate)),1);(LIC(1:LICinfo(a,2)+(tEvent*sampRate)+1))];
        elseif LICinfo(a,2)-(tEvent*sampRate)>0 && (LICinfo(a,2)+(tEvent*sampRate))>length(LIC)
            LICeventMat(a,:)=[LIC(LICinfo(a,2)-(tEvent*sampRate):end);nan((LICinfo(a,2)+(tEvent*sampRate))-length(LIC),1)];
        end
        
    end
else
    LICeventMat=[];
end
if ~isempty(RICinfo)
    RICeventMat=nan(length(RICinfo(:,2)),(tEvent*2*sampRate)+1);
    for a=1:length(RICinfo(:,2))
        if RICinfo(a,2)-(tEvent*sampRate)>0 && (RICinfo(a,2)+(tEvent*sampRate))<=length(RIC)
            RICeventMat(a,:)=RIC(RICinfo(a,2)-(tEvent*sampRate):RICinfo(a,2)+(tEvent*sampRate)); %taking 40 frames as arbitrary
        elseif RICinfo(a,2)-(tEvent*sampRate)<=0 && (RICinfo(a,2)+(tEvent*sampRate))<=length(RIC)
            RICeventMat(a,:)=[nan(abs(RICinfo(a,2)-(tEvent*sampRate)),1);(RIC(1:RICinfo(a,2)+(tEvent*sampRate)+1))];
        elseif RICinfo(a,2)-(tEvent*sampRate)>0 && (RICinfo(a,2)+(tEvent*sampRate))>length(RIC)
            RICeventMat(a,:)=[RIC(RICinfo(a,2)-(tEvent*sampRate):end);nan((RICinfo(a,2)+(tEvent*sampRate))-length(RIC),1)];
        end
    end
else
    RICeventMat=[];
end
if ~isempty(LICinfo) && ~isempty(RICinfo)
    eventTime=[-(tEvent*sampRate):(tEvent*sampRate)];
    figure; h1=shadedErrorBar(eventTime,nanmean(LICeventMat,1),nanstd(LICeventMat,[],1)/sqrt(size(LICeventMat,1)),'lineprops','-b','transparent',1);
    hold on; h2=shadedErrorBar(eventTime,nanmean(RICeventMat,1),nanstd(RICeventMat,[],1)/sqrt(size(RICeventMat,1)),'lineprops','-r','transparent',1);
    legend([h1.mainLine,h2.mainLine],'Left IC','Right IC')
%     [fitobject,gof,output] = fit(eventTime',nanmean(LICeventMat,1)','gauss1');
%     plot(fitobject,'k');title(num2str(gof.adjrsquare))
elseif ~isempty(LICinfo)
    eventTime=[-(tEvent*sampRate):(tEvent*sampRate)];
    figure; h1=shadedErrorBar(eventTime,nanmean(LICeventMat,1),nanstd(LICeventMat,[],1)/sqrt(size(LICeventMat,1)),'lineprops','-b','transparent',1);
    legend([h1.mainLine],'Left IC')
elseif ~isempty(RICinfo)
    eventTime=[-(tEvent*sampRate):(tEvent*sampRate)];
    figure; h1=shadedErrorBar(eventTime,nanmean(RICeventMat,1),nanstd(RICeventMat,[],1)/sqrt(size(RICeventMat,1)),'lineprops','-r','transparent',1);
    legend([h1.mainLine],'Right IC')
else
    figure;
end

%% plotting
if(plotFlag)
    %% plot signal with peaks
    lt_org = [255, 166 , 38]/255;
    dk_org = [255, 120, 0]/255;
    lt_blue = [50, 175, 242]/255;
    dk_blue = [0, 13, 242]/255;
    figure('Position',[50 100 1200 500])
    %     %     if median(LIC)<0
    %     %         plot(time/10,LIC-median(LIC),'Color',dk_org);
    %     %     else
    %     %         plot(time/10,LIC+median(LIC),'Color',dk_org);
    %     %     end
    %     %     hold on;
    %     %     if median(LIC)<0
    %     %         plot(time/10,RIC-median(RIC),'Color',dk_blue);
    %     %     else
    %     %         plot(time/10,RIC+median(RIC),'Color',dk_blue);
    %     %     end
    %     sampRate=10;
    %     plot(time/sampRate,LIC,'Color',lt_org,'LineWidth',1.5);
    %     hold on;
    %     plot(time/sampRate,RIC,'Color',lt_blue,'LineWidth',1.5);
    %     plot(time/sampRate,ctx,'y');
    %     if ~isempty(RICinfo)
    %         plot(RICinfo(:,2)/sampRate,RICinfo(:,1),'LineStyle','none','Marker','o','MarkerFaceColor',lt_blue,'MarkerEdgeColor',lt_blue,'MarkerSize',8);
    %     end
    %     if ~isempty(LICinfo)
    %         plot(LICinfo(:,2)/sampRate,LICinfo(:,1),'LineStyle','none','Marker','o','MarkerFaceColor',lt_org,'MarkerEdgeColor',lt_org,'MarkerSize',8);
    %     end
    %     ylabel('Amplitude (dF/F)')
    %     xlabel('Time (s)');
    %     xlim([0 ceil(length(LIC)/sampRate)])
    %     legend('Left','Right','Ctx')
    sampRate=10;
    plot(time,LIC,'Color',lt_org,'LineWidth',1.5);
    hold on;
    plot(time,RIC,'Color',lt_blue,'LineWidth',1.5);
    plot(time,ctx,'y');
    if ~isempty(RICinfo)
        plot(RICinfo(:,2),RICinfo(:,1),'LineStyle','none','Marker','o','MarkerFaceColor',lt_blue,'MarkerEdgeColor',lt_blue,'MarkerSize',8);
    end
    if ~isempty(LICinfo)
        plot(LICinfo(:,2),LICinfo(:,1),'LineStyle','none','Marker','o','MarkerFaceColor',lt_org,'MarkerEdgeColor',lt_org,'MarkerSize',8);
    end
    ylabel('Amplitude (dF/F)')
    xlabel('Time (frames)');
    xlim([0 ceil(length(LIC))])
    legend('Left','Right','Ctx')
    
    %     legend('Left','Right')
    %     %     graphEvents(pkData,[filePath,'timevert_',analysisName],length(LIC),cellAns);
    %% lollipop plots
    h=figure;
    ax = gca;
    %     %     if ~isempty(LICinfo)
    %     %     hline1=line([LICinfo(:,2)'; LICinfo(:,2)'],[(LICinfo(:,1))'; zeros(size(LICinfo,1),1)'],'Color',lt_org,'LineWidth',1.25);
    %     %     end
    %     %     hold on
    %     %     if ~isempty(RICinfo)
    %     %     h2=line([RICinfo(:,2)'; RICinfo(:,2)'],[(RICinfo(:,1))'; zeros(size(RICinfo,1),1)'],'Color',lt_blue,'LineWidth',1.25);
    %     %     end
    %     %     xlabel('Frame number')
    %     %     ylabel('dF/F amplitude')
    %     %     xlim([0 (length(LIC))])
    %     %     if ~isempty(LICinfo) && ~isempty(RICinfo)
    %     %     legend([hline1(1), h2(1)],'Left IC','Right IC')
    %     %     end
    %
    %     if ~isempty(LICinfo)
    %         line([(-LICinfo(:,1))'; zeros(size(LICinfo,1),1)'],[-LICinfo(:,2)'; -LICinfo(:,2)'],'Color',lt_org,'LineWidth',1.5);
    %         hold on
    %         scatter(-pkData((pkData(:,7)==1),2),-pkData((pkData(:,7)==1),1),pkData(pkData(:,7)==1,5)*200,dk_org,'Marker','o','MarkerFaceColor', dk_org);
    %     end
    %     hold on;
    %     if ~isempty(RICinfo)
    %         line([zeros(size(RICinfo,1),1)'; (RICinfo(:,1))'; ],[-RICinfo(:,2)'; -RICinfo(:,2)'],'Color',lt_blue,'LineWidth',1.5);
    %         scatter(pkData((pkData(:,7)==2),4),-pkData((pkData(:,7)==2),3),pkData(pkData(:,7)==2,5)*200,dk_blue,'Marker','o','MarkerFaceColor', dk_blue);
    %     end
    if ~isempty(pkData)
        line([(-pkData(:,2))'; zeros(size(pkData,1),1)'],[-pkData(:,1)'; -pkData(:,1)'],'Color',lt_org,'LineWidth',2);
        hold on;
        line([zeros(size(pkData,1),1)'; (pkData(:,4))'; ],[-pkData(:,3)'; -pkData(:,3)'],'Color',lt_blue,'LineWidth',2);
        
        scatter(pkData((pkData(:,7)==2),4),-pkData((pkData(:,7)==2),3),pkData(pkData(:,7)==2,5)*800,dk_blue,'MarkerFaceColor', dk_blue);
        scatter(-pkData((pkData(:,7)==1),2),-pkData((pkData(:,7)==1),1),pkData(pkData(:,7)==1,5)*800,dk_org,'MarkerFaceColor', dk_org);
    end
    line([0 0]',[-length(LIC) 50]','Color','black');
    %     scatter(pkData((pkData(:,7)==2),4),-pkData((pkData(:,7)==2),3),pkData(pkData(:,7)==2,5)*200,dk_blue,'filled','MarkerEdgeColor',dk_blue); %'MarkerEdgeColor',dk_blue);
    %     scatter(-pkData((pkData(:,7)==1),2),-pkData((pkData(:,7)==1),1),pkData(pkData(:,7)==1,5)*200,dk_org,'filled','MarkerEdgeColor',dk_org); %,'MarkerEdgeColor',dk_org);
    set(h,'Position',[500,300,350,500]);
    if strcmp(cellAns,'Astrocyte')
        ax.XLim = [-0.05 0.05]; %.5
    else
        ax.XLim = [-0.5 0.5];
    end
    ax.YLim = [-length(LIC) 50];
    ylabel('Frames (reversed)')
    xlabel('dF/F')
    
    %     if ~isempty(LICinfo) && ~isempty(RICinfo)
    %         legend([hline1(1), h2(1)],'Left IC','Right IC')
    %     end
    %     %     xticks=get(gca,'XTick')/10;
    %     %     set(gca,'XTick',[xticks(1):xticks(end)])
%     %% plot histograms for amplitude and hw
%         if ~isempty(RICinfo)
%             Rpks = RICinfo(:,1);
%         else
%             Rpks=[];
%         end
%         if ~isempty(LICinfo)
%             Lpks = LICinfo(:,1);
%         else
%             Lpks=[];
%         end
%         hold off;
%         if strcmp(cellAns,'Astrocyte')
% %             [Rcounts, edges]= histcounts(Rpks,[0:.01:.42]); %.42 %0.1
% %             [Lcounts,edges] =histcounts(Lpks,[0:.01:.42]); %[0:.06:.42 100]
%             [Rcounts, edges]= histcounts(Rpks,[0:.005:.05]); %.42 %0.1
%             [Lcounts,edges] =histcounts(Lpks,[0:.005:.05]); %[0:.06:.42 100]
%         else
%             [Rcounts, edges]= histcounts(Rpks,[0:.01:.42]); %.42 %.03
%             [Lcounts,edges] =histcounts(Lpks,[0:.01:.42]); %[0:.06:.42 100]
%         end
%     
%         %         binY = [.03:.06:.45];
%         figure;
%         h=barh(edges(1:end-1),Rcounts,.9);
%         hold on;
%         barh(edges(1:end-1),-Lcounts,.9,'FaceColor',lt_org,'EdgeColor','none');
%         %barh(Lbins,-Lcounts,.9,'FaceColor',lt_org,'EdgeColor',lt_org);
%         h.FaceColor = lt_blue;
%         h.EdgeColor = 'none';
%         %xlim([-100 100]);
%         xlim([-20 20]);
%     %     xticks([-10 0 10]);
%         %xticks([-100 -50 0 50 100]);
%         if strcmp(cellAns,'Astrocyte')
%             ylim([0 .30]);%.5 %0.1
%         else
%             ylim([0 .40])
%         end
%         %             yticks([0 .25 .50]);
%         box off;
%         ylabel('dF/F peak amplitude')
%         xlabel('Number of events')
%     
    
    %% plot histograms of event widths and amplitudes - fix
    if ~isempty(LICinfo) && ~isempty(RICinfo)
        Amps=[RICinfo(:,1);LICinfo(:,1)];
        Wids=[RICinfo(:,3);LICinfo(:,3)];
    elseif ~isempty(RICinfo)
        Amps=RICinfo(:,1);
        Wids=RICinfo(:,3);
    elseif ~isempty(LICinfo)
        Amps=LICinfo(:,1);
        Wids=LICinfo(:,3);
    else
        Amps=[];
        Wids=[];
    end
    hold off;
    if strcmp(cellAns,'Astrocyte')
        [counts, edges]= histcounts(Amps,[0:.01:.42]); %.42 %0.1
    else
        [counts, edges]= histcounts(Amps,[0:.01:.42]); %.42
    end
    
    figure;
    h=bar(edges(1:end-1),counts,.9);
    h.FaceColor = lt_blue;
    h.EdgeColor = 'none';
    
    if strcmp(cellAns,'Astrocyte')
        xlim([0 .30]);%.5 %0.1
    else
        xlim([0 .50])
    end
    box off;
    xlabel('dF/F peak amplitude')
    ylabel('Number of events')
    title('Whole ROI')
    
    [Wcounts,Wedges]=histcounts(Wids/sampRate,[0:0.3:10]);
    figure;
    h=bar(Wedges(1:end-1),Wcounts,.9);
    h.FaceColor = lt_blue;
    h.EdgeColor = 'none';
    xlim([0 10])
    box off;
    xlabel('Half width of detected peaks (s)')
    ylabel('Event count')
    title('Whole ROI')
    %% put info into structure
    wholeROIinfo=struct;
    wholeROIinfo.LIC=LICinfo;
    wholeROIinfo.RIC=RICinfo;
    wholeROIinfo.frmNum=length(time);
    if ~isempty(LICeventMat)
        wholeROIinfo.avgEventL=nanmean(LICeventMat);
        wholeROIinfo.stdEventL=nanstd(LICeventMat,[],1);
        wholeROIinfo.nEventL=size(LICinfo,1);
    else
        wholeROIinfo.avgEventL=[];
        wholeROIinfo.stdEventL=[];
        wholeROIinfo.nEventL=[];
    end
    if ~isempty(RICeventMat)
        wholeROIinfo.avgEventR=nanmean(RICeventMat);
        wholeROIinfo.stdEventR=nanstd(RICeventMat,[],1);
        wholeROIinfo.nEventR=size(RICinfo,1);
    else
        wholeROIinfo.avgEventR=[];
        wholeROIinfo.stdEventR=[];
        wholeROIinfo.nEventR=[];
    end
    
end
end

function [indexToKeep,indxToRmv] = ctx_PkRemoval(ctxSignal, ICsignal, IClocs)
% time = [1:1:size(ctxSignal,1)]';

%[pks,locs,w] = findpeaks(msbackadj(time,ctxSignal),'MinPeakHeight',0.03,'Annotate','extents');

ctxBright = ctxSignal > ICsignal;

indexToKeep = [];
indxToRmv = [];
for i=1:size(IClocs)
    if sum(ctxBright(IClocs(i)-1:IClocs(i)+1)) == 0
        indexToKeep = [indexToKeep i];
    else
        indxToRmv=[indxToRmv i];
    end
end
end


function [pkData] = ICcompare2(LIC, RIC, master_pks,master_locs,pks,locs)
windowforhit = 10;

pkData = []; %col 1: LIC loc, 2: LIC pk, 3:RIC loc, 4:RIC pk, 5:delta, 6: peak type (1=matched, 2=LIC only, 3=RIC only) 7: which is bigger (1=LIC, 2=RIC)
locs = [locs(:) zeros(size(locs,1),1)];
for i=1:size(master_pks,1)
    match = 0;
    j = 1;
    if ~isempty(locs)
        while match == 0
            if abs(master_locs(i)-locs(j,1)) <= windowforhit && locs(j,2) ~= 1
                pkData(i, 1:6) = [master_locs(i) master_pks(i) locs(j,1) pks(j) abs(pks(j)-master_pks(i)) 1];
                locs(j,2) = 1;
                match = 1;
            else
                j = j+1;
            end
            
            if j > size(pks,1)
                match = 1;
                pkData(i, 1:6) = [master_locs(i) master_pks(i) master_locs(i) RIC(master_locs(i)) abs(master_pks(i)-RIC(master_locs(i))) 2];
                %pkData(i, 1:5) = [master_locs(i) master_pks(i) master_locs(i) 0 pks(i)];
            end
        end
    else
        pkData(i, 1:6) = [master_locs(i) master_pks(i) 0 0 abs(0-master_pks(i)) 1];
    end
end

%insert unmatched pks
unmatchedPks = pks(locs(:,2)==0);
unmatchedLocs = locs(locs(:,2) == 0);
for i=1:size(unmatchedPks,1)
    loc_under = 1;
    j=1;
    while loc_under == 1
        if j > size(pkData,1)
            temp_pkdata = [unmatchedLocs(i) LIC(unmatchedLocs(i)) unmatchedLocs(i) unmatchedPks(i) abs(unmatchedPks(i)-LIC(unmatchedLocs(i))) 3];
            pkData = [pkData; temp_pkdata];
            loc_under = 0;
        end
        if pkData(j, 3) < unmatchedLocs(i)
            j = j + 1;
        else
            temp_pkdata = [unmatchedLocs(i) LIC(unmatchedLocs(i)) unmatchedLocs(i) unmatchedPks(i) abs(unmatchedPks(i)-LIC(unmatchedLocs(i))) 3];
            %temp_pkdata = [unmatchedLocs(i) 0 unmatchedLocs(i) unmatchedPks(i) unmatchedPks(i)];
            pkData = [pkData(1:j-1,:); temp_pkdata; pkData(j:end,:)];
            loc_under = 0;
        end
    end
end

%compute which side is bigger
for i=1:size(pkData,1)
    if pkData(i,2) > pkData(i,4)
        pkData(i,7) = 1;
    else
        pkData(i,7) = 2;
    end
end
end
