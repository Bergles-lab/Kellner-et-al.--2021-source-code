%Kellner et al., 2021
%dual color imaging of astrocytes and neurons 
%analysis using grid ROIs
close all

%% load files 
[fn,dname1] = uigetfile('*.czi;*.lsm;*.tif','Open gcamp file');
[pathstr, name, ext] = fileparts([dname1 fn]);
bf = bfopen([dname1 fn]);
[m,n] = size(bf{1}{1});
t = size(bf{1},1);
img = zeros(m,n,t,'int16');
for i=1:t
    img(:,:,i) = bf{1}{i};
end
clear bf;

[fn,dname1] = uigetfile('*.czi;*.lsm;*.tif','Open RGECO file');
[pathstr, name, ext] = fileparts([dname1 fn]);
bf = bfopen([dname1 fn]);
[m,n] = size(bf{1}{1});
t = size(bf{1},1);
imgRG = zeros(m,n,t,'int16');
for i=1:t
    imgRG(:,:,i) = bf{1}{i};
end
clear bf;

imgCrop=img(13:end,13:end,:);
imgRGCrop=imgRG(13:end,13:end,:);

imgdA=imgCrop(:,:,1:600);
imgdN=imgRGCrop(:,:,1:600);

%% Grid ROIS
figure; imshow(mean(imgdA,3)/mean(max(mean(imgdA,3))));
widthImg = size(imgdA,2);
heightImg = size(imgdA,1);
%0.8303 microns per pixel - 425x425 microns sq=25
sizeSq = 25; %the diagonal line in the grid comes out to 58um
[positiveIndices] = getGrid(widthImg,heightImg,sizeSq);
T=size(imgdA,3);
for i=1:size(positiveIndices,1)
    hold on;
    tempX=positiveIndices(i,1:2:end);
    tempY=positiveIndices(i,2:2:end);
    plot(tempX,tempY,'Color','g');
    if positiveIndices(i,6)<widthImg
        xs = [positiveIndices(i,2):positiveIndices(i,6)];
    else
        xs=[positiveIndices(i,2):widthImg];
    end
    if positiveIndices(i,3)<heightImg
        ys = [positiveIndices(i,1):positiveIndices(i,3)];
    else
        ys=[positiveIndices(i,1):heightImg];
    end
    
    temp = squeeze(mean(mean(imgdA(xs,ys,:),2),1));
    roisAst(:,i) = (msbackadj([1:T]',temp,'WindowSize',45,'StepSize',45,'Showplot',0));

end
title('Astrocytes')

figure; imshow(mean(imgdN,3)/mean(max(mean(imgdN,3))));
for i=1:size(positiveIndices,1)
    hold on;
    tempX=positiveIndices(i,1:2:end);
    tempY=positiveIndices(i,2:2:end);
    plot(tempX,tempY,'Color','g');
    
    if positiveIndices(i,6)<widthImg
        xs = [positiveIndices(i,2):positiveIndices(i,6)];
    else
        xs=[positiveIndices(i,2):widthImg];
    end
    if positiveIndices(i,3)<heightImg
        ys = [positiveIndices(i,1):positiveIndices(i,3)];
    else
        ys=[positiveIndices(i,1):heightImg];
    end
    
    temp = squeeze(mean(mean(imgdN(xs,ys,:),2),1));
    roisNeur(:,i) = (msbackadj([1:T]',temp,'WindowSize',45,'StepSize',45,'Showplot',0));
end
title('Neurons')

%% plot
[t,n] = size(roisAst);
numToShow = n;
figure('Position',[200,150,1200,800]); plot(roisAst(:,1:numToShow) - .5*repmat(1:numToShow,t,1),'Color','k'); 
ylim([-128 0]); 
title('Astrocytes')

[t,n] = size(roisNeur);
numToShow = n;
figure('Position',[200,150,1200,800]); plot(roisNeur(:,1:numToShow) - .5*repmat(1:numToShow,t,1),'Color','k');ylim([-128 0]); title('Neurons')

%% find peaks
[t,n] = size(roisAst);
dataAstMat=[]; %pks, frms,hw,roi
dataNeurMat=[]; %pks, frms,hw,roi
for cc=1:size(roisNeur,2) %go through each ROI
    threshAst=median(roisAst(:,cc))+3*mad(roisAst(:,cc));    
    [pksAst,locsAst,wAst] = findpeaks((double(roisAst(:,cc))),1:t,'MinPeakHeight',...
        threshAst,'MinPeakProminence', threshAst/4,'WidthReference','halfprom');  
    if ~isempty(pksAst)
    figure; findpeaks((double(roisAst(:,cc))),1:t,...
        'MinPeakHeight',threshAst,'MinPeakProminence', threshAst/4,...
        'WidthReference','halfprom','Annotate','extents');%halfheight
    ax = gca;
    widthHandle = findobj(ax, 'Tag', 'HalfProminenceWidth');%HalfHeightWidth
    widthLocA = [widthHandle.XData(1:3:end)', widthHandle.XData(2:3:end)'];
    else
        widthLocA=[];
    end
    dataAstMat=[dataAstMat;[pksAst,locsAst',wAst',repmat(cc,length(locsAst),1),widthLocA]];
    close
    threshNeur=median(roisNeur(:,cc))+ 1*mad(roisNeur(:,cc));
    [pksNeur,locsNeur,wNeur] = findpeaks((double(roisNeur(:,cc))),1:t,'MinPeakHeight',...
        threshNeur,'WidthReference','halfprom');
        dataNeurMat=[dataNeurMat;[pksNeur,locsNeur',wNeur',repmat(cc,length(locsNeur),1)]];
end

%% find common events and which grids are associated with them - Astrocytes
dataAMatSrt=sortrows(dataAstMat,2); %sort according to peak location (frames)
dataAMatSrtNew=dataAMatSrt;
[C,ia,ic] = unique(dataAMatSrtNew(:,2)); %find the unique values
a_countsA = accumarray(ic,1); %count how many times the unique values appear
value_countsA = [C(1:end-1),a_countsA(1:end-1)];
eLocsA=value_countsA(find(value_countsA(:,2)>=10),1);
eventAinfo=struct;
countA=0;
frmRate=2; %hz

%% user defines events first
uiDefEventsA=[];
for aa=1:10:length(eLocsA)
    h=figure('Position',[10 -10 1500 1000]);
    subplot(3,4,1)
    tempInd=find(ismember(dataAMatSrtNew(:,2),eLocsA(aa)));
    tempW=median(dataAMatSrtNew(tempInd,5:6));
    imagesc(squeeze(mean(imgCrop(:,:,floor(tempW(1)):ceil(tempW(2))),3)));        
    colormap gray
    for bb=1:length(tempInd)
        hold on;
        tempX=positiveIndices(dataAMatSrtNew(tempInd(bb),4),1:2:end);
        tempY=positiveIndices(dataAMatSrtNew(tempInd(bb),4),2:2:end);
        pgon = polyshape(tempX,tempY);
        plot(tempX,tempY,'Color',[1 1 1]);
    end
    axis off 
    title(num2str(eLocsA(aa)))
    for a=1:9
        if aa+a<length(eLocsA)
        subplot(3,4,a+1)        
        tempInd=find(ismember(dataAMatSrtNew(:,2),eLocsA(aa+a)));
        imagesc((imgCrop(:,:,eLocsA(aa+a)))); %[5 15]
        colormap gray
        for bb=1:length(tempInd)
            hold on;
            tempX=positiveIndices(dataAMatSrtNew(tempInd(bb),4),1:2:end);
            tempY=positiveIndices(dataAMatSrtNew(tempInd(bb),4),2:2:end);
            pgon = polyshape(tempX,tempY);
            plot(tempX,tempY,'Color',[1 1 1]);            
        end
        axis off 
        title(num2str(eLocsA(aa+a)))
        end
    end
    uiAns=inputdlg('Enter unique event frames');
    uiDefEventsA=[uiDefEventsA,str2num(uiAns{1})];
    close
end
uiDefEventsA=unique(uiDefEventsA); 

%% Now gather all the data based on those events
eventATraceTm=nan(length(uiDefEventsA),21);
 for aa=1:length(uiDefEventsA)   
    tempInd=find(ismember(dataAMatSrtNew(:,2),uiDefEventsA(aa)));    
    eventAinfo(aa).pkFrm=uiDefEventsA(aa);
    eventAinfo(aa).pkAmp=median(dataAMatSrtNew(tempInd,1));
    eventAinfo(aa).pkHW=median(dataAMatSrtNew(tempInd,3))/frmRate; %hw in seconds
    roiInd=dataAMatSrtNew(tempInd,4); %roi index
    eventAinfo(aa).Grds=roiInd;
    tempEventTrace=nan(length(roiInd),21);
    for bb=1:length(roiInd)
        if uiDefEventsA(aa)-10<=0
            tempEventTrace(bb,1:uiDefEventsA(aa)+10)=roisAst(1:uiDefEventsA(aa)+10,roiInd(bb));
        elseif uiDefEventsA(aa)+10>size(imgCrop,3)
            tempEventTrace(bb,1:length(uiDefEventsA(aa)-10:size(imgCrop,3)))=roisAst(uiDefEventsA(aa)-10:end,roiInd(bb));
        else
            tempEventTrace(bb,1:end)=roisAst(uiDefEventsA(aa)-10:uiDefEventsA(aa)+10,roiInd(bb));
        end
    end
    eventATrace(aa,:)=nanmean(tempEventTrace,1);
    
    if uiDefEventsA(aa)-10<=0
        eventATraceTm(aa,1:uiDefEventsA(aa)+10)=1:uiDefEventsA(aa)+10;
    elseif uiDefEventsA(aa)+10>size(imgCrop,3)
        eventATraceTm(aa,1:length(uiDefEventsA(aa)-10:size(imgCrop,3)))=uiDefEventsA(aa)-10:size(imgCrop,3);
    else
        eventATraceTm(aa,1:length(uiDefEventsA(aa)-10:uiDefEventsA(aa)+10))=uiDefEventsA(aa)-10:uiDefEventsA(aa)+10;
    end
    tempInd=[];
    roiInd=[];
 end

%% find common events and which grids are associated with them - neurons
dataNMatSrt=sortrows(dataNeurMat,2); %sort according to peak location (frames)
dataNMatSrtNew=dataNMatSrt;
[C,ia,ic] = unique(dataNMatSrtNew(:,2)); %find the unique values
a_counts = accumarray(ic,1); %count how many times the unique values appear
value_countsN = [C(1:end-1),a_counts(1:end-1), [abs(diff(a_counts))]];
eLocsN=value_countsN(find(value_countsN(:,2)>=10),1);
eventNinfo=struct;
eventNTraceTm=nan(length(eLocsN),21);
countN=0;
frmRate=2; %hz

%% user defines events first
uiDefEvents=[];
for aa=1:10:length(eLocsN)
    h=figure('Position',[10 -10 1500 1000]);
    subplot(3,4,1)
    tempInd=find(ismember(dataNMatSrtNew(:,2),eLocsN(aa)));
    imagesc((imgRGCrop(:,:,eLocsN(aa))),[5 15]);
    colormap gray
    for bb=1:length(tempInd)
        hold on;
        tempX=positiveIndices(dataNMatSrtNew(tempInd(bb),4),1:2:end);
        tempY=positiveIndices(dataNMatSrtNew(tempInd(bb),4),2:2:end);
        pgon = polyshape(tempX,tempY);
        plot(tempX,tempY,'Color',[1 1 1]);
    end
    axis off 
    title(num2str(eLocsN(aa)))
    for a=1:9
        if aa+a<length(eLocsN)
        subplot(3,4,a+1)
        tempInd=find(ismember(dataNMatSrtNew(:,2),eLocsN(aa+a)));
        imagesc((imgRGCrop(:,:,eLocsN(aa+a))),[5,15]);
        colormap gray
        for bb=1:length(tempInd)
            hold on;
            tempX=positiveIndices(dataNMatSrtNew(tempInd(bb),4),1:2:end);
            tempY=positiveIndices(dataNMatSrtNew(tempInd(bb),4),2:2:end);
            pgon = polyshape(tempX,tempY);
            plot(tempX,tempY,'Color',[1 1 1]);            
        end
        axis off 
        title(num2str(eLocsN(aa+a)))
        end
    end
    uiAns=inputdlg('Enter unique event frames');
    uiDefEvents=[uiDefEvents,str2num(uiAns{1})];
    close
end
uiDefEvents=unique(uiDefEvents);  

%% Now gather all the data based on those events
eventNTraceTm=nan(length(uiDefEvents),21);
 for aa=1:length(uiDefEvents)   
    tempInd=find(ismember(dataNMatSrtNew(:,2),uiDefEvents(aa)));    
    eventNinfo(aa).pkFrm=uiDefEvents(aa);
    eventNinfo(aa).pkAmp=median(dataNMatSrtNew(tempInd,1));
    eventNinfo(aa).pkHW=median(dataNMatSrtNew(tempInd,3))/frmRate; %hw in seconds
    roiInd=dataNMatSrtNew(tempInd,4); %roi index
    eventNinfo(aa).Grds=roiInd;
    tempEventTrace=nan(length(roiInd),21);
    for bb=1:length(roiInd)
        if uiDefEvents(aa)-10<=0
            tempEventTrace(bb,1:uiDefEvents(aa)+10)=roisNeur(1:uiDefEvents(aa)+10,roiInd(bb));
        elseif uiDefEvents(aa)+10>size(imgdN,3)
            tempEventTrace(bb,1:length(uiDefEvents(aa)-10:size(imgdN,3)))=roisNeur(uiDefEvents(aa)-10:end,roiInd(bb));
        else
            tempEventTrace(bb,1:end)=roisNeur(uiDefEvents(aa)-10:uiDefEvents(aa)+10,roiInd(bb));
        end
    end
    eventNTrace(aa,:)=nanmean(tempEventTrace,1);
    
    if uiDefEvents(aa)-10<=0
        eventNTraceTm(aa,1:uiDefEvents(aa)+10)=1:uiDefEvents(aa)+10;
    elseif uiDefEvents(aa)+10>size(imgdN,3)
        eventNTraceTm(aa,1:length(uiDefEvents(aa)-10:size(imgdN,3)))=uiDefEvents(aa)-10:size(imgdN,3);
    else
        eventNTraceTm(aa,1:length(uiDefEvents(aa)-10:uiDefEvents(aa)+10))=uiDefEvents(aa)-10:uiDefEvents(aa)+10;
    end
 end

%% put into struct
dualstruct(f).Name=name;
dualstruct(f).Raw.AstStats=eventAinfo;
dualstruct(f).Raw.eventATrace=eventATrace;
dualstruct(f).Raw.eventATraceTm=eventATraceTm;
dualstruct(f).Raw.newAstInds=uiDefEventsA;
dualstruct(f).Raw.AllFrmAstAmp=aAmp; %3 frames after neuron peak
dualstruct(f).Raw.NeurStats=eventNinfo;
dualstruct(f).Raw.eventNTrace=eventNTrace;
dualstruct(f).Raw.eventNTraceTm=eventNTraceTm;
dualstruct(f).Raw.newNeurInds=uiDefEvents;
dualstruct(f).Raw.AllFrmNAmp=nPk;

%% save
save(['dualStruct.mat'],'dualstruct','-v7.3');