%Kellner et al., 2021

%time line for AC astrocyte activity P4,P7,P11,P14-15
close all
clearvars -except ACstruct
f=1;
%% load file
[fn,dname1] = uigetfile('*.czi','Open raw movie');
[pathstr, name, ext] = fileparts([dname1 fn]);
switch ext
    case '.czi'
        bf = bfopen([dname1 fn]);
        tic;
        [m,n] = size(bf{1}{1});
        t = size(bf{1},1);
        img = zeros(m,n,t,'int16');
        for i=1:t
            img(:,:,i) = bf{1}{i};
        end
        clear bf;
        toc;
        %% take only 10 min (choose manually)
        if size(img,3)>6000
        img2=img(:,:,1:6000);
        else
            img2=img;
        end
        %% bleach correction
        sampRate = 10; %sampling rate in Hz
        [imgBC] = bleachCorrect(img2, sampRate);
    case '.tif'
        imgBC = loadTif([dname1 fn],16);
        figure;
        imagesc(mean(imgBC,3)); axis image; colormap gray
        clear img
end
%% crop image for grid analysis
figure;
imagesc(mean(imgBC,3)); axis image; colormap gray
brect = imrect(gca,[0,0,249,249]);%250
setResizable(brect,1);
wait(brect);
pos = getPosition(brect);
pos = int16(round(pos));
crimg=imgBC((pos(2)):(pos(2)+pos(4)),pos(1):(pos(1)+pos(3)),:);
close
imgCrop=crimg; clear crimg

%% get dFF
[dFoF,Fo]=normalizeImg(double(imgCrop));

%% Try grid ROIS
figure; imshow(mean(imgCrop,3)/mean(max(mean(imgCrop,3))));
widthImg = size(imgCrop,2);
heightImg = size(imgCrop,1);
%7.69 microns per pixel for P11 
sizeSq = 50; %divid into 25
% sizeSq = 30; %divid into 25

[positiveIndices] = getGrid(widthImg,heightImg,sizeSq);
T=size(imgCrop,3);
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
    
    temp = squeeze(mean(mean(dFoF(xs,ys,:),2),1));
    rois(:,i) = smooth(msbackadj([1:T]',temp,'WindowSize',20,'StepSize',20,'Showplot',0));
end

%% find movement times (from Travis / Vered)
imgCrop=double(imgCrop);
norms = [];
regto = fft2(mean((imgCrop),3));
t = size(imgCrop,3);
parfor j = 1:t
    [output, Greg] = dftregistration(regto,fft2(imgCrop(:,:,j)),5) %the last number will determine how sensitive this is
    norms(j) = norm([1 0 output(3); 0 1 output(4); 0 0 1],2);
end
tempNorms=diff(norms);

hit = zeros(size(norms));
for j = 1:t
    if norms(j) > 1
        %look 50 ahead
        if j>5 && j < t-30
            if sum(norms(j+1:j+30)-1) > 0
                hit(j-5:j+30) = 1;
            elseif hit(j-5) == 1
                %                 hit(i:i+5) = 1;
                %                 hit(i+6:i+56)=0;
                hit(j:j+1) = 1;
                hit(j+2:j+30)=0;
            else
                hit(j) = 0;
            end
        end
    end
end
figure; plot(norms-1); hold on; plot(hit)
mvmInd=find(hit>0);
ACstruct(f).Mvm=mvmInd;

%% plot grids
[t,n] = size(rois);
numToShow = 25;
figure('Position',[200,150,1200,800]); plot(rois(:,1:numToShow) - .1*repmat(1:numToShow,t,1),'Color','k'); 

%% remove movement
roisNMVM=rois;
roisNMVM(mvmInd,:)=0;
roisNMVM2=rois;
roisNMVM2(mvmInd,:)=nan;

%% find peaks
[t,n] = size(roisNMVM);
dataMat=[]; %pks, frms,hw,roi

for cc=1:size(roisNMVM,2) %go through each ROI
    threshAst(cc)=nanmedian(roisNMVM2(:,cc))+3*mad(roisNMVM2(:,cc));
   [pksAst,locsAst,wAst,pAst] = findpeaks((double(roisNMVM(:,cc))),1:t,...
       'MinPeakProminence', threshAst(cc),'MinPeakWidth',10,'WidthReference','halfheight');
    dataMat=[dataMat;[pksAst,locsAst',wAst',pAst,repmat(cc,length(locsAst),1)]];

end
dataMatSrt=sortrows(dataMat,2); %sort according to peak location (frames)
tempDiff=diff(dataMatSrt(:,2));
tempDiff=[0;tempDiff]; %to make this the correct n
%% find the common events and average over them
tempSame=[];
count=0;
dataMatNew=zeros(1,size(dataMatSrt,2));
for n=1:length(tempDiff)
    if tempDiff(n)<20
        tempSame=[tempSame;dataMatSrt(n,:)];
    else
        count=count+1;
        if tempDiff(n-1)>20
            dataMatNew(count,:)=(tempSame);
            tempSame=[];
            tempSame=[tempSame;dataMatSrt(n,:)];            
        else
            dataMatNew(count,:)=median(tempSame,1);
            tempSame=[];
            tempSame=[tempSame;dataMatSrt(n,:)];
        end
    end
end
%% save to structure
ACstruct(f).Name=fn;
sensorList={'GCaMP3','GCaMP6s','RCAMP','RGECO','Virus-GCaMP6s','GCaMP6short'};
sensorAns=listdlg('ListString',sensorList,'PromptString','Sensor:');
ACstruct(f).Sensor=sensorList{sensorAns};
promoterList={'GLAST-CreER','SNAP25','Aldh1l1-CreER','Pax2Cre','GFAP','GLAST-KI','hSyn','Thy1'};
promoterAns=listdlg('ListString',promoterList,'PromptString','Promoter:');
ACstruct(f).Promoter=promoterList{promoterAns};
ageList={'4','5','6','7','8','9','10','11','12','13','14','15'};
ageAns=listdlg('ListString',ageList,'PromptString','Age (days):');
ACstruct(f).Age=ageList{ageAns};
sexList={'Male','Female','?'};
sexAns=listdlg('ListString',sexList,'PromptString','Sex:');
ACstruct(f).Sex=sexList{sexAns};
tmxList={'None','Once','Twice','3x Mom','3x Mom + 2xpup','4HT twice'};
tmxAns=listdlg('ListString',tmxList,'PromptString','Tamoxifen injection:');
ACstruct(f).Tamoxifen=tmxList{tmxAns};
cmntList={'Good','OK','Slightly cloudy','Pretty cloudy','Bad','1024x1024','No/low activity','Animal moving a lot','pinworms'};
cmntAns=listdlg('ListString',cmntList,'PromptString','Brain state');
ACstruct(f).comment=cmntList{cmntAns};
%%
ACstruct(f).Trace=ACeventMat;
ACstruct(f).Events=dataMatNew; %amp,frames,FWHM,prominence,grid
% ACstruct(f).Events10frW=dataMatNew; %amp,frames,FWHM,prominence,grid

ACstruct(f).frmNum=length(ACfilt)-length(mvmInd);
ACstruct(f).Thresh=median(threshAst);

%% save
save(['ACtimelineGrid.mat'],'ACstruct','-v7.3');