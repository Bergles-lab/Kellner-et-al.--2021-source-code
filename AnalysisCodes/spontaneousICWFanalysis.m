%Kellner et al., 2021
%% analysis file of spontaneous IC activity in neurons and astrocytes
% adapted from Travis Babola's codes by Vered May 2017
close all

%% load file
[fn,dname1] = uigetfile('*.tif','Open BC movie');
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
        img = imrotate(img,180);
        toc;
%         img2=imresize(img,0.5);
%         img=img2; clear img2
%         [m,n,t] = size(img);
        %% crop image
        h = figure;
        imagesc(mean(img,3)); axis image; colormap gray
        brect = imrect(gca,[0,0,size(img,2)-12,150]);%150
        setResizable(brect,1);
        wait(brect);
        pos = getPosition(brect);
        pos = int16(round(pos));
        crimg=img((pos(2)):(pos(2)+pos(4)),:,:);
        close
        figure; imagesc(mean(crimg,3)); axis image; colormap gray
        img2=crimg; clear crimg       
        
        
        %% bleach correction
        sampRate = 10; %sampling rate in Hz
        [imgBC] = bleachCorrect(img2, sampRate);
         clear img2
    case '.tif'
        img = loadTif([dname1 fn],16);
        imgBC=img;
        figure;
        imagesc(mean(img,3)); axis image; colormap gray
        clear img
end

%% crop section for background activity in SC
h2 = figure;
imagesc(mean(imgBC,3)); axis image; colormap gray
brect = imrect(gca,[0,0,100,89]);%40 %this corresponds to the number of pixels in LICmask and RICmask
setResizable(brect,1);
wait(brect);
pos = getPosition(brect);
pos = int16(round(pos));
crimg2=imgBC((pos(2)):(pos(2)+pos(4)),pos(1):(pos(1)+pos(3)),:);
close
figure; imagesc(mean(crimg2,3)); axis image; colormap gray
imgArt=crimg2; clear crimg2

%% Find times of movement automatically using dft registration (Guizar-Sicairos et al., 2008)
imgBC=double(imgBC);
temp = {};
tforms = [];
norms = [];
regd =[];
regto = fft2(mean((imgBC),3));
t = size(imgBC,3);
parfor i = 1:t
    [output, Greg] = dftregistration(regto,fft2(imgBC(:,:,i)),1) %the last number will determine how sensitive this is
    norms(i) = norm([1 0 output(3); 0 1 output(4); 0 0 1],2);
end

hit = zeros(size(norms));
for i = 1:t
    if norms(i) > 1
        %look 50 ahead       
        if i>5 && i < t-30
            if sum(norms(i+1:i+30)-1) > 0
                hit(i-5:i+30) = 1;
            elseif hit(i-5) == 1 
%                 hit(i:i+5) = 1;
%                 hit(i+6:i+56)=0;
                hit(i:i+1) = 1;
                hit(i+2:i+30)=0;
            else
                hit(i) = 0;
            end
        end
    end
end

figure; plot(norms-1); hold on; plot(hit)
mvmInd=find(hit>0);

%% dFF
 [dFoF2, Fo] = normalizeImg2(imgBC); %using median as Fo
[dFoFArt,FoArt]=normalizeImg2(double(imgArt));
 
%% get L/R/ctx masks and signals
[LICmask, RICmask, ctxmask] = getROImasks(imgBC);

%% Whole IC analysis
LICmask=double(LICmask);
LICmask(LICmask==0)=nan;
LICsignal = double(dFoF2).*LICmask;
LICsignal2 = squeeze(median(median(LICsignal,2,'omitnan'),'omitnan'));

RICmask=double(RICmask);
RICmask(RICmask==0)=nan;
RICsignal = double(dFoF2).*RICmask;
RICsignal2 = squeeze(median(median(RICsignal,2,'omitnan'),'omitnan'));

artsignal = squeeze(median(median((dFoFArt),2,'omitnan'),'omitnan'));

%% Remove baseline
sampRate=10;
win=5;
t=[0:1/sampRate:(length(RICsignal2)/sampRate)-1/sampRate];
RICfilt=msbackadj(t',RICsignal2,'SHOWPLOT',0,'WindowSize',win,'StepSize',win);
msbackadj(t',RICsignal2,'SHOWPLOT',1,'WindowSize',win,'StepSize',win);
figure; plot(RICsignal2); hold on; plot(RICfilt,'r'); title(['RIC Win=',num2str(win)])
LICfilt=msbackadj(t',LICsignal2,'SHOWPLOT',0,'WindowSize',win,'StepSize',win);
figure; plot(LICsignal2); hold on; plot(LICfilt,'r'); title(['LIC Win=',num2str(win)])
artfilt = msbackadj(t',(artsignal),'SHOWPLOT',0,'WindowSize',win,'StepSize',win);
figure; plot(artsignal); hold on; plot(artfilt,'r'); title(['Artifact Win=',num2str(win)])

%% remove movement indices
LICfilt(mvmInd)=0;
RICfilt(mvmInd)=0;
artfilt(mvmInd)=0;

%% Find IC peaks and plot

cellAns=questdlg('Cell type?','','Neuron','Astrocyte','Astrocyte');
[wholeROIinfo,pkData] = findICpeaksdFoFVK_new([LICfilt RICfilt artfilt],1,cellAns);
%pkData: col 1: LIC loc, 2: LIC pk, 3:RIC loc, 4:RIC pk, 5:delta, 6: peak type (1=matched, 2=LIC only, 3=RIC only) 7: which is bigger (1=LIC, 2=RIC)
wholeROIinfo.frmNum=length(LICfilt)-length(mvmInd);

%% look at bilaterality of peaks
%Dominance ratio: smaller divided by larger
if ~isempty(pkData)
matchPkInd=find(pkData(:,6)==1);
domRatio=[];
for cc=1:length(matchPkInd)
    if pkData(matchPkInd(cc),2)>pkData(matchPkInd(cc),4) %left larger than right
        domRatio(cc)=pkData(matchPkInd(cc),4)/pkData(matchPkInd(cc),2); 
    else
        domRatio(cc)=pkData(matchPkInd(cc),2)/pkData(matchPkInd(cc),4); 
    end
end
leftDomInd=find(pkData(:,6)==2);
rightDomInd=find(pkData(:,6)==3);
end
%correlation
[rho,pcorr]=corr(LICfilt,RICfilt,'tail','right');

% put into bilaterality structure:
bilatInfo=struct;
if ~isempty(pkData)
    bilatInfo.dominance=domRatio;
    bilatInfo.matched=length(domRatio)/size(pkData,1);
    bilatInfo.LeftDom=length(leftDomInd)/size(pkData,1);
    bilatInfo.RightDom=length(rightDomInd)/size(pkData,1);
else
    bilatInfo.dominance=[];
    bilatInfo.matched=[];
    bilatInfo.LeftDom=[];
    bilatInfo.RightDom=[];
end
bilatInfo.Corr=[rho,pcorr];
