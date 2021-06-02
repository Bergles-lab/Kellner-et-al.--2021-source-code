%Kellner et al., 2021

%Analyze auditory cortex WF data
close all
clearvars -except ACICstruct
f=1;

%% load file

[fn,dname1] = uigetfile('*.tif','Open tiff movie');
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
        %         img = imrotate(img,180);
        toc;
        %% take only 10 min (choose manually)
        img2=img(:,:,6000:12000);
        %% bleach correction
        sampRate = 10; %sampling rate in Hz
        [imgBC] = bleachCorrect(img2, sampRate);
        imgBC=imrotate(imgBC,180);
    case '.tif'
        img = loadTif([dname1 fn],16);
        %% bleach correction
        sampRate = 10; %sampling rate in Hz
        [imgBC] = bleachCorrect(img, sampRate);
        figure;
        imagesc(mean(imgBC,3)); axis image; colormap gray
        clear img
end

%% crop image for movement analysis
figure;
imagesc(mean(imgBC,3)); axis image; colormap gray
brect = imrect(gca,[0,0,150,150]);%150
setResizable(brect,1);
wait(brect);
pos = getPosition(brect);
pos = int16(round(pos));
crimg=imgBC((pos(2)):(pos(2)+pos(4)),pos(1):(pos(1)+pos(3)),:);
close
figure; imagesc(mean(crimg,3)); axis image; colormap gray
imgCrop=crimg; clear crimg

%% find movement times (from Travis / Vered)
imgCrop=double(imgCrop);
norms = [];
regto = fft2(mean((imgCrop),3));
t = size(imgCrop,3);
parfor j = 1:t
    [output, Greg] = dftregistration(regto,fft2(imgCrop(:,:,j)),1) %the last number will determine how sensitive this is
    norms(j) = norm([1 0 output(3); 0 1 output(4); 0 0 1],2);
end

hit = zeros(size(norms));
for j = 1:t
    if norms(j) > 1
        %look 50 ahead
        if j>1 && j < t-30
            if sum(norms(j+1:j+30)-1) > 0
                hit(j:j+30) = 1;
            elseif hit(j-1) == 1
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

%% down sample image
imgBCdn=imresize(imgBC,0.5); %down sample by 50% to make more manageable

%% get dFF
[dFoF,Fo]=normalizeImg(double(imgBCdn));

%% define ROI around IC
[m,n,t]=size(imgBCdn);
figure; h_im=imagesc(squeeze(mean(imgBCdn,3)));
h = imellipse();
wait(h);
indices = find(h.createMask);
ICsignal2 = zeros(1,t);
for i=1:t
    Xw = dFoF(:,:,i);
    ICsignal2(1,i) = squeeze(mean(Xw(indices)));
end

%% filter the IC trace
sampRate=10;
win=20;
T=[0:1/sampRate:(t/sampRate)-1/sampRate];
ICfilt = msbackadj(T',smooth(ICsignal2),'SHOWPLOT',0,'WindowSize',win,'StepSize',win);

%% remove movement artifacts
ICsignalNMVM=ICfilt;
ICsignalNMVM(mvmInd)=0;
ICsignalNMVM2=ICfilt;
ICsignalNMVM2(mvmInd)=nan;

%% run the correlation
reX = single(reshape(imgBCdn,m*n,t));
corrmat = corr(ICsignalNMVM, reX');
corrmat = reshape(corrmat,m,n);
figure; imagesc(corrmat); colorbar
truesize; axis off

%% get ROI around AC based on correlation
figure(5);
AC = imellipse(gca,[40,40,60,60]);
setResizable(AC,1);
wait(AC);
indices2 = find(AC.createMask);
ACsignal = zeros(1,t);
for i=1:t
    Xw = dFoF(:,:,i);
    ACsignal(1,i) = squeeze(mean(Xw(indices2)));
end

ACfilt = msbackadj(T',smooth(ACsignal),'SHOWPLOT',0,'WindowSize',win,'StepSize',win);
ACsignalNMVM=ACfilt;
ACsignalNMVM(mvmInd)=0;
ACsignalNMVM2=ACfilt;
ACsignalNMVM2(mvmInd)=nan;


%% overlay the traces
figure('Position',[500,900,1500,400]); plot(ICsignalNMVM); hold on; plot(ACsignalNMVM,'r'); 
legend('IC','AC')
xlim([0 6001])

%% calculate correlations (this is on the whole signal)
ICACcorr=corr(ICsignalNMVM2,ACsignalNMVM2,'rows','complete');

%% cross correlations (this is on the whole signal)
[acICcor,acIClag] = xcorr(ICsignalNMVM-mean(ICsignalNMVM),ACsignalNMVM-mean(ACsignalNMVM),'normalized'); %need to use zeros here for the movement related events - doesn't work with nans

[~,I] = max(abs(acICcor));
lagDiff = acIClag(I)
timeDiff = lagDiff/sampRate

figure
plot(acIClag,acICcor)

%% colorful correlation map for IC and AC with seeds - possibly also traces
%% crop AC image 
figure(5);
brect = imrect(gca,[0,0,60,60]);
setResizable(brect,1);
wait(brect);
pos = getPosition(brect);
pos = int16(round(pos));
crimg=dFoF((pos(2)):(pos(2)+pos(4)),pos(1):(pos(1)+pos(3)),:);
figure; imagesc(mean(crimg,3)); axis image; colormap gray
imgACCrop=crimg; clear crimg

%% dfof on cropped IC image
[imgICCrop,Fo]=normalizeImg((imgCrop));

%% filter
% filter AC movie
s = size(imgACCrop);
ACmovieFilt = zeros(s);
for j = 1:s(1)
    parfor k = 1:s(2)
        ACmovieFilt(j,k,:) = msbackadj(T',squeeze(imgACCrop(j,k,:)),'SHOWPLOT',0,'WindowSize',win,'StepSize',win);
    end
end

ACmovieFiltMS = zeros(s);
parfor j = 1:s(3)
    t2 = imgaussfilt(ACmovieFilt(:,:,j), 1);
    ACmovieFiltMS(:,:,j) = t2 - mean(mean(t2));
end
ACmovieFiltMS = max(ACmovieFiltMS,0);

% filter IC movie - blur + mean subtraction: pick out spatially
% isolated events
s2 = size(imgICCrop);
ICmovieFilt = zeros(s2);
parfor j = 1:s2(3)
    t2 = imgaussfilt(imgICCrop(:,:,j), 1);
    ICmovieFilt(:,:,j) = t2 - mean(mean(t2));
end

%% run PCA on IC, AC movies
figure;
ICmoviePCA = runPCA(ICmovieFilt.*2000);
close
figure;
ACmoviePCA = runPCA(ACmovieFiltMS.*2000);
close
%% filter PCA movies
s = size(ACmoviePCA);
ACmovief = zeros(s);
for j = 1:s(1)
    parfor k = 1:s(2)
        ACmovief(j,k,:) = msbackadj(T',squeeze(ACmoviePCA(j,k,:)),'SHOWPLOT',0,'WindowSize',win,'StepSize',win);
    end
end

s2 = size(ICmoviePCA);
ICmovief = zeros(s2);
for j = 1:s2(1)
    parfor k = 1:s2(2)
        ICmovief(j,k,:) = msbackadj(T',squeeze(ICmoviePCA(j,k,:)),'SHOWPLOT',0,'WindowSize',win,'StepSize',win);
    end
end

%% calculate series of IC traces and AC correlation images
%place seeds on IC image
figure; imagesc(squeeze(sum(ICmovief,3)));
[Xic,Yic]=ginput(4);
ICtraces = zeros(4,s(3));
ACcorrImgs = zeros(s(1),s(2),4);
ICcorrImgs = zeros(s2(1),s2(2),4);
ACmoviep = permute(ACmovief, [3 1 2]);
ICmoviep = permute(ICmovief, [3 1 2]);
for j = 1:4
    t2 = squeeze(ICmovief(round(Xic(j)), round(Yic(j)),:));
    t2 = msbackadj(T',t2,'SHOWPLOT',0,'WindowSize',win,'StepSize',win);
    ICtraces(j,:) = t2;
    t2(mvmInd)=nan;
    t3 = [];
    parfor k = 1:s(2)
        t3 = [t3 corr(t2, ACmoviep(:,:,k),'rows','complete')'];
    end
    ACcorrImgs(:,:,j) = t3;
    t3 = []; 
    parfor k = 1:s2(2)
        t3 = [t3 corr(t2, ICmoviep(:,:,k),'rows','complete')'];
    end
    ICcorrImgs(:,:,j) = t3;
end

tempAC1=ACcorrImgs(:,:,1);
tempAC2=ACcorrImgs(:,:,2);
tempAC3=ACcorrImgs(:,:,3);
tempAC=cat(3,tempAC1,tempAC2,tempAC3);
figure; imagesc(rescale(tempAC,-0.2,0.9))
axis off; axis image;

tempIC1=ICcorrImgs(:,:,1);
tempIC2=ICcorrImgs(:,:,2);
tempIC3=ICcorrImgs(:,:,3);
tempIC=cat(3,tempIC1,tempIC2,tempIC3);
figure; imagesc(rescale(tempIC,-0.2,0.9))
axis off; axis image; 
hold on; plot(Xic(1:3),Yic(1:3),'ok','MarkerFaceColor','k')

figure; plot(ICtraces(1,:)'); hold on; plot(ICtraces(2,:)-20');plot(ICtraces(3,:)-40');

%% save to structure
ACICstruct(f).Name=fn;
ACICstruct(f).Corr.Img=corrmat;
ACICstruct(f).Corr.ICAC=ICACcorr;
ACICstruct(f).Corr.ICVC=ICVCcorr;
ACICstruct(f).Corr.ICSC=ICSCcorr;
ACICstruct(f).Corr.Xcorr=acICcor;
ACICstruct(f).Corr.XcorrLag=acIClag;
ACICstruct(f).Corr.ICACPk=ICACcorrPks;
ACICstruct(f).Corr.ICVCPk=ICVCcorrPks;
ACICstruct(f).Corr.ICSCPk=ICSCcorrPks;
ACICstruct(f).Corr.XcorrPk=acICcorPk;
ACICstruct(f).Corr.XcorrLagPk=acIClagPk;
ACICstruct(f).Event.IC=ICeventMat;
ACICstruct(f).Event.AC=ACeventMat;
ACICstruct(f).Seeds.seed=[Xic,Yic];
ACICstruct(f).Seeds.ICtraces=ICtraces;
ACICstruct(f).Seeds.ACcorr=ACcorrImgs;
ACICstruct(f).Seeds.ICcorr=ICcorrImgs;

%% save
save([path 'ACICsummary.mat'],'ACICstruct','-v7.3');
