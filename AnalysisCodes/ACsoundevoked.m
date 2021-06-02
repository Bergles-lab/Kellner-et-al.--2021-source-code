%% Auditory cortex widefield imaging
% 
% Kellner et al. Neuron, 2021
%
% Calvin Kersbergen, 20210601
% 
% Load widefield movie and pure tone data file, align and average normalized 
% images to pure tone onset, calculate dF/Fo in ROI over AC to tones. 

[fname, dname] = uigetfile('F:\Calvin\Sound evoked\AC widefield astrocytes\*.tif','Multiselect','on');

[~,~,ext] = fileparts([dname, fname]);

%%tif performance 
if strcmp(ext,'.tif') 
    tic;
    img = loadTif([dname fname],16);
    toc;
elseif strcmp(ext,'.czi')
    tic;
    img = bfLoadTif([dname fname]);
    img = imrotate(img,180);
    toc;
else
    disp('Invalid image format.')
end

%%
toneImg = {};
imgorig = img;
img = normalizeImg(img,10,1);

%% Slow trigger on Axiozoom microscope
offset = 8;

%% load params
[fn2 dname2] = uigetfile(dname);
load([dname2 fn2]);
starttone = (params.baselineDur)/100; %5 s silence before first tone
toneISI = params.stimInt/100; % inter-tone interval - 4 s
toneDur = params.stimDur/100; % duration of tone - 1 s
timeBetweenStart = toneISI + toneDur;
before = 10; %frames before tone to analyze 
[freqSort, order] = sort(params.freqs);
toneImg = {};
for i = 1:params.repeats
    for j = 1:params.numFreqs
       startImg = offset + starttone + timeBetweenStart*(j-1) + params.numFreqs * timeBetweenStart * (i-1) - before
       endImg = startImg + timeBetweenStart + before;
       toneImg{i,j} = img(:,:,startImg:endImg);
    end

    %sort order based on random ordering of frequencies presented
    startInd = 1 + (i-1)*params.numFreqs;
    endInd = startInd + params.numFreqs - 1;
    [freqSort, order] = sort(params.freqs(startInd:endInd));
    toneImg(i,:) = toneImg(i,order);
end

avgToneImg = {};
totalImg = [];
[mm,nn,tt] = size(toneImg{1,1});
ACsig = [];
% Get AC ROI based on local vasculature
 [ACmask] = getACmask(img);
 leftInd = find(ACmask);
for i = 1:params.numFreqs
    totalImg = zeros(mm,nn,tt,params.repeats);
    for j = 1:params.repeats
        wimg = toneImg{j,i};
        totalImg(:,:,:,j) = wimg;
        for k = 1:size(wimg,3)
            tempImg = wimg(:,:,k);
            ACsig(j,k,i) = mean(tempImg(leftInd));
        end
    end
    avgToneImg{1,i} = mean(totalImg,4);
end


%makes Brady-bunch style splaying of images versus frequency presented for
% 16 pure tones 

concat = [];
for i = 1:4
   concat = [concat; avgToneImg{1,(4*(i-1)+1)} avgToneImg{1,(4*(i-1)+2)} avgToneImg{1,(4*(i-1)+3)} avgToneImg{1,(4*(i-1)+4)}];
end
implay(concat) % play movie of sound-evoked aligned responses to 16 pure tones

%% plot dFoFo signals from normalized, aligned, averaged movies

figure;
lt_org = [255, 166 , 38]/255;
dk_org = [255, 120, 0]/255;
lt_blue = [50, 175, 242]/255;
dk_blue = [0, 13, 242]/255;
for i = 1:16
    subplot(4,4,i)
    wimg = avgToneImg{1,i};
     for j = 1:size(wimg,3)
         tempImg = wimg(:,:,j);
        % avgAC(i,j) = mean(tempImg(leftInd));
     end
    imagesc(mean(wimg(:,:,15:25),3));
    colormap('default');
    caxis([-0.1 .3]);
   
end

% Get AC ROI based on temporally projected dFoF signal to all frequencies
 [ACmask] = getACmask(img);
 leftInd = find(ACmask);
 for i = 1:16
    wimg = avgToneImg{1,i};
    for j = 1:size(wimg,3)
        tempImg = wimg(:,:,j);
         avgAC(i,j) = mean(tempImg(leftInd));
    end
 end
% save mean dFoF responses for whole AC ROI
defaultDir = 'F:\Calvin\Sound evoked\AC widefield astrocytes';
save([defaultDir '\413_2gcamp_0dB.mat'],'avgAC')

% Plot averaged dFoF ROI responses to pure tones
figure
for i = 1:16
    subplot(1,16,i);
    plot(avgAC(i,:),'Color',lt_blue,'LineWidth',2);
    patch([10 10 20 20], [1 -6 -6 1],'k','EdgeColor','none','FaceAlpha',0.2);
    ylim([0 0.4])
end

copyACsig = ACsig;
for j = 1:size(ACsig,1)
    copyACsig(j,:,:) = copyACsig(j,:,:) - (j) * 0.4;
end

% plot individual traces (not averages) relative to tone timing
lt_org = [255, 166 , 38]/255;
dk_org = [255, 120, 0]/255;
lt_blue = [50, 175, 242]/255;
dk_blue = [0, 13, 242]/255;
sorted = sort(params.freqs);
fig = figure;
for i = 1:16
    subplot(1,16,i);
    plot(copyACsig(:,:,i)','Color','k');
    hold on;
    %plot(copySigR(:,:,i)','Color',lt_blue); 
    plot(avgAC(i,:),'Color',lt_blue,'LineWidth',2);
    ylim([-10*.4-.1 0.5]);
    xlim([0 60]);
    patch([10 10 20 20], [1 -6 -6 1],'k','EdgeColor','none','FaceAlpha',0.2);
    yticklabels('');
    title([sprintf('%0.3f',sorted(i*params.repeats)/1000) ' kHz']);
    axis off;
end

fig.Units = 'inches';
fig.Position = [2 2 12 8];


