%Kellner et al., 2021
%slice analysis using Travis Babola's grid code
load('sliceTimes.mat') %this is just data copied from excel about how long each condition was
load('sliceSize.mat') %this is just data copied from excel about the size of the image
%%
close all
clearvars -except sliceStruct manipTMat Size
f=1;

%% load files
[fn,dname1] = uigetfile('*.tif','Open slice file');
[pathstr, name, ext] = fileparts([ dname1 fn]);
bf = bfopen([dname1 fn]);
[m,n] = size(bf{1}{1});
t = size(bf{1},1);
img = zeros(m,n,t,'int16');
for i=1:t
    img(:,:,i) = bf{1}{i};
end
clear bf;
% %crop image (this is because of registration issues)
imgCrop=img(13:end,13:end,:);
imgCropBC = bleachCorrect(imgCrop,1);
manipNum=length(find(~isnan(manipTMat(f,:))));

%% make grid ROIs
figure; imshow(mean(imgCropBC,3)/mean(max(mean(imgCropBC,3))));
widthImg = size(imgCropBC,2);
heightImg = size(imgCropBC,1);
switch Size(f)
    case {118,121}
        sizeSq = 100;
    case 425
        sizeSq = 25;
    case {170,143,212,141,184}
        sizeSq = 50;
end

positiveIndices = getGrid(widthImg,heightImg,sizeSq);
for i=1:size(positiveIndices,1)
    hold on;
    plot(positiveIndices(i,1:2:end),positiveIndices(i,2:2:end),'Color','w');
end
% mic2px=0.2306; %0.2306: 118x118 - sq=100 %0.8303 - 425x425 sq=25 %0.3321 170x170 - sq=50
% 143X143 0.2798 - sq = 50; 121x121 0.2372 sq=100 ; 212 x 212 0.4151 sq=50
%141x141 0.2768 sq=50 ; 184x184 0.361 sq=50
%1 pixel per X microns
%try to get 15-25 micron squares

%% get dFoF and for each ROI
rois = normalizeGridImg(imgCropBC,10,positiveIndices);

%% plot
[t,n] = size(rois);
manipTime(1)={1:manipTMat(f,1)};
countT=manipTMat(f,1);
for a=2:manipNum
    manipTime(a)={countT+1:countT+manipTMat(f,a)};
    countT=countT+manipTMat(f,a);
end
numToShow = round(n/4);
jump=2;
figure;
l=1;
plot(manipTime{l},rois(manipTime{l},1:jump:numToShow) - 0.05*repmat(1:jump:numToShow,length(manipTime{l}),1),'Color','k');
for aa=l+1:manipNum
    hold on;
    plot(manipTime{aa}+5*aa,rois(manipTime{aa},1:jump:numToShow) - 0.05*repmat(1:jump:numToShow,length(manipTime{aa}),1),'Color','k');
end

%%
%% put data into structure
sliceStruct(f).dir=dname1;
sliceStruct(f).fileName=fn;
sensorList={'GCaMP3','GCaMP6s'};
sensorAns=listdlg('ListString',sensorList,'PromptString','Sensor:');
sliceStruct(f).Sensor=sensorList{sensorAns};
promoterList={'GLAST-CreER','Aldh1l1-CreER'};
promoterAns=listdlg('ListString',promoterList,'PromptString','Promoter:');
sliceStruct(f).Promoter=promoterList{promoterAns};
ageList={'6','7','8','9','10','11','12','13','14','15'};
ageAns=listdlg('ListString',ageList,'PromptString','Age (days):');
sliceStruct(f).Age=ageList{ageAns};
genotypeList={'Control','mGluR5','mGluR3'};
genoAns=listdlg('ListString',genotypeList,'PromptString','Genotype:');
sliceStruct(f).Genotype=genotypeList{genoAns};
alleleList={'+/+','-/-','+/-','fl/+','fl/fl'};
alleleAns=listdlg('ListString',alleleList,'PromptString','Allele:');
sliceStruct(f).Allele=alleleList{alleleAns};
animalAns=inputdlg('Animal number');
sliceStruct(f).Animal=animalAns{:};
sliceAns=inputdlg('Slice number');
sliceStruct(f).SliceNum=sliceAns{:};
manips={'aCSF','TTX','mGluR5 Agonist','mGluR5 Agonist+Antagonist','mGluR3 Agonist','mGluR3 Agonist+Antagonist','Thapsigargin + TTX','Thaps + mGluR3 Agonist','NE'};
manipAns=listdlg('ListString',manips,'PromptString','Choose the manipulations');
sliceStruct(f).Manipulation=manips(manipAns);
areaList={'ECIC','CIC'};
areaAns=listdlg('ListString',areaList,'PromptString','Area:');
sliceStruct(f).Area=areaList{areaAns};

%% calculate area under curve
agInd=find(strcmp(sliceStruct(f).Manipulation,manips{3}) | strcmp(sliceStruct(f).Manipulation,manips{5}));
neInd=find(strcmp(sliceStruct(f).Manipulation,manips{9}));

figure; plot(rois(manipTime{agInd},1:25))
onAns=inputdlg('Input On time','',1,{'250:350'});
if isempty(neInd)
    for aa=1:manipNum
        auc{aa}=trapz(rois(manipTime{aa}(1)+str2num(onAns{:}),:));
    end
else
    for aa=1:manipNum-1
        auc{aa}=trapz(rois(manipTime{aa}(1)+str2num(onAns{:}),:));
    end
    auc{manipNum}=trapz(rois(manipTime{manipNum}(1)+(1:100),:));
end
%% find peaks
respROIManip=cell(1,manipNum);
pksROIManip=cell(1,manipNum);
eventManip=cell(1,manipNum);
riseEvent=cell(1,manipNum);
fallEvent=cell(1,manipNum);
durPkManip=cell(1,manipNum);
% tEvent=8;%for 32 frames total
% sampRate=2;
for cc=1:size(rois,2) %go through each ROI
    thresh=median(rois(:,cc))+5*mad(rois(:,cc));
    %     figure; plot(rois(:,cc)); hold on; plot([1:size(rois,1)],repmat(thresh,1,size(rois,1)))
    [pks{cc},locs{cc},w{cc}] = findpeaks(double(rois(:,cc)),1:t,'MinPeakHeight',thresh,...
        'MinPeakProminence',thresh/2,'WidthReference','halfheight','Annotate','Extent');
    if ~isempty(pks{cc})
%         avgPksTemp=nanmean(pks{cc});
        for d=1:length(locs{cc})%go through each peak and check which manipulation it belongs to            
            for dd=1:size(manipTime,2) %go through each manipulation time stamp
                if ismember(locs{cc}(d),manipTime{dd})
                    respROIManip{dd}=[respROIManip{dd}, cc];
                    pksROIManip{dd}=[pksROIManip{dd}; pks{cc}(d)];
                    durPkManip{dd}=[durPkManip{dd}; w{cc}(d)];
                    tEvent=ceil(w{cc}(d))+2; %use the half width for the timing
                    if tEvent>20
                        tEvent=20; %cap at 20
                    end
                    if (locs{cc}(d)-(tEvent))>0 && (locs{cc}(d)+(tEvent))<=t
                        tempEvent=double(rois(locs{cc}(d)-(tEvent):locs{cc}(d)+(tEvent),cc));
                    elseif (locs{cc}(d)-(tEvent))<=0
                        tempEvent=double(rois(1:locs{cc}(d)+(tEvent),cc));
                    elseif (locs{cc}(d)+(tEvent))>t
                        tempEvent=double(rois(locs{cc}(d)-(tEvent):t,cc));
                    end
%                     figure; plot(tempEvent)
                    
                    R=risetime(tempEvent,'StateLevels',[median(tempEvent), pks{cc}(d)]);
                    riseEvent{dd}=[riseEvent{dd};R];
                    F=falltime(tempEvent,'StateLevels',[median(tempEvent), pks{cc}(d)]);
                    fallEvent{dd}=[fallEvent{dd};F];
%                     if length(R)>1
%                         title(['R=',num2str(R(2)),' F=',num2str(F(2))])
%                     else
%                         title(['R=',num2str(R),' F=',num2str(F)])
%                     end
                    if length(tempEvent)<41
                        tempEvent=[tempEvent;nan(41-length(tempEvent),1)];
                    end
                    eventManip{dd}=[eventManip{dd},tempEvent]; %taking 20 frames as arbitrary
                    %                     elseif locs{cc}(d)-(tEvent*sampRate)<=0
                    %                         eventManip{dd}(d,:)=[nan(abs(LICinfo(a,2)-(tEvent*sampRate)),1);(LIC(1:LICinfo(a,2)+(tEvent*sampRate)+1))];
                    %                     elseif LICinfo(a,2)-(tEvent*sampRate)>0 && (LICinfo(a,2)+(tEvent*sampRate))>length(LIC)
                    %                         eventManip(a,:)=[LIC(LICinfo(a,2)-(tEvent*sampRate):end);nan((LICinfo(a,2)+(tEvent*sampRate))-length(LIC),1)];
                    %                     end
                    tempEvent=[];
                end
            end
        end
    end
end

%% define responsive ROIs
%respROIManip gives several peaks for each situation. But I only need to
%know if responded once
for a=1:manipNum
    respGridROI{a}=unique(respROIManip{a},'stable');
%     pksGridROI{a}=unique(pksROIManip{a},'stable');
end

%% Add to slice Struct
sliceStruct(f).gridROI.auc=auc;
sliceStruct(f).gridROI.Responsive=respGridROI;
% sliceStruct(f).gridROI.Amplitudes=pksROIManip;
sliceStruct(f).gridROI.Amplitudes=pksROIManip;
sliceStruct(f).gridROI.EventTrace=eventManip;
sliceStruct(f).gridROI.RiseTime=riseEvent;%frames
sliceStruct(f).gridROI.FallTime=fallEvent;%frames
sliceStruct(f).gridROI.HW=durPkManip; %frames
%%
save(['sliceSummary.mat'],'sliceStruct','-v7.3');
