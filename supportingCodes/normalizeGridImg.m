%Kellner et al., 2021
function [dFoF, Fo, img] = normalizeGridImg(img, percentile, positiveIndices)
%normalizeImg Normalizes image based on percentile chosen after bleach correction. 
%   On a block by block basis, Fo is created by taking the pixel value at
%   the xth percentile. This is subtracted off of the original image; the
%   resulting image is then divided by Fo. 

    [m,n,T] = size(img);
    
    for i = 1:size(positiveIndices,1)
        if positiveIndices(i,6)>m
            xs = [positiveIndices(i,2):m];
        else
            xs = [positiveIndices(i,2):positiveIndices(i,6)];
        end
        if positiveIndices(i,3)>n
            ys = [positiveIndices(i,1):n];
        else
            ys = [positiveIndices(i,1):positiveIndices(i,3)];
        end
        temp = squeeze(mean(mean(img(xs,ys,:),2),1));
        rois(:,i) = temp;
    end
    
    %Normalize by taking Xth percentile
    Fo = prctile(rois,percentile,1);
%     Fo=squeeze(median((rois),1)); %normalize by median
    dFoF = (single(rois) - single(Fo)) ./ single(Fo);
    
    %correct any baseline drift
    for i = 1:size(dFoF,2)
        dFoF(:,i) = smooth(msbackadj([1:T]',dFoF(:,i),'WindowSize',45,'StepSize',45,'Showplot',0));
    end
end


