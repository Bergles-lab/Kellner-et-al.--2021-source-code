%Kellner et al., 2021
function [dFoF, Fo] = normalizeImg(img, percentile)
% function [dFoF, Fo] = normalizeImg(img)
% %normalizeImg Normalizes image based on percentile chosen after bleach correction. 
% %   On a pixel by pixel basis, Fo is created by taking the pixel value at
% %   the xth percentile. This is subtracted off of the original image; the
% %   resulting image is then divided by Fo. 

% %     sampRate = 10; %sampling rate in Hz
    [m,n,T] = size(img);
% %     imgBC = bleachCorrect(img,sampRate);
% %     clear img;
% %     disp('Bleach correction finished. Subtracting baseline...');
    
    %%Normalize by taking Xth percentile
    Xreshape = reshape(img,m*n,T);
    Fo = prctile(double(Xreshape),percentile,2);
    Fo = reshape(Fo,m,n);
%     Fo=squeeze(median((img),3)); %normalize by median
    dFoF = ((img) - Fo) ./ (Fo);
end


