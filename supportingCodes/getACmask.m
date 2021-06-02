%Kellner et al., 2021
function [ACmask] = getACmask(X)
    % Select ROI for AC dF/Fo analysis

    Xmean = mean(X,3);
    [m,n] = size(Xmean);
    
    h = figure;
    h_im = imagesc(Xmean);
     
    AC = drawcircle('Center',[270 320],'Radius',150,'Color','r');
    setResizable(AC,0);
    wait(AC);
    ACmask = createMask(AC, h_im);
    
end