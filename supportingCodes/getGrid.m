%Kellner et al., 2021
function [indices] = getGrid(imgWidth,imgHeight,sizeGrid)
    index = 1;
    indices = [];
    for i = 1:sizeGrid:imgWidth
        for j = 1:sizeGrid:imgHeight
            indices(index,:) = [i j i+sizeGrid-1 j i+sizeGrid-1 j+sizeGrid-1 i j+sizeGrid-1];
            index = index + 1;
        end
    end
end