%New generalized path
function [xway,yway] = pathfinder_general(ysk,xsk,skel,wsize);
%     wsize = floor(window/2);
    Local = double(skel(ysk-wsize:ysk+wsize,xsk-wsize:xsk+wsize));
    [row,col,v] = find(Local);
    xway = [];
    yway = [];
    for i = 1:length(col)
        yway(i) = ysk + row(i) - wsize-1;
        xway(i) = xsk + col(i) - wsize-1;
    end
    xway = xway.';
    yway = yway.';
    
end

