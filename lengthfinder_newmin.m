
function [x1,y1,x2,y2,len] = lengthfinder_newmin(xline,pline,xsk,ysk,bx,by,bw,UB);
%     len1 = 0;
%     len2 = 0;
    count = 0;
    x1 = 0;
    y1 = 0;
    x2 = 0;
    y2 =0;
    possibleDist = [];
    for line = floor(length(xline)/2):-1:1
        
        xval = xline(line);
        yval = pline(line);
        if (xval ~= 1 && yval ~= 1 && (yval < size(bw,1)) && (xval <size(bw,2)))
        for bound = 1:length(bx)
            xbound = bx(bound);
            ybound = by(bound);
           
            
            if (xbound ~= 1 && ybound ~=1 && (ybound < size(bw,1)) && (xbound <size(bw,2)))
            if (xval == xbound && yval == ybound)
               
                count = count + 1;
%                 len1 = len1 + sqrt((xbound-xsk)^2+(ybound-ysk)^2);
%                 if count ==1
%                     x1 = x1+xval;
%                     y1 = y1+yval;
                    possibleDist(count,1) = xval;
                    possibleDist(count,2) = yval;
                    possibleDist(count,3) = sqrt((xval-xsk).^2+(yval-ysk).^2);
                    possibleDist = possibleDist(any(possibleDist,2),:);
                   
%                 end
                
            end
            

        end
        end
        end
    end
    if ~isempty(possibleDist)
     [minDistance, indexOfMin] = min(possibleDist(:,3));
     if minDistance < UB
                    x1 = possibleDist(indexOfMin,1);
                    y1 = possibleDist(indexOfMin,2);
     end
    end
    count = 0;
    possibleDist = [];
    for line = floor(length(xline)/2):length(xline)
        xval = xline(line);
        yval = pline(line);
        
        if (xval ~= 1 && yval ~= 1 && (yval < size(bw,1)) && (xval <size(bw,2)))
        for bound = 1:length(bx)
            xbound = bx(bound);
            ybound = by(bound);
            if (xbound ~= 1 && ybound ~=1 && (ybound < size(bw,1)) && (xbound <size(bw,2)))
            if (xval == xbound && yval == ybound)
               
                count = count + 1;
%                 len2 = len2 + sqrt((xbound-xsk)^2+(ybound-ysk)^2);
%                 if count ==1
%                     x2 = x2+xval;
%                     y2 = y2+yval;
% 
%                 end
                    possibleDist(count,1) = xval;
                    possibleDist(count,2) = yval;
                    possibleDist(count,3) = sqrt((xval-xsk).^2+(yval-ysk).^2);
                    possibleDist = possibleDist(any(possibleDist,2),:);
                   
                
            end
            

        end
        end
        end
    end
    if ~isempty(possibleDist)
     [minDistance, indexOfMin] = min(possibleDist(:,3));
     if minDistance < UB
                    x2 = possibleDist(indexOfMin,1);
                    y2 = possibleDist(indexOfMin,2);
     end
    end
%     len = len1 + len2;
    len = sqrt((x1-x2).^2+(y1-y2).^2);
end