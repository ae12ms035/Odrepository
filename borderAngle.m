%calcutes the orientation of the ODC borders with respect to the V1/v2
%border
function [borderAng] = borderAngle(bw, odBorder, wsize, Theta0)
CC = bwconncomp(bw);
ObjNum = CC.NumObjects;
[B,L,n,A] = bwboundaries(bw,'holes');
borderTheta0 = Theta0;
ct = 0;
for ob = 1:ObjNum
    enclosed_boundaries = find(A(:,ob));
                   
    BB = B{ob,1};
    if ~isempty(enclosed_boundaries)
                        
        cat = vertcat(BB,B{enclosed_boundaries(:),1});
        BB = cat;
    end
    by = BB(:,1);            
    bx = BB(:,2);
    for k = 1:length(bx)
        yBorder = by(k);
        xBorder = bx(k);
    
    if (xBorder>wsize && yBorder>wsize && yBorder <size(bw,1)-wsize && xBorder < size(bw,2)-wsize)
                ct = ct+1;
                [xway,yway] = pathfinder_general(yBorder,xBorder,odBorder,wsize);
                if length(xway)>2 
                 [p] = linortfit2(xway, yway);
                  m1 = atan(p(1));
                
%                 fig2 = figure(30);
               
                borderTheta(ct,1) = yBorder;
                borderTheta(ct,2) = xBorder;
                invtan = -m1;%-
                if invtan < 0
                    borderTheta(ct,3) = pi + invtan-borderTheta0;% i-(atan(m1)*180/pi);%-borderThetaV12(n1,2));
                     if borderTheta(ct,3) < 0
                        borderTheta(ct,3) = pi + borderTheta(ct,3);
                    end
                else
                    borderTheta(ct,3) = invtan-borderTheta0;%180-(atan(m1)*180/pi);%-borderThetaV12(n1,2));
                     if borderTheta(ct,3) < 0
                        borderTheta(ct,3) = pi + borderTheta(ct,3);
                    end
                end
                end
    end
    end
end
borderAng = borderTheta(any(borderTheta,2),:);
end