function [xline,pline,m1] = OrthoLine2(xway,yway,UB);
xm0 = mean(xway);
ym0 = mean(yway);
% ex = length(unique(xway));
% verticality = std(yway)/std(xway);
% if abs(verticality) > 1.5
%     xwayv = yway;
%     ywayv= xway;
%     xm0v = mean(xwayv);
%     ym0v = mean(ywayv);
%     [p,S] = polyfit(xwayv,ywayv,2);
% % if p(1) ~= inf
% m = atan(p(1)*2*xm0v+p(2));
% if verticality > 0
% PLm = -m;
% elseif verticality < 0
%     PLm = m;
% end%tan( m + pi/2); 
% m1 = atan(m) + pi/2;
% PLc = ym0 - PLm*xm0;
% xlinev = linspace(xm0-UB,ym0+UB);
% plinev = PLm*xlinev + PLc;
% xline = round(xlinev);
% pline = round(plinev);
% %     [p,S] = polyfit(yway,xway,1);
% %    xm0 = mean(yway);
% % ym0 = mean(xway);
% %     m = atan(p(1));
% % PLm = 1/m;%tan( m +pi/2); 
% % m1 = m+pi/2;
% % PLc = xm0 - PLm*ym0;
% % PLc = -PLc/m;
% % xline = linspace(ym0-UB,ym0+UB);
% % pline = PLm*xline - PLc;
% % xline = round(xline);
% % pline = round(pline);
% %     xline = xm0*ones(1,round(2*UB));
% %     pline = ym0*linspace(xm0-UB,xm0+UB,length(xline));
% %     xline = round(xline);
% %     pline = round(pline);
% else

% [p] = polyfit(xway,yway,2);
% % if p(1) ~= inf
% m = atan(p(1)*2*xm0+p(2));
 [p] = linortfit2(xway, yway);
 m = atan(p(1));
PLm = tan( m +pi/2); 
m1 = m;
PLc = ym0 - PLm*xm0;
xline = linspace(xm0-UB,xm0+UB,150);
pline = PLm*xline + PLc;
xline = round(xline);
pline = round(pline);
if m1 == 0
    xline = round(xm0*ones(1,150));
    pline = round(linspace(ym0-UB,ym0+UB,150));
end
% end
% end
% end
% plot(xline,pline);
end
