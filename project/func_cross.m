function [solute_fraction,stages] = func_cross(F,S,xcf,ycs,stages)
    

%plotting equilateral triangle
ABx = 0:0.1:1;
ABy = zeros(1,length(ABx));
ACx = zeros(1,length(ABx));
ACy = zeros(1,length(ABx));
BCx = zeros(1,length(ABx));
BCy = zeros(1,length(ABx));

figure(1)
plot([0 1],[0 0],'ko-','linewidth',1.5);hold on;
plot([0 0.5],[0 sqrt(3)/2],'ko-','linewidth',1.25)
plot([0.5 1],[sqrt(3)/2 0],'ko-','linewidth',1.25)
title('interpolated data - overall process schematic in equilateral triangle system')
ylim([0 sqrt(3)/2])
xlabel('xB,xB');ylabel('xC,yC');
set(gca,'YtickLabel',ABx(1:1:end))
set(gca,'Ytick',0:0.0866:0.866)

text(0,0,'A')
text(1,0,'B')
text(0.5,sqrt(3)/2,'C')

slope1 = -sqrt(3);
slope2 = sqrt(3);

syms x y
for i = 2:10;    
    [ACx(i),ACy(i)] = vpasolve([y == sqrt(3)*x,y == ABy(i) + slope1*(x - ABx(i))],[x,y]);
    [BCx(i),BCy(i)] = vpasolve([y == -sqrt(3)*(x-1),y == ABy(i) + slope2*(x - ABx(i))],[x,y]);
end

for i = 1:11;
    plot([ABx(i) ACx(i)],[ABy(i) ACy(i)],'o-','color',[0.75,0.75,0.75],'linewidth',0.5);hold on;
    plot([ABx(i) BCx(i)],[ABy(i) BCy(i)],'o-','color',[0.75,0.75,0.75],'linewidth',0.5)
    plot([ACx(i) BCx(12-i)],[ACy(i) BCy(12-i)],'o-','color',[0.75,0.75,0.75],'linewidth',0.5)
end

%% plotting raffinate and extract curves
Br = [0 0.0155 0.017 0.025 0.038 0.06 0.122 0.225];
Be = [0.229 0.354 0.472 0.609 0.735 0.857 0.9533 0.9788 1];

Cr = [0 0 0.0285 0.117 0.205 0.262 0.328 0.346];
Ce = [0.344 0.336 0.308 0.246 0.173 0.089 0.0187 0.0 0];

Bdata = [Br Be];
Cdata = [Cr Ce];
Adata =  Bdata + Cdata;

len = length(Adata);
len1 = length(Br);
Cx = zeros(1,len);
Cy = zeros(1,len);

syms x y
for i = 1:len;    
    [Cx(i),Cy(i)] = vpasolve([y == -sqrt(3)*(x-Adata(i)),y ==  sqrt(3)*(x - Bdata(i))],[x,y]);
end
plot(Cx,Cy,'bo-')

%fitting raffinate curve
p1 = polyfit(Cx(1:len1),Cy(1:len1),3);
f1= polyval(p1,Cx(1:len1));
%fitting extract curve
p2 = polyfit(Cx(len1+1:end),Cy(len1+1:end),2);
f2= polyval(p2,Cx(len1+1:end));

%tielines
tiexc = [0.0285 0.205 0.328];
tiexb = [0.017 0.038 0.122];

tiexa = tiexc + tiexb;

tieyc = [0.0187 0.173 0.308];
tieyb = [0.9533 0.735 0.472];

tieya = tieyc + tieyb;
len_tie = length(tiexa);

C1x = zeros(1,len_tie);
C1y = zeros(1,len_tie);

len_tie = length(tieya);

C2x = zeros(1,len_tie);
C2y = zeros(1,len_tie);

syms x y
for i = 1:len_tie;    
    [C1x(i),C1y(i)] = vpasolve([y == -sqrt(3)*(x-tiexa(i)),y ==  sqrt(3)*(x - tiexb(i))],[x,y]);
    [C2x(i),C2y(i)] = vpasolve([y == -sqrt(3)*(x-tieya(i)),y ==  sqrt(3)*(x - tieyb(i))],[x,y]);
end

%plotting tie lines
for i = 1:length(tiexc);
   plot( [C1x(i)  C2x(i)],[C1y(i)  C2y(i)],'r-.','linewidth',1.0)
end
%calculating tie lines slope
tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc);
    tie_slope(i) = (C2y(i) - C1y(i))/(C2x(i) - C1x(i)) ;
end

% F = 1000;
% S = 1300;
xbf =0;
% xcf = 0.35;
ybs =1;
% ycs = 0;

xcN = 0.03;

B1 = zeros(1,10001);
C1 = zeros(1,10001);
for i = 1:10001;
    B1(i) = (i-1)/10000;
    C1(i) =interp1(Bdata,Cdata,B1(i));
end

index = find(abs(C1-xcN) <= 0.001,1,'first');
xbN = B1(index);

%plotting FS
xaN = xbN + xcN;
xaf = (xbf + xcf);
yas = (ybs + ycs);

[Cfx,Cfy] = vpasolve([y == -sqrt(3)*(x-xaf),y ==  sqrt(3)*(x - xbf)],[x,y]);
[Csx,Csy] = vpasolve([y == -sqrt(3)*(x-yas),y ==  sqrt(3)*(x - ybs)],[x,y]);
[Cnx,Cny] = vpasolve([y == -sqrt(3)*(x-xaN),y ==  sqrt(3)*(x - xbN)],[x,y]);

text(Cfx,Cfy,'F');
text(Csx,Csy-0.02,'S');
plot(Cfx,Cfy,'ro');
plot(Csx,Csy,'ro');
plot([Cfx,Csx],[Cfy,Csy],'LineWidth',1);
%slope of F AND S
slope_FS = (Cfy - Csy)/(Cfx-Csx);

plot(Cnx, Cny, 'ko');
text(Cnx, Cny,'RN');
% stages = 3;
% %calculating M 
% M = F + S;
% My = ((xcf*F) + (ycs*S))/M;
% Mx = (F*xbf + S*ybs)/M;
% 
% Ma = Mx+ My;
% [Cmx,Cmy] = vpasolve([y == -sqrt(3)*(x-Ma),y ==  sqrt(3)*(x - Mx)],[x,y]);
% 
% %plotting M
% plot(double(Cmx),double(Cmy),'ko');
% text(double(Cmx),double(Cmy),'M');
% 
% %slope of RN and M to find E1
% % slope_RNM = (Cmy - Cny)/(Cmx - Cnx);
% 
% 
% %to calculate R1
% if ((0 < Mx) && (Mx <= 0.0162))
%    slope = 0 + (Mx - 0)*tie_slope(1)/(0.0162);
% elseif((0.0162 < Mx) && (Mx <= 0.1498))
%    slope = tie_slope(1) + (Mx - 0.0162)*(tie_slope(2) - tie_slope(1))/(0.1498 - 0.0162);
% elseif((0.1498 < Mx) && (Mx <= 0.2667))
%    slope = tie_slope(2) + (Mx - 0.1498)*(tie_slope(3) - tie_slope(2))/(0.2667 - 0.1498);
% elseif((Mx > 0.2667))
%    slope = tie_slope(3);
% end
% 
% syms x y
% eqn1 = y == slope*(x - Mx) + My;
% eqn2 = y == poly2sym(p2);
% [ybe(1), yce(1)] = vpasolve([eqn1, eqn2], [x, y],[0.4 1;0 0.4]);
% 
% syms x y
% eqn2 = y == poly2sym(p1);
% [xbr(1), xcr(1)] = vpasolve([eqn1, eqn2], [x, y],[0 0.4;0 0.4]);
% 
% plot(ybe(1), yce(1), 'ko');
% text(ybe(1), yce(1)+0.01, 'E1');
% plot(xbr(1), xcr(1), 'ko');
% text(xbr(1), xcr(1)+0.01, 'R1');
% plot([xbr, ybe(1)], [xcr, yce(1)],'m-');
% 
% %slope of S and R1
% slope_SR1 = (xcr(1) - ycs)/(xbr(1) - ybs);

% %slope of RN and S
% slope_RNS = (Cny - Csy)/(Cnx - Csx);
% 
% eqn3 = y == slope_FE1*(x - Cfx) + Cfy;
% eqn4 = y == slope_RNS*(x - Cnx) + Cny;
% 
% %finding P >> delta point
% [dx, dy] = vpasolve([eqn3, eqn4], [x, y]);
% plot(dx, dy, 'ko');
% text(dx + 0.05, dy, 'd');
% 
% plot([Cfx dx], [Cfy, dy], 'k');
% plot([Cnx dx], [Cny, dy], 'k');
% 
% %to calculate R1
% if ((0 < yce(1)) && (yce(1) <= 0.0162))
%    slope = 0 + (yce(1) - 0)*tie_slope(1)/(0.0162);
% elseif((0.0162 < yce(1)) && (yce(1) <= 0.1498))
%    slope = tie_slope(1) + (yce(1) - 0.0162)*(tie_slope(2) - tie_slope(1))/(0.1498 - 0.0162);
% elseif((0.1498 < yce(1)) && (yce(1) <= 0.2667))
%    slope = tie_slope(2) + (yce(1) - 0.1498)*(tie_slope(3) - tie_slope(2))/(0.2667 - 0.1498);
% elseif((yce(1) > 0.2667))
%    slope = tie_slope(3);
% end
% 
% eqn5 = y == slope*(x - ybe(1)) + yce(1);
% 
% eqn6 = y == poly2sym(p1);
% [xbr(1) xcr(1)] = vpasolve([eqn5, eqn6], [x, y],[0 0.4;0 0.4]);
% 
% plot([xbr(1), ybe(1)], [xcr(1),yce(1)], 'm',LineWidth=1);
% plot(xbr(1), xcr(1), 'ko');
% text(xbr(1), xcr(1), 'R1');
% 
% %finding number of stages
% i = 1;
% while(xcr(i) > xcN)
%     plot([xbr(i), dx], [xcr(i), dy], 'k');
%     
%     %slope of R(i) and delta pont to get E(i+1)
%     slope_Rid = (dy - xcr(i))/(dx - xbr(i));
%    
%     %solving R(i)_d and extract curve to get E(i+1)
%     eqn7 = y == slope_Rid*(x - xbr(i)) + xcr(i);
%     eqn8 = y == poly2sym(p2);
%     
%     [ybe(i+1), yce(i+1)] = vpasolve([eqn7, eqn8], [x, y], [0.4 1; 0.0 0.4]);
%     plot(ybe(i+1), yce(i+1), 'ko');
%     text(ybe(i+1), yce(i+1), strcat('E',num2str(i + 1)));
% 
%     %interpolation of tieline wrt E(i+1) to find R(i+1)
%     if ((0 < yce(i+1)) && (yce(i+1) <= 0.0162))
%        slope = 0 + (yce(i+1) - 0)*tie_slope(1)/(0.0162);
%     elseif((0.0162 < yce(i+1)) && (yce(i+1) <= 0.1498))
%        slope = tie_slope(1) + (yce(i+1) - 0.0162)*(tie_slope(2) - tie_slope(1))/(0.1498 - 0.0162);
%     elseif((0.1498 < yce(i+1)) && (yce(i+1) <= 0.2667))
%        slope = tie_slope(2) + (yce(i+1) - 0.1498)*(tie_slope(3) - tie_slope(2))/(0.2667 - 0.1498);
%     elseif((yce(i+1) > 0.2667))
%        slope = tie_slope(3);
%     end
%     %%solving eqn of line passing thriugh E(i+1) and Raffinate curve to get R(i+1)
%    
% 
%     eqn9 = y == slope*(x - ybe(i+1)) + yce(i+1);
%     eqn10 = y == poly2sym(p1);
%     [xbr(i+1), xcr(i+1)] = vpasolve([eqn9, eqn10], [x, y], [0 0.1; 0.0 0.2]);
%     if(xcr(i+1) > Cny)
%         plot([xbr(i+1), ybe(i+1)], [xcr(i+1),yce(i+1)], 'm',LineWidth=1);
%         plot(xbr(i+1), xcr(i+1), 'ko');
%         text(xbr(i+1) - 0.03, xcr(i+1), strcat('R',num2str(i + 1)));
%     end
%     
%     i = i + 1;
% end
% stages = i;
% 
% Ybe = ones(1, stages);
% Yce = ones(1, stages);
% Xbr = ones(1, stages);
% Xcr = ones(1, stages);
% stage = ones(1, stages);
% for j = 1 : stages
%     stage(j) = j;
%     Xcr(j) = xcr(j);
% end
% figure(2)
% plot(stage, Xcr, 'bo-')
% grid on
% 
% for j = 1 : stages
%     stage(j) = j;
%     Yce(j) = yce(j);
% end
% hold on
% plot(stage, Yce, 'mo-')
% grid on
% title('intermediate stage vs Rxc,Eyc conc profile');grid on;
% xlabel('intermediate stage')
% ylabel('fraction of acetic acid')
% legend('raffinate','extract')
% 
% R = ones(1, stages);
% E = ones(1, stages);
% syms x y
% eqn11 = x == F + S - y;
% eqn12 = x == (1/Xcr(stages))*((S*ycs) + (F*xcf) - (y*Yce(1)));
% [x y] = vpasolve([eqn11, eqn12], [x, y]);
% R(stages) = x;
% E(1) = y;
% frac = (R(stages)*Xcr(stages)) - (S*ycs);
% c = R(stages) - S;
% for k = 2: stages
% syms x y
%     eqn13 = (x*Xcr(k-1)) - (y*Yce(k)) == frac;
%     eqn14 = x - y == c;
%     [R(k-1), E(k)] = vpasolve([eqn13, eqn14], [x, y]);
% end
% percentageRemoved = ones(1, stages);
% for m = 1 : stages 
%         percentageRemoved(m) = ((F*xcf) - (R(m)*Xcr(m)))/(F*xcf);
% end
% 
% figure(3)
% plot(stage, percentageRemoved*100, 'bo-')
% grid on
% title('intermediate stage vs acetic acid removal');grid on;
% xlabel('intermediate stage')
% ylabel('percentage of acetic acid removed')
% 


%%To be calculated-----------------
r = ones(1,20);
e = ones(1,20);
xbr = ones(1,20);
xcr = ones(1,20);
ybe = ones(1,20);
yce = ones(1,20);
M = ones(1, 20);
Mx = ones(1,20);
My = ones(1, 20);
i=1;

M(i) = F + S;
Mx(i) = (F*xbf + S*ybs)/M(i);
My(i) = (F*xcf + S*ycs)/M(i);
plot(Mx(i),My(i),'o','MarkerSize',5,'MarkerFaceColor','yellow')
text(Mx(i),My(i),'M')

if ((0 < My(i)) && (My(i) <= tiexc(1)))
    slope = 0 + (My(i)-0)/(tiexc(1)-0)*tie_slope(1);
elseif((tiexc(1) < My(i)) && (My(i) <= tiexc(2)))
    slope= tie_slope(1) + (My(i)-tiexc(1))/(tiexc(2)-tiexc(1))*tie_slope(2);
elseif((tiexc(2) < My(i)) && (My(i) <= tiexc(3)))
    slope= tie_slope(2) + (My(i)-tiexc(2))/(tiexc(3)-tiexc(2))*tie_slope(3);
elseif((tiexc(3) < My(i)) && (My(i) <= tiexc(4)))
    slope= tie_slope(3) + (My(i)-tiexc(3))/(tiexc(4)-tiexc(3))*tie_slope(4);
elseif((tiexc(4) < My(i)) && (My(i) <= tiexc(5)))
    slope= tie_slope(4) + (My(i)-tiexc(4))/(tiexc(5)-tiexc(4))*tie_slope(5);
elseif((tiexc(5) < My(i)) && (My(i) <= tiexc(6)))
    slope= tie_slope(5) + (My(i)-tiexc(5))/(tiexc(6)-tiexc(5))*tie_slope(6);
elseif((tiexc(6) < My(i)) && (My(i) <= tiexc(7)))
    slope= tie_slope(6) + (My(i)-tiexc(6))/(tiexc(7)-tiexc(6))*tie_slope(7);
elseif((My(i) > tiexc(7)))
    slope= tie_slope(7) + (My(i)-tiexc(7))/(0.5-tiexc(7))*1;
end

%%tie line and the polynomial solving together
syms x y
[xbr(i),xcr(i)] = vpasolve([y==poly2sym(p1), y==(My(i) +slope*(x - Mx(i)) )],[x,y],[0 0.24;0 0.346]);
plot(xbr(i),xcr(i),'o','MarkerSize',5,'MarkerFaceColor','yellow',Color='black');
text(xbr(i),xcr(i),'R');
syms x y
[ybe(i),yce(i)] = vpasolve([y == poly2sym(p2),y == My(i) + slope*(x - Mx(i))],[x,y],[0.24 1; 0 0.346]);
plot(ybe(i),yce(i),'o','MarkerSize',5,'MarkerFaceColor','yellow',Color='black')
text(ybe(i),yce(i),'E')

plot([xbr(i) ybe(i)],[xcr(i) yce(i)],Color='black',LineWidth=1.0) %plotting the the required tie line.

R = [xbr(i) xcr(i)];
E = [ybe(i) yce(i)];

syms E R
[e(i),r(i)] = vpasolve(E == M(i)-R, E == R*((My(i)-xcr(i))/(yce(i)-My(i))),[E,R]); %calculation of E and R using equation 8.3 and 8.7 BKD 
disp(r(i))
disp(e(i))

% solute_fraction = (yce(i)*e - ycs*S)/(F*xcf); %fraction of solute saparated.

F = r(i);
xbf = xbr(i);
xcf = xcr(i);
plot([xbf ybs],[xcf ycs])

%%finding of the of required tie line-------------------
for i = 2:stages
    M(i) = F + S;
    Mx(i) = (F*xbf + S*ybs)/M(i);
    My(i) = (F*xcf + S*ycs)/M(i);
    plot(Mx(i),My(i),'o','MarkerSize',5,'MarkerFaceColor','yellow')
    text(Mx(i),My(i),'M')

    if ((0 < My(i)) && (My(i) <= tiexc(1)))
        slope = 0 + (My(i)-0)/(tiexc(1)-0)*tie_slope(1);
    elseif((tiexc(1) < My(i)) && (My(i) <= tiexc(2)))
        slope= tie_slope(1) + (My(i)-tiexc(1))/(tiexc(2)-tiexc(1))*tie_slope(2);
    elseif((tiexc(2) < My(i)) && (My(i) <= tiexc(3)))
        slope= tie_slope(2) + (My(i)-tiexc(2))/(tiexc(3)-tiexc(2))*tie_slope(3);
    elseif((tiexc(3) < My(i)) && (My(i) <= tiexc(4)))
        slope= tie_slope(3) + (My(i)-tiexc(3))/(tiexc(4)-tiexc(3))*tie_slope(4);
    elseif((tiexc(4) < My(i)) && (My(i) <= tiexc(5)))
        slope= tie_slope(4) + (My(i)-tiexc(4))/(tiexc(5)-tiexc(4))*tie_slope(5);
    elseif((tiexc(5) < My(i)) && (My(i) <= tiexc(6)))
        slope= tie_slope(5) + (My(i)-tiexc(5))/(tiexc(6)-tiexc(5))*tie_slope(6);
    elseif((tiexc(6) < My(i)) && (My(i) <= tiexc(7)))
        slope= tie_slope(6) + (My(i)-tiexc(6))/(tiexc(7)-tiexc(6))*tie_slope(7);
    elseif((My(i) > tiexc(7)))
        slope= tie_slope(7) + (My(i)-tiexc(7))/(0.5-tiexc(7))*1;
    end

    %%tie line and the polynomial solving together
    syms x y
    [xbr(i),xcr(i)] = vpasolve([y==poly2sym(p1), y==(My(i) +slope*(x - Mx(i)) )],[x,y],[0 0.24;0 0.346]);
    plot(xbr(i),xcr(i),'o','MarkerSize',5,'MarkerFaceColor','yellow',Color='black');
    text(xbr(i),xcr(i),'R');
    syms x y
    [ybe(i),yce(i)] = vpasolve([y == poly2sym(p2),y == My(i) + slope*(x - Mx(i))],[x,y],[0.24 1; 0 0.346]);
    plot(ybe(i),yce(i),'o','MarkerSize',5,'MarkerFaceColor','yellow',Color='black')
    text(ybe(i),yce(i),'E')

    plot([xbr(i) ybe(i)],[xcr(i) yce(i)],Color='black',LineWidth=1.0) %plotting the the required tie line.

    R = [xbr(i) xcr(i)];
    E = [ybe(i) yce(i)];

    syms E R
    [e(i),r(i)] = vpasolve(E == M(i)-R, E == R*((My(i)-xcr(i))/(yce(i)-My(i))),[E,R]); %calculation of E and R using equation 8.3 and 8.7 BKD 
    disp(r(i))
    disp(e(i))

    

    F = r(i);
    xbf = xbr(i);
    xcf = xcr(i);
    plot([xbf ybs],[xcf ycs])
    
    %plot([double(xbr(i)) double(ybe(i))],[double(xcr(i)) double(yce(i))],'o-','Color',[0,0.25,1],'linewidth',0.35)
    %plot(double(Mx(i)),double(My(i)),'bo')

end
fprintf('The minimum number of stages is %d\n',i);

ax = gca;
exportgraphics(ax,'plot_cross.png','Resolution',720)

solute_fraction = (yce(i)*e - ycs*S)/(F*xcf); %fraction of solute saparated.



try
    writematrix([solute_fraction stages],"out_cross.xlsx");
%     writecell(count-2,"out_counter.xlsx");
    disp("Done")
    catch
    disp("Can't write excel")
end

end