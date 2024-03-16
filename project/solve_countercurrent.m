clc
clear all
B = [ 0.0155 0.017 0.025 0.038 0.06 0.122 0.225 0.354 0.472 0.609 0.735 0.857 0.9533 0.9788]
C = [ 0 0.0285 0.117 0.205 0.262 0.328 0.346 0.336 0.308 0.246 0.173 0.089 0.0187 0.0]
A = 1- B - C


A= 1- B - C;

B1 = zeros(1,10001);
C1 = zeros(1,10001);
for i = 1:10001;
    B1(i) = (i-1)/10000;
    C1(i) =interp1(B,C,B1(i));
end
figure(1)
plot(B1,C1,'b','linewidth',1.00);grid on;
hold on;
plot([0 1 0 0],[0 0 1 0],'k-.','linewidth',1.25)
xlabel('xB');ylabel('xC,yC');title('interpolated data - overall process schematic for 4% DPH')

%fitting the raffinate data
p1 = polyfit(B(1:7),C(1:7),5);
f1= polyval(p1,B(1:7));

%fitting the extract data
p2 = polyfit(B(8:14),C(8:14),6);
f2= polyval(p2,B(8:14));

%tielines
tiexc = [0.0285 0.205 0.328];
tieyc = [0.0187 0.173 0.308];
tiexb = [0.017 0.038 0.122];
tieyb = [0.9533 0.735 0.472];

tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc);
    tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i)) ;
end

%% Tie line drawing for n tie lines where n= length(tiexc)
for i = 1:length(tiexc);
   plot(tiexb(i),tiexc(i),'ro-')
   plot(tieyb(i),tieyc(i),'ro-')
   plot( [tiexb(i)  tieyb(i)],[tiexc(i)  tieyc(i)],'m-.','linewidth',1.0)
end


F = 1000;
S = 1300
xbf =0;
xcf = 0.35;
ybs =1;
ycs = 0;

xcN = 0.03;
index = find(abs(C1-xcN) <= 0.001,1,'first');
xbN = B1(index);
%plotting FS
text(xbf, xcf,'F');
text(ybs,ycs,'S');
plot([xbf,ybs],[xcf,ycs],'LineWidth',1);
slope_FS = (xcf-ycs)/(xbf-ybs);

plot(xbN, xcN, 'ko');
text(xbN, xcN,'RN');

%calculating M 
M = F + S;
My = ((xcf*F) + (ycs*S))/M;
Mx = (F*xbf + S*ybs)/M;
plot(double(Mx),double(My),'ko')
text(double(Mx),double(My),'M')

%slope of RN and M to find E1
slope_RNM = (My - xcN)/(Mx - xbN);
syms x y
eqn1 = y == slope_RNM*(x - xbN) + xcN;
eqn2 = y == poly2sym(p2);
sol12 = vpasolve([eqn1, eqn2], [x, y]);
if(sol12.x(5) < 1 && sol12.y(5) < 1 && sol12.x(5) >= 0 && sol12.y(5) >= 0)
    ybe(1) = sol12.x(5); yce(1) = sol12.y(5);
else
    ybe(1) = sol12.x(6); yce(1) = sol12.y(6);
end
plot(ybe(1), yce(1), 'ko');
text(ybe(1), yce(1)+0.01, 'E1');
plot([xbN, ybe(1)], [xcN, yce(1)])

%slope of F and E1
slope_FE1 = (yce(1) - xcf)/(ybe(1) - xbf);
%slope of RN and S
slope_RNS = (ycs - xcN)/(ybs - xbN);

eqn3 = y == slope_FE1*(x - xbf) + xcf;
eqn4 = y == slope_RNS*(x - xbN) + xcN;

%finding P 
[dx, dy] = vpasolve([eqn3, eqn4], [x, y]);
plot(dx, dy, 'ko');
text(dx + 0.05, dy, 'd');

plot([xbf dx], [xcf, dy], 'k');
plot([xbN(1) dx], [xcN(1), dy], 'k');

%to calculate R1
if ((0 < yce(1)) && (yce(1) <= 0.0187))
   slope = 0 + (yce(1) - 0)*tie_slope(1)/(0.0187);
elseif((0.0187 < yce(1)) && (yce(1) <= 0.173))
   slope = tie_slope(1) + (yce(1) - 0.0187)*(tie_slope(2) - tie_slope(1))/(0.173 - 0.0187);
elseif((0.173 < yce(1)) && (yce(1) <= 0.308))
   slope = tie_slope(2) + (yce(1) - 0.173)*(tie_slope(3) - tie_slope(2))/(0.308 - 0.173);
elseif((yce(1) > 0.308))
   slope = tie_slope(3);
end

eqn5 = y == slope*(x - ybe(1)) + yce(1);

eqn6 = y == poly2sym(p1);
sol12 = vpasolve([eqn5, eqn6]);
if(sol12.x(4) < 1 && sol12.y(4) < 1 && sol12.x(4) >= 0 && sol12.y(4) >= 0)
    xbr(1) = sol12.x(4); xcr(1) = sol12.y(4);
else
    xbr(1) = sol12.x(5); xcr(1) = sol12.y(5);
end
plot([xbr(1), ybe(1)], [xcr(1),yce(1)], 'g',LineWidth=1);
plot(xbr(1), xcr(1), 'ko');
text(xbr(1), xcr(1), 'R1');

i = 1;
while(xcr(i) > xcN)
    plot([xbr(i), dx], [xcr(i), dy], 'k');
    
    slope_Rid = (dy - xcr(i))/(dx - xbr(i));
   
    eqn7 = y == slope_Rid*(x - xbr(i)) + xcr(i);
    eqn8 = y == poly2sym(p2);
    
    [ybe(i+1), yce(i+1)] = vpasolve([eqn7, eqn8], [x, y], [0.8 1; 0.0 0.3])
    plot(ybe(i+1), yce(i+1), 'ko');
    text(ybe(i+1), yce(i+1), strcat('E',num2str(i + 1)));

    if ((0 < yce(i+1)) & (yce(i+1) <= 0.0187))
        slope = 0 + (yce(i+1) - 0)*tie_slope(1)/(0.0187);
    elseif((0.0187 < yce(i+1)) & (yce(i+1) <= 0.173))
        slope = tie_slope(1) + (yce(i+1) - 0.0187)*(tie_slope(2) - tie_slope(1))/(0.173 - 0.0187);
    elseif((0.173 < yce(i+1)) & (yce(i+1) <= 0.308))
        slope = tie_slope(2) + (yce(i+1) - 0.173)*(tie_slope(3) - tie_slope(2))/(0.308- 0.173);
    elseif((yce(i+1) > 0.308))
        slope = tie_slope(3);
    end

    eqn9 = y == slope*(x - ybe(i+1)) + yce(i+1);
    eqn10 = y == poly2sym(p1);
    [xbr(i+1), xcr(i+1)] = vpasolve([eqn9, eqn10], [x, y], [0 0.1; 0.0 0.2]);
    if(xcr(i+1) > xcN)
        plot([xbr(i+1), ybe(i+1)], [xcr(i+1),yce(i+1)], 'g',LineWidth=1);
        plot(xbr(i+1), xcr(i+1), 'ko');
        text(xbr(i+1) - 0.05, xcr(i+1), strcat('R',num2str(i + 1)));
    end
    
    i = i + 1;
end
stages = i

Ybe = ones(1, stages);
Yce = ones(1, stages);
Xbr = ones(1, stages);
Xcr = ones(1, stages);
stage = ones(1, stages);
for j = 1 : stages
    stage(j) = j;
    Xcr(j) = xcr(j);
end
figure(2)
plot(stage, Xcr, 'bo-')
grid on

for j = 1 : stages
    stage(j) = j;
    Yce(j) = yce(j);
end
hold on
plot(stage, Yce, 'mo-')
grid on
title('intermediate stage vs Rxc,Eyc conc profile');grid on;
xlabel('intermediate stage')
ylabel('fraction of acetic acid')
legend('raffinate','extract')

R = ones(1, stages);
E = ones(1, stages);
syms x y
eqn11 = x == F + S - y;
eqn12 = x == (1/Xcr(stages))*((S*ycs) + (F*xcf) - (y*Yce(1)));
[x y] = vpasolve([eqn11, eqn12], [x, y]);
R(stages) = x;
E(1) = y;
frac = (R(stages)*Xcr(stages)) - (S*ycs);
c = R(stages) - S;
for k = 2: stages
syms x y
    eqn13 = (x*Xcr(k-1)) - (y*Yce(k)) == frac;
    eqn14 = x - y == c;
    [R(k-1), E(k)] = vpasolve([eqn13, eqn14], [x, y]);
end
percentageRemoved = ones(1, stages);
for m = 1 : stages 
        percentageRemoved(m) = ((F*xcf) - (R(m)*Xcr(m)))/(F*xcf);
end

figure(3)
plot(stage, percentageRemoved*100, 'bo-')
grid on
title('intermediate stage vs acetic acid removal');grid on;
xlabel('intermediate stage')
ylabel('percentage of acetic acid removed')





