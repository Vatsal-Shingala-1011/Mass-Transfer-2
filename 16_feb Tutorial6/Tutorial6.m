clc;
clear;
B = [0.04 0.05 0.07 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.993 1];
C = [0 0.11 0.26 0.375 0.474 0.487 0.468 0.423 0.356 0.274 0.185 0.09 0.001 0];
A = 1 - B - C;

figure(1)
plot(B,C,'bo-');grid on;
hold on;

plot([0 1 0 0],[0 0 1 0], 'k-.', 'linewidth',1.25);
xlabel('Xb');ylabel('Xc');
title('raw data');

tiexc = [0.1 0.245 0.426];
tiexb = [0.048 0.065 0.133];
tieyc = [0.098 0.242 0.409];
tieyb = [0.891 0.736 0.523];

%% Plotting the tie lines
for i = 1:length(tiexc)
plot([tiexb(i) tieyb(i)], [tiexc(i) tieyc(i)],'mo:','linewidth',1.25);
end

%% Slopes of the tie lines 
tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc)
tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i));
end

%% Draw the feed line (FS line)
F = [0, 0.45];
S = [0.96, 0.04];
plot(S,F , 'g', 'linewidth', 0.35)
text(0,  0.45,'F');
text(0.96,  0.04,'S');

p1 = polyfit(B,C,6);
X_=0.04001:0.02:1;
f1 = polyval(p1,X_);
plot(X_,f1,'black')

S = 2500;
F = 2000;

xbf = 0;
xcf = 0.45;
Rny = 0.04; %for Rnx need to find x point from LLE curve
%finding Rnx(ybs) point
ep=0.01;
for i=0:0.01:0.5 
    if polyval(p1, i)<=.04+ep && polyval(f1, i)>=.04-ep
%         disp(i);
        Rnx=i;
        text(i,0.04 ,'Rn');
        break;
    end
end

%finding Mx,My points
M = F + S;
ycs = Rny;
ybs = 1-Rny;
My = (F*xcf + S*ycs)/M; % x coordinate of the mixture point
Mx = (F*xbf + S*ybs)/M; % y coordinate of the mixture point
text(Mx,  My,'M');

% drawing the line RnM
plot([Rnx, Mx], [Rny, My], 'o'); 
hold on;
m = (My - Rny) / (Mx - Rnx);% Calculate the slope (m) and y-intercept (b) of the line y = mx + b
b = Rny - m * Rnx;
x_val = linspace(Rnx, 1, 100); 
y_val = m * x_val + b;
plot(x_val, y_val, 'r');

%find intercept of RnM and LLE curve to get E point
% Extract 
syms x y;
[ybe1,yce1] = vpasolve([y == poly2sym(p1), y == m * x + b],[x,y],[.5 1;0 .7]); %range=[x1 x2;y1 y2]
text(ybe1,  yce1,'E');


%Finding FE and RnS lines
%for FE line F(xbf,xcf) E(ybe,yce)
m2 = (yce1 - xcf) / (ybe1 - xbf);
b2 = xcf - m2 * xbf; %intercept
x_val2 = linspace(0, 1.8, 200); 
y_val2 = m2 * x_val2 + b2;
plot(x_val2, y_val2, 'b');
%for RnS line Rn(xbn,xcn)(Rnx,Rny) S(ybs,ycs)
m3 = (ycs - Rny) / (ybs - Rnx);
b3 = Rny - m3 * Rnx;
x_val3 = linspace(0, 1.8, 200); 
y_val3 = m3 * x_val3 + b3;
plot(x_val3, y_val3, 'b'); 


%Finding delta point by solving FE and RnS lines
syms x y;
eq1 = y == m2*x + b2;
eq2 = y == m3*x + b3;
[Px, Py] = solve([eq1, eq2], [x, y]);
text(Px,  Py,'P');


%after gettig E1 we calculate R1(using tie) and then E2(using R1P line) and so on..
stages=3;
xbr = ones(1,stages);xcr = ones(1,stages);
ybe = ones(1,stages);yce = ones(1,stages);
yce(1)=yce1;ybe(1)=ybe1;
syms x y
%Finding the R1 point by using tie line
%find eqation of tieline
for i=1:stages
    if (0 < yce(i)) && (yce(i) <= 0.098) %setting <yce >yce condition by using graph 
    slope = 0 + (yce(i) - 0) * 0 / 0.098; % 0 
    elseif (0.098 < yce(i)) && (yce(i) <= 0.249)
    slope = tie_slope(1) + (yce(i) - 0.098) * (tie_slope(2) - tie_slope(1)) / (0.249 - 0.098);
    elseif (0.249 < yce(i)) && (yce(i) <= 0.409)
    slope = tie_slope(2) + (yce(i) - 0.249) * (tie_slope(3) - tie_slope(2)) / (0.409 - 0.249);
    end
    
    x_vals = linspace(0, ybe(i), 10);
    y_vals = slope * (x_vals - ybe(i) ) + yce(i);
    intercept = (yce(i) - slope * ybe(i)); %y at x_vals=0
    plot(x_vals, y_vals, 'r', 'linewidth', .05);
    
    %find equation R1 point by intercection of tie line and LLE curve
    syms x y;
    [xbr(i),xcr(i)] = vpasolve([y == poly2sym(p1), y == slope * x + intercept],[x,y],[0 .5;0 .7]); %range=[x1 x2;y1 y2]
    text(xbr(i),xcr(i),['R', num2str(i)]);
    
    
    
    %find R1P line and find point E2 by intersection of LLE & R1P
    %lineR1P   R1(xbr1,xcr1) P(Px,Py)
    m4 = (Py - xcr(i)) / (Px - xbr(i));
    b4 = xcr(i) - m4 * xbr(i);
    x_val4 = linspace(0, Px, 200); 
    y_val4 = m4 * x_val4 + b4;
    plot(x_val4, y_val4, 'b'); 
    %find intercept of R1P and LLE curve to get E2 point 
    [ybe(i+1),yce(i+1)] = vpasolve([y == poly2sym(p1), y == m4 * x + b4],[x,y],[.5 1;0 .7]); %range=[x1 x2;y1 y2]
    text(ybe(i+1),  yce(i+1),['E', num2str(i+1)]);
    
end


%to find number of stages need to plot graph of conc oc c in raffinate atevery stage
figure(2)
plot(xcr, '-bo');
xlabel("Number of stages");
ylabel("Solute Fraction in raffinate");
title("Solute Fraction vs no of stages");

disp("from the graph 0.04 fraction of c remain in raffinat after 10 stages")

