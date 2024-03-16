%tutorial4_e
%cross->counter
clc;
clear;

%dufi
B = [0 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0021 0.0031 0.0058 0.0136 0.0372 0.05 0.08 0.09 0.102 0.11  0.1259 0.1376 0.1508 0.16 0.3357 0.4751 0.608 0.7698 0.8872 0.9982 1];
C = [0 0 0.005 0.01 .015 .02 .03 .04 0.045 .05 0.052 .055 .06 .065 .07 .09  0.1 0.2 0.3 0.4 0.5 0.5137 0.546 0.557 0.57 0.5821 0.6 0.6058 0.6107 0.6093 0.5819 0.4944 0.3748 0.2223 0.1079 0 0];
figure(1)
plot(B,C,'bo-');
grid on;
hold on;
plot([0 1 0 0],[0 0 1 0], 'k-.', 'linewidth',1.25);
xlabel("xB, yB");
ylabel("xC, yC");
title('raw data');

tiexc = [0, 0.2, 0.4, 0.6];
tieyc = [0.0, 0.2223, 0.4944, 0.6107];
tiexb = [0.0011, 0.0031, 0.0136, 0.1259];
tieyb = [0.9982, 0.7698, 0.4751, 0.1508];

for i = 1:length(tiexc)
plot([tiexb(i) tieyb(i)], [tiexc(i) tieyc(i)],'mo:','linewidth',1.25);
end

tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc)
tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i));
end

F = [0, 0.5];
S = [0.995, 0.005];
plot(S,F,'g','linewidth',0.35)
text(0, 0.5,'F');
text(0.995, 0.005,'S');

p1 = polyfit(B(1:28),C(1:28),3);
p2 = polyfit(B(28:37),C(28:37),3);
xbr = [0: 0.001: 0.17];
xcr = [0.17:0.001:1];
f1 = polyval(p1, xbr);
f2 = polyval(p2, xcr);

plot(xbr, f1, 'r');
hold on;
plot(xcr, f2, 'r');

first=1;

stages = 20;
filename = 'Counter.xlsx';
% for F =400:5:2000
f=1000;  
Feed=F;
        S =F;
        xbf = 0;
        xcf = 0.5;
        xc_in_feed=xcf;
        ybs = 0.995;
        ycs = 0.005;
        
        r = ones(1,stages);
        e = ones(1,stages);
        xbr = ones(1,stages);xcr = ones(1,stages);
        ybe = ones(1,stages);yce = ones(1,stages);
        M = ones(1,stages);
        Mx = ones(1,stages);
        My = ones(1,stages);
        syms R E
        syms x y
        
%%%%%
        %finding Mx,My points
        M = F + S;
        My = (F*xcf + S*ycs)/M; % y coordinate of the mixture point
        Mx = (F*xbf + S*ybs)/M; % x coordinate of the mixture point
        text(double(Mx),double(My),'M');
        plot(double(Mx),double(My),'bo');
              
%         for Rny = 0.02 :0.001 : 0.2
            Rny=0.04;
            Rnx=0.04;%const
            % drawing the line RnM
            plot([Rnx, Mx], [Rny, My], 'o'); 
            hold on;
            m = (My - Rny) / (Mx - Rnx);% Calculate the slope (m) and y-intercept (b) of the line y = mx + b
            b = Rny - m * Rnx;
            x_val = linspace(Rnx, 1, 100); 
            y_val = m * x_val + b;
            plot(x_val, y_val, 'r');
% 
% 
%             %find intercept of RnM and LLE curve to get E point
%             % Extract 
%             syms x y;
%             [ybe1,yce1] = vpasolve([y == poly2sym(p2), y == m * x + b],[x,y],[.5 1;0 .7]); %range=[x1 x2;y1 y2]
%             text(ybe1,  yce1,'E');
%             
%             %Finding FE and RnS lines
%             %for FE line F(xbf,xcf) E(ybe,yce)
%             m2 = (yce1 - xcf) / (ybe1 - xbf);
%             b2 = xcf - m2 * xbf; %intercept
%             x_val2 = linspace(0, 1.8, 200); 
%             y_val2 = m2 * x_val2 + b2;
%             plot(x_val2, y_val2, 'b');
%             %for RnS line Rn(xbn,xcn)(Rnx,Rny) S(ybs,ycs)
%             m3 = (ycs - Rny) / (ybs - Rnx);
%             b3 = Rny - m3 * Rnx;
%             x_val3 = linspace(0, 1.8, 200); 
%             y_val3 = m3 * x_val3 + b3;
%             plot(x_val3, y_val3, 'b'); 
% 
%             %Finding delta point by solving FE and RnS lines
%             syms x y;
%             eq1 = y == m2*x + b2;
%             eq2 = y == m3*x + b3;
%             [Px, Py] = solve([eq1, eq2], [x, y]);
%             text(Px,  Py,'P');
% 
%             %after gettig E1 we calculate R1(using tie) and then E2(using R1P line) and so on..
%             stages=3;
%             xbr = ones(1,stages);xcr = ones(1,stages);
%             ybe = ones(1,stages);yce = ones(1,stages);
%             yce(1)=yce1;ybe(1)=ybe1;
%             syms x y
%             %Finding the R1 point by using tie line
%             %find eqation of tieline
%             for i=1:stages
%                 if (0 < yce(i)) && (yce(i) <= 0.2) %setting <yce >yce condition by using graph 
%                 slope = 0 + (yce(i) - 0) * 0 / 0.2; % 0 
%                 elseif (0.2 < yce(i)) && (yce(i) <= 0.4)
%                 slope = tie_slope(1) + (yce(i) - 0.2) * (tie_slope(2) - tie_slope(1)) / (0.4 - 0.2);
%                 elseif (0.4 < yce(i)) && (yce(i) <= 0.6)
%                 slope = tie_slope(2) + (yce(i) - 0.4) * (tie_slope(3) - tie_slope(2)) / (0.6 - 0.4);
%                 end
%                 
%                 x_vals = linspace(0, ybe(i), 10);
%                 y_vals = slope * (x_vals - ybe(i) ) + yce(i);
%                 intercept = (yce(i) - slope * ybe(i)); %y at x_vals=0
%                 plot(x_vals, y_vals, 'r', 'linewidth', .05);
%                 
%                 %find equation R1 point by intercection of tie line and LLE curve
%                 syms x y;
%                 [xbr(i),xcr(i)] = vpasolve([y == poly2sym(p1), y == slope * x + intercept],[x,y],[0 .5;0 .7]); %range=[x1 x2;y1 y2]
%                 text(xbr(i),xcr(i),['R', num2str(i)]);
%                 
%                 %find R1P line and find point E2 by intersection of LLE & R1P
%                 %lineR1P   R1(xbr1,xcr1) P(Px,Py)
%                 m4 = (Py - xcr(i)) / (Px - xbr(i));
%                 b4 = xcr(i) - m4 * xbr(i);
%                 x_val4 = linspace(0, Px, 200); 
%                 y_val4 = m4 * x_val4 + b4;
%                 plot(x_val4, y_val4, 'b'); 
%                 %find intercept of R1P and LLE curve to get E2 point 
%                 [ybe(i+1),yce(i+1)] = vpasolve([y == poly2sym(p1), y == m4 * x + b4],[x,y],[.5 1;0 .7]); %range=[x1 x2;y1 y2]
%                 text(ybe(i+1),  yce(i+1),['E', num2str(i+1)]);
%                 
%             end


       