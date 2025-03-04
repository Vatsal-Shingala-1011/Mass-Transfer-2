clc;
clear;
B = [0 0.0011 0.0021 0.0031 0.0058 0.0136 0.0372 0.11 0.1259 0.1376 0.1508 0.16 0.3357 0.4751 0.608 0.7698 0.8872 0.9982 1];
C = [0 0 0.1 0.2 0.3 0.4 0.5 0.5821 0.6 0.6058 0.6107 0.6093 0.5819 0.4944 0.3748 0.2223 0.1079 0 0];
figure(1)
plot(B,C,'bo-');
grid on;
hold on;
plot([0 1],[1 0],'k-','linewidth',1.25)
plot([0 0],[1 0],'k-','linewidth',1.25)
plot([0 1],[0 0],'k-','linewidth',1.25)
xlabel("xB, yB");
ylabel("xC, yC");
title("raw data");

tiexc = [0.0 0.2 0.4 0.6];
tieyc = [0.0 0.2223 0.4944 0.6107];
tiexb = [0.0011 0.0031 0.0136 0.1259];
tieyb = [0.9982 0.7698 0.4751 0.1508];

for i = 1:length(tiexc)
plot([tiexb(i) tieyb(i)], [tiexc(i) tieyc(i)], 'mo:', 'linewidth', 1.25);
end

tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc)
    tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i));
end

F=1000;
S=800;
xcf=0.5;
xbf=0;
ybs = 0.995;
ycs = 0.005;

M = F + S;
My = (F*xcf + S*ycs)/M; % y coordinate of the mixture point
Mx = (F*xbf + S*ybs)/M; % x coordinate of the mixture point
if (0 < My) && (My <= 0.2)
slope = 0 + (My - 0) * tie_slope(1) / 0.2;
elseif (0.2 < My) && (My <= 0.4)
slope = tie_slope(2) + (My - 0.2) * (tie_slope(3) - tie_slope(2)) / (0.4 - 0.2);
elseif (0.4 < My) && (My <= 0.6)
slope = tie_slope(3) + (My - 0.4) * (tie_slope(4) - tie_slope(3)) / (0.6 - 0.4);
elseif (My > 0.6)
slope = tie_slope(4) + (My - 0.6) * 0.25 / (My - 0.6);
end

x_vals = linspace(0, 1, 10);
y_vals = slope * (x_vals - Mx) + My;
plot(x_vals, y_vals, 'r-', 'linewidth', 1.5);
plot(Mx, My, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'g');

p1 = polyfit(B(1:10), C(1:10), 5);
f1 = polyval(p1, B(1:10));

p2 = polyfit(B(11:19), C(11:19), 3);
f2 = polyval(p2, B(11:19));

syms x y
% raffinate
[xbr,xcr] = vpasolve([y == poly2sym(p1),y == My + slope*(x - Mx)],[x,y],[0 0.3; 0 0.65]); %xlimit ,ylimit
% Extract 
[ybe,yce] = vpasolve([y == poly2sym(p2), y == My + slope*(x-Mx)],[x,y],[0.3 1; 0 0.65]);

plot([double(xbr) double(ybe)],[double(xcr) double(yce)],'o-','Color',[0,0.25,1],'linewidth',0.35);
plot(double(Mx),double(My),'bo');
text(double(xbr-0.05),double(xcr),['R1']);
text(double(ybe),double(yce),[' - ','E1']);
text(double(Mx),double(My+0.02),['M1']);

syms R E
[r,e] = solve([R + E - M,((R*xcr) + (E*yce) - (M*My)),R>0,E>0],[R,E]);  
e=double(e);
r=double(r);
yce=double(yce);
mass_c_ex=e*yce;
