%1.1. Construct the LLE curve using right angled triangular co-ordinate (GREEN) and
%show the given tie lines using a different color (RED) [Tutorial one]
%% Data
clc;clear;
B = [0, 0.0011, 0.0021, 0.0031, 0.0058, 0.0136, 0.0372, 0.11, 0.1259, 0.1376, 0.1508, 0.16, 0.3357, 0.4751, 0.608, 0.7698, 0.8872, 0.9982, 1];
C = [0, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5821, 0.6, 0.6058, 0.6107, 0.6093, 0.5819, 0.4944, 0.3748, 0.2223, 0.1079, 0, 0];

%2. Show the feed composition and solvent composition as two points on your graphical
%system (BLACK). Draw the FS line. Suppose two tie lines are given with the
%following co-ordinates.
figure(1)
plot(B,C,'bo-');grid on;
hold on;

plot([0 1 0 0],[0 0 1 0], 'k-.', 'linewidth',1.25);
xlabel('Xb');ylabel('Xc');
title('raw data');

tiexc = [0, 0.2, 0.4, 0.6];
tieyc = [0.0, 0.2223, 0.4944, 0.6107];
tiexb = [0.0011, 0.0031, 0.0136, 0.1259];
tieyb = [0.9982, 0.7698, 0.4751, 0.1508];

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
F = [0, 0.5];
S = [0.995, 0.005];
plot(S,F , 'g', 'linewidth', 0.35)
text(0,  0.5,'F');
text(1,  0,'S');

%3. Calculate the amount of extract and raffinate phase after 3
%stage-extraction.

%Step 1: Fit a polynomial for the raffinate and extract phases.
p1 = polyfit(B(1:10),C(1:10),5);
p2 = polyfit(B(10:19),C(10:19),3);

stages = 3;
S = [800 800 800];
F = 1000;

xbf = 0;
xcf = 0.5;
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

for i = 1:3
    M(i) = F + S(i);
    My(i) = (F*xcf + S(i)*ycs)/M(i); % y coordinate of the mixture point
    Mx(i) = (F*xbf + S(i)*ybs)/M(i); % x coordinate of the mixture point
    
    if (0.1 < My(i)) && (My(i) <= 0.3)
    slope = tie_slope(1) + ( My(i)-0.1) * (tie_slope(2) - tie_slope(1)) / (0.3 - 0.1 );
    elseif (0.3 < My(i)) && (My(i) <= 0.5)
    slope = tie_slope(2) + (My(i) - 0.3) * (tie_slope(3) - tie_slope(2)) / (0.5 - 0.3);
    elseif (0.5 < My(i)) && (My(i)<=0.6)
    slope = tie_slope(3) + (My(i) - 0.3) * (tie_slope(4) - tie_slope(3)) / (0.5 - 0.3);
    elseif (My(i)>0.6)
    slope = tie_slope(4) + (My(i) - 0.3) * 0.12 / (0.6107 - 0.5);
    end
    
    % raffinate
    [xbr(i),xcr(i)] = vpasolve([y == poly2sym(p1),y == My(i) + slope*(x - Mx(i))],[x,y],[0 0.1508; 0 0.6107]); %xlimit ,ylimit
    % Extract 
    [ybe(i),yce(i)] = vpasolve([y == poly2sym(p2), y == My(i) + slope*(x-Mx(i))],[x,y],[0.16 1; 0.6093 0]);
    
    plot([double(xbr(i)) double(ybe(i))],[double(xcr(i)) double(yce(i))],'o-','Color',[0,0.25,1],'linewidth',0.35);
    plot(double(Mx(i)),double(My(i)),'bo');
    text(double(xbr(i)-0.05),double(xcr(i)),['R',num2str(i),' - ']);
    text(double(ybe(i)),double(yce(i)),[' - ','E',num2str(i)]);
    text(double(Mx(i)),double(My(i)+0.02),['M',num2str(i)]);

    [r(i),e(i)] = solve([R + E - M(i),((R*xcr(i)) + (E*yce(i)) - (M(i)*My(i))),R>0,E>0],[R,E]);
    F = r(i); %This is the feed in the subsequent stage
    xcf = xcr(i); % Respective concentration
end

figure(2)
xbr = [0: 0.001: 0.17];
xcr = [0.17:0.001:1];

f1 = polyval(p1, xbr);
f2 = polyval(p2, xcr);

plot(xbr, f1, 'g');
hold on;
plot(xcr, f2, 'r');


figure(3)
plot(yce, '-bo');
xlabel("number of stages");
ylabel("solute fraction");
title("Solute fraction vs no. of stages");

figure(4)
solute_level = yce.*e;
plot(solute_level, '-bo');
xlabel("Number of Stages");
ylabel("Solute level");
title("Solute level vs no of stages");