%LLE data
 B = [0 0.004 0.006 0.01 0.02 0.03 0.036 0.07 0.13 0.16 0.19 0.23 0.26 0.5 0.63 0.71 0.78 0.84 0.9 0.95 ];
 C = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.4025 0.405 0.402 0.4 0.35 0.3 0.25 0.2 0.15 0.1 0.05 ];
 %FS line
 F1 = [0 0.995];
 S1 = [0.5 0.005];
xbf = 0;
xcf = 0.35;
ybs = 0.98;
ycs = 0.02;

 figure(1)
 plot(B,C,'bo-');grid on;
 hold on;

 plot(F1,S1,'g^-','linewidth',0.35)
 text(0,0.5,'F')
 text(0.995,0.005,'- S')

 
plot([0 1 0 0],[0 0 1 0], 'k-.', 'linewidth',1.25);
xlabel('Xb');ylabel('Xc');
title('raw data');

 tie_slope=[-0.0052   -0.0162   -0.0443   -0.0849   -0.1758];
  %tie line intervals
  Mx=0.45;
  My=0.19;
 if  ((0 < My) && (My <= 0.04));
        slope = 0 + (My - 0)*tie_slope(1)/(0.04);
 elseif((0.04< My) && (My <= 0.083));
        slope = tie_slope(1) + (My - 0.04)*tie_slope(2)/(0.083 - 0.04);
 elseif ((0.083 < My) && (My <= 0.13));
        slope = tie_slope(2) + (My - 0.083) * (tie_slope(3) - tie_slope(2)) / (0.13 - 0.083);
 elseif ((0.13 < My) && (My <= 0.215));
        slope = tie_slope(3) + (My - 0.13) * (tie_slope(4) - tie_slope(3)) / (0.215 - 0.13);
 elseif ((0.215 < My) && (My <= 0.395));
        slope = tie_slope(4) + (My - 0.215) * (tie_slope(4) - tie_slope(3)) / (0.395 - 0.215);
 elseif((My > 0.395));
        slope = tie_slope(5) + (My - 0.395)*(-0.155)/(0.4 - 0.395);
 end

    % raffinate
    [xbr,xcr] = vpasolve([y == poly2sym(p1),y == My + slope*(x - Mx)],[x,y],[0 0.1508; 0 0.6107]);    
    % Extract 
    [ybe,yce] = vpasolve([y == poly2sym(p2), y == My + slope*(x - Mx)],[x,y],[0.16 1; 0.6093 0]);
    
    plot([double(xbr) double(ybe)],[double(xcr) double(yce)],'o-','Color',[0,0.25,1],'linewidth',0.35);
    plot(double(Mx),double(My),'bo');
    [r,e] = solve([R + E - M,((R*xcr) + (E*yce) - (M*My)),R>0,E>0],[R,E]);
    F = r; %This is the feed in the subsequent stage
    xcf = xcr; % Respective concentration
    %display(r);
    %display(e);

x_vals = linspace(0, 1, 10);
y_vals = slope * (x_vals - Mx) + My;
plot(x_vals, y_vals, 'r-', 'linewidth', 1.5);
plot(Mx, My, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'g');



p1 = polyfit(B(1:10), C(1:10), 3);
f1 = polyval(p1, B(1:10));

p2 = polyfit(B(11:19), C(11:19), 2);
f2 = polyval(p2, B(11:19));

figure(2)
plot(B(1:10), f1, 'g');
hold on;
plot(B(11:19), f2, 'r');
hold on;
scatter(B,C,'o');
legend("raffinate", "extract", "experimental data");
title("curve fitting");
xlabel("xB, yB");
ylabel("xC, yC");

figure(3)
xcr = [0.003:0.001:0.13];
f1 = polyval(p1, xcr);
plot(xcr,f1 , 'r'); hold on;
scatter(B(1:10),C(1:10),'b');



%part e
yce3 = zeros(1,80);
k = 1;
for j = 600:10:840
    stages = 3;
    S = [j j j];
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

        [xbr(i),xcr(i)] = vpasolve([y == poly2sym(p1),y == My(i) + slope*(x - Mx(i))],[x,y],[0 0.1508; 0 0.6107]); % for the raffinate part

        [ybe(i),yce(i)] = vpasolve([y == poly2sym(p2), y == My(i) + slope*(x-Mx(i))],[x,y],[0.16 1; 0.6093 0]); % for the extract part

        [r(i),e(i)] = solve([R + E - M(i),((R*xcr(i)) + (E*yce(i)) - (M(i)*My(i))),R>0,E>0],[R,E]);
        F = r(i); % This is the feed in the subsequent stage
        xcf = xcr(i); % Respective concentration
        yce3(k) = yce(3);
    end
    k = k+1;
end

figure(4)
x = [400:10:1200];
plot(x,yce3, '-bo');
xlabel("amount of solvent");
ylabel("Solute level");
title("Solvent level vs no of stages");