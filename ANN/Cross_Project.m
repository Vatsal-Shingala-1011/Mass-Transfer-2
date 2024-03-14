%tutorial4_e
clc;
clear;

B = [0 0.0011 0.0021 0.0031 0.0058 0.0136 0.0372 0.11 0.1259 0.1376 0.1508 0.16 0.3357 0.4751 0.608 0.7698 0.8872 0.9982 1];
C = [0 0 0.1 0.2 0.3 0.4 0.5 0.5821 0.6 0.6058 0.6107 0.6093 0.5819 0.4944 0.3748 0.2223 0.1079 0 0];
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

p1 = polyfit(B(1:10),C(1:10),5);
p2 = polyfit(B(10:19),C(10:19),3);
xbr = [0: 0.001: 0.17];
xcr = [0.17:0.001:1];
f1 = polyval(p1, xbr);
f2 = polyval(p2, xcr);

plot(xbr, f1, 'r');
hold on;
plot(xcr, f2, 'r');

first=1;

stages = 20;
filename = 'Crrrrrrrross.xlsx';
for F =400:5:2000
        Feed=F;
        S = repmat(F, 1, stages);
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
        
    for i = 1:stages
        M(i) = F + S(i);
        My(i) = (F*xcf + S(i)*ycs)/M(i); % y coordinate of the mixture point
        Mx(i) = (F*xbf + S(i)*ybs)/M(i); % x coordinate of the mixture point
        text(double(Mx(i)),double(My(i)+0.02),['M',num2str(i)]);
        plot(double(Mx(i)),double(My(i)),'bo');
        
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
        [xbr(i),xcr(i)] = vpasolve([y == poly2sym(p1),y == My(i) + slope*(x - Mx(i))],[x,y],[-0.1 0.1508; -1 0.6107]); %xlimit ,ylimit
        % Extract 
        [ybe(i),yce(i)] = vpasolve([y == poly2sym(p2), y == My(i) + slope*(x-Mx(i))],[x,y],[-0.1 1; 0.6093 -1]);
        if xcr(i)<=0
            break;
        end
        plot([double(xbr(i)) double(ybe(i))],[double(xcr(i)) double(yce(i))],'o-','Color',[0,0.25,1],'linewidth',0.35);
        
        text(double(xbr(i)-0.05),double(xcr(i)),['R',num2str(i),' - ']);
        text(double(ybe(i)),double(yce(i)),[' - ','E',num2str(i)]);
    
        ff = [xbr(i),xcr(i)];
        ss = [0.995, 0.005];
        plot(ss,ff , 'g', 'linewidth', 0.35)
    
        [r(i),e(i)] = solve([R + E - M(i),((R*xcr(i)) + (E*yce(i)) - (M(i)*My(i))),R>0,E>0],[R,E]);
        F = r(i); %This is the feed in the subsequent stage
        
        solute_removed_pr(i) = ( (xc_in_feed-xcr(i)) /xc_in_feed) * 100;
        xcf = xcr(i); % Respective concentration
    end
    number_of_stages=1:i-1;
    Feed = Feed*ones(1,i-1);
    if first==1
        data = [Feed',number_of_stages', solute_removed_pr'];
        xlswrite(filename, data);
        first=2;
    end
    existing_data = xlsread(filename);
    new_data = [Feed',number_of_stages', solute_removed_pr'];
    combined_data = [existing_data; new_data];
    xlswrite(filename, combined_data);
    
end

figure(3)
plot(solute_removed_pr, '-bo');
xlabel("Number of stages");
ylabel("Solute removed percentage");
title("Solute removed percentage vs no of stages");

