clc;
clear;

B = [0.04 0.05 0.07 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.993 1];
C = [0 0.11 0.26 0.375 0.474 0.487 0.468 0.423 0.356 0.274 0.185 0.09 0.001 0];
A = 1 - B - C;

figure;
plot(B,C,'bo-');
grid on;
hold on;

plot([0 1 0 0],[0 0 1 0], 'k-.', 'linewidth',1.25);
xlabel('Xb');
ylabel('Xc');
title('Raw data');

tiexc = [0.1 0.245 0.426];
tiexb = [0.048 0.065 0.133];
tieyc = [0.098 0.242 0.409];
tieyb = [0.891 0.736 0.523];

% Plotting the tie lines
for i = 1:length(tiexc)
    plot([tiexb(i) tieyb(i)], [tiexc(i) tieyc(i)],'mo:','linewidth',1.25);
end

% Slopes of the tie lines 
tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc)
    tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i));
end

% Draw the feed line (FS line)
F = [0, 0.5];
S = [0.995, 0.005];
plot(S,F,'g','linewidth',0.35)
text(0, 0.5,'F');
text(0.995, 0.005,'S');

p1 = polyfit(B,C,6);
X_=0.04001:0.02:1;
f1 = polyval(p1,X_);
plot(X_,f1,'black')

xbf = 0;
xcf = 0.5;
xc_in_feed=xcf;

filename = 'CCOOuu.xlsx';
if exist(filename, 'file') == 0
    headers = {'Feed', '%C_removed', 'Stages'};
    xlswrite(filename, headers, 'Sheet1', 'A1');
end

% Loop over different feed values
F_values = 800;
for F = F_values
    disp(F);
    Feed = F;
    S = F;
%     Rnx_values = 0.0401:0.0008:0.066;
Rnx_values=0.0601;
    for Rnx = Rnx_values     
        Rny = polyval(p1,Rnx);
        M = F + S;
        ycs = Rny;
        ybs = 1-Rny;
        My = (F*xcf + S*ycs)/M;
        Mx = (F*xbf + S*ybs)/M;
        text(Mx, My,'M');
        
        plot([Rnx, Mx], [Rny, My], 'o'); 
        hold on;
        m = (My - Rny) / (Mx - Rnx);
        b = Rny - m * Rnx;
        x_val = linspace(Rnx, 1, 100); 
        y_val = m * x_val + b;
        plot(x_val, y_val, 'r');
        
        syms x y;
        [ybe1, yce1] = vpasolve([y == poly2sym(p1), y == m * x + b],[x, y],[.5 1;0 .7]);
        text(ybe1, yce1,'E');
        
        m2 = (yce1 - xcf) / (ybe1 - xbf);
        b2 = xcf - m2 * xbf;
        x_val2 = linspace(0, 1.8, 200); 
        y_val2 = m2 * x_val2 + b2;
        plot(x_val2, y_val2, 'b');
        
        m3 = (ycs - Rny) / (ybs - Rnx);
        b3 = Rny - m3 * Rnx;
        x_val3 = linspace(0, 1.8, 200); 
        y_val3 = m3 * x_val3 + b3;
        plot(x_val3, y_val3, 'b');
        
        eq1 = y == m2 * x + b2;
        eq2 = y == m3 * x + b3;
        [Px, Py] = solve([eq1, eq2], [x, y]);
        text(Px, Py,'P');
        
        stages = 100;
        xbr = ones(1,stages);
        xcr = ones(1,stages);
        ybe = ones(1,stages);
        yce = ones(1,stages);
        yce(1) = yce1;
        ybe(1) = ybe1;
        ii = 1;
        for i = 1:stages
            disp(i);
            if (0 < yce(i)) && (yce(i) <= 0.098)
                slope = 0 + (yce(i) - 0) * 0 / 0.098; 
            elseif (0.098 < yce(i)) && (yce(i) <= 0.249)
                slope = tie_slope(1) + (yce(i) - 0.098) * (tie_slope(2) - tie_slope(1)) / (0.249 - 0.098);
            elseif (0.249 < yce(i)) && (yce(i) <= 0.409)
                slope = tie_slope(2) + (yce(i) - 0.249) * (tie_slope(3) - tie_slope(2)) / (0.409 - 0.249);
            end
            
            x_vals = linspace(0, ybe(i), 10);
            y_vals = slope * (x_vals - ybe(i)) + yce(i);
            intercept = (yce(i) - slope * ybe(i)); 
            plot(x_vals, y_vals, 'r', 'linewidth', .05);
            
            [xbr(i), xcr(i)] = vpasolve([y == poly2sym(p1), y == slope * x + intercept],[x, y],[0 .5;0 .7]);
            text(xbr(i), xcr(i),['R', num2str(i)]);
            
            m4 = (Py - xcr(i)) / (Px - xbr(i));
            b4 = xcr(i) - m4 * xbr(i);
            x_val4 = linspace(0, Px, 200); 
            y_val4 = m4 * x_val4 + b4;
            plot(x_val4, y_val4, 'b'); 
            
            [ybe(i+1), yce(i+1)] = vpasolve([y == poly2sym(p1), y == m4 * x + b4],[x, y],[.5 1;0 .7]);
            text(ybe(i+1), yce(i+1),['E', num2str(i+1)]);
            
            if xcr(i) <= Rny || xcr(i) < Rny+0.005
                solute_removed_pr = ((xc_in_feed-Rny) / xc_in_feed) * 100;
                existing_data = xlsread(filename);
                next_row = size(existing_data, 1) + 1;
                data = [Feed, solute_removed_pr, ii];
                xlswrite(filename, data, 'Sheet1', ['A' num2str(next_row)]);
                break;
            end
            if i==100
                solute_removed_pr = ((xc_in_feed-Rny) / xc_in_feed) * 100;
                existing_data = xlsread(filename);
                next_row = size(existing_data, 1) + 1;
                data = [Feed, solute_removed_pr, 100];
                xlswrite(filename, data, 'Sheet1', ['A' num2str(next_row)]);
                break;
            end
            ii = ii+1;
        end
    end
end

% Plotting solute fraction vs number of stages
figure(2)
plot(xcr, '-bo');
xlabel("Number of stages");
ylabel("Solute Fraction in raffinate");
title("Solute Fraction vs Number of Stages");
