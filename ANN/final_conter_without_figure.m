%final 
clc;
clear;

B = [0.04 0.05 0.07 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.993 1];
C = [0 0.11 0.26 0.375 0.474 0.487 0.468 0.423 0.356 0.274 0.185 0.09 0.001 0];

tiexc = [0.1 0.245 0.426];
tiexb = [0.048 0.065 0.133];
tieyc = [0.098 0.242 0.409];
tieyb = [0.891 0.736 0.523];

% Slopes of the tie lines 
tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc)
    tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i));
end

% Draw the feed line (FS line)
F = [0, 0.5];
S = [0.995, 0.005];
xF = F(1);
yF = F(2);
xS = S(1);
yS = S(2);
mm=(S(2)-F(2))/(S(1)-F(1));

p1 = polyfit(B,C,6);
xbf = 0;
xcf = 0.5;
xc_in_feed=xcf;

filename = 'conter_final.xlsx';
if exist(filename, 'file') == 0
    headers = {'Feed', '%C_removed', 'Stages'};
    xlswrite(filename, headers, 'Sheet1', 'A1');
end

% Loop over different feed values
F_values = 1600:10:2000;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% manas=400:10:800 asit=800:10:1200 aditya=1200:10:1600 vatsal 1600:10:2000
for F = F_values
    disp(F);
    Feed = F;
    S = F;
    Rnx_values = 0.0401:0.0008:0.07;
    for Rnx = Rnx_values     
        Rny = polyval(p1,Rnx);
        M = F + S;
        ycs = Rny;
        ybs = 1-Rny;
%         My = (F*xcf + S*ycs)/M;
        Mx = (F*xbf + S*ybs)/M;
        My= mm*Mx + .5;

        m = (My - Rny) / (Mx - Rnx);
        b = Rny - m * Rnx;
        
        syms x y;
        [ybe1, yce1] = vpasolve([y == poly2sym(p1), y == m * x + b],[x, y],[.5 1;0 .7]);
        
        m2 = (yce1 - xcf) / (ybe1 - xbf);
        b2 = xcf - m2 * xbf;
        
        m3 = (ycs - Rny) / (ybs - Rnx);
        b3 = Rny - m3 * Rnx;
        
        eq1 = y == m2 * x + b2;
        eq2 = y == m3 * x + b3;
        [Px, Py] = solve([eq1, eq2], [x, y]);
        
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
            
            intercept = (yce(i) - slope * ybe(i)); 
            
            [xbr(i), xcr(i)] = vpasolve([y == poly2sym(p1), y == slope * x + intercept],[x, y],[0 .5;0 .7]);

            
            m4 = (Py - xcr(i)) / (Px - xbr(i));
            b4 = xcr(i) - m4 * xbr(i); 
            
            [ybe(i+1), yce(i+1)] = vpasolve([y == poly2sym(p1), y == m4 * x + b4],[x, y],[.5 1;0 .7]);
            
            if xcr(i) <= Rny || xcr(i) < Rny+0.005
                solute_removed_pr = ((xc_in_feed-Rny) / xc_in_feed) * 100;
                existing_data = xlsread(filename);
                next_row = size(existing_data, 1) + 1;
                data = [Feed, solute_removed_pr, ii];
                xlswrite(filename, data, 'Sheet1', ['A' num2str(next_row)]);
                break;
            end

            ii = ii+1;
        end
    end
end

