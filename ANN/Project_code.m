clc;
clear;

B = [0 0.0011 0.0021 0.0031 0.0058 0.0136 0.0372 0.11 0.1259 0.1376 0.1508 0.16 0.3357 0.4751 0.608 0.7698 0.8872 0.9982 1];
C = [0 0 0.1 0.2 0.3 0.4 0.5 0.5821 0.6 0.6058 0.6107 0.6093 0.5819 0.4944 0.3748 0.2223 0.1079 0 0];

tiexc = [0.0 0.2 0.4 0.6];
tieyc = [0.0 0.2223 0.4944 0.6107];
tiexb = [0.0011 0.0031 0.0136 0.1259];
tieyb = [0.9982 0.7698 0.4751 0.1508];

tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc)
    tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i));
end

% Define ranges for F, ycs, and xcf
F_range = 100:30:2000;
ycs_range = 0:0.02:0.1;
ybs = 1-ycs_range;
xcf_range = 0.3:0.05:0.8;

% Initialize an empty cell array to store results
results = cell(1, 5);
results{1} = {'F', 'ycs', 'xcf', 'mass_c_ex', 'yce'};

row = 2; % Starting row for storing data in Excel

for F = F_range
    for ycs = ycs_range
        for xcf = xcf_range
            S = F;
            M = F + S;
            My = (F*xcf + S*ycs)/M; % y coordinate of the mixture point
            Mx = (F*0 + S*1)/M; % x coordinate of the mixture point
            
            % Calculate slope based on My
            if (0 < My) && (My <= 0.2)
                slope = 0 + (My - 0) * tie_slope(1) / 0.2;
            elseif (0.2 < My) && (My <= 0.4)
                slope = tie_slope(2) + (My - 0.2) * (tie_slope(3) - tie_slope(2)) / (0.4 - 0.2);
            elseif (0.4 < My) && (My <= 0.6)
                slope = tie_slope(3) + (My - 0.4) * (tie_slope(4) - tie_slope(3)) / (0.6 - 0.4);
            elseif (My > 0.6)
                slope = tie_slope(4) + (My - 0.6) * 0.25 / (My - 0.6);
            end
            slope=slope+0.0001;
            % Calculate y_vals using slope
            x_vals = linspace(0, 1, 10);
            y_vals = slope * (x_vals - Mx) + My;

            % Polynomial fitting for Raffinate and Extract regions
            p1 = polyfit(B(1:10), C(1:10), 5);
            p2 = polyfit(B(11:19), C(11:19), 3);
            
            % Solve equations for Raffinate and Extract
            syms x y
            [xbr, xcr] = vpasolve([y == poly2sym(p1), y == My + slope * (x - Mx)], [x, y], [0, 0.1508; 0, 0.6107]);
            [ybe, yce] = vpasolve([y == poly2sym(p2), y == My + slope * (x - Mx)], [x, y], [0.16, 1; 0.6093, 0]);
                       
            % Solve for R and E
            syms R E
            [r, e] = solve([R + E - M, ((R * xcr) + (E * yce) - (M * My)), R > 0, E > 0], [R, E]);
            e = double(e);
            r = double(r);
            yce = double(yce);
            mass_c_ex = e * yce;

            % Store the results in the cell array
            results{row, 1} = F;
            results{row, 2} = ycs;
            results{row, 3} = xcf;
            results{row, 4} = mass_c_ex;
            results{row, 5} = yce;
            row=row+1;
            % Output the results
            fprintf('F = %d, ycs = %.2f, xcf = %.2f: mass_c_ex = %.4f, yce = %.4f, e=%d\n', F, ycs, xcf, mass_c_ex, yce,e);
        end
    end
end

% Write the cell array to an Excel file
xlswrite('hi.xlsx', results);
