clc;
clear;
xb = [0.33 0.302 0.272 0.242 0.213 0.1817 0.1492 0.1148];
 xc = [0 0.0336 0.0682 0.1039 0.1419 0.1817 0.224 0.268];
 xa = 1-xb-xc;
 yc = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7];
 for i = 1: length(yc);
    yb(i) = -yc(i) + 1;
 end
 ya = 1-yb-yc;

 %plotting in right triangular system
 figure(1)
 plot([0 1 0 0], [0 0 1 0], 'k-', 'linewidth', 1.25); grid on; hold on;
 plot(xc, xb, 'bo-'); grid on; hold on;

 plot(yc, yb, 'bo-'); grid on; hold on;
 xlabel('xC,yC'); ylabel('xB,yB');
 title('Right-Angled Triangular System');
 
 %plotting FS line
 F1 = [0.3 0.02];
 S1 = [0.02 1];
 plot(F1, S1, 'go-');grid on; hold on;
 text(0.3, 0.02, 'F');
 text(0.02, 1, 'S');

 tie_xb = xb;
 tie_xc = xc;
 tie_yb = yb;
 tie_yc = yc;
%plotting the tie lines
for i = 1:length(tie_xb)
   plot([tie_xc(i) tie_yc(i)],[tie_xb(i) tie_yb(i)], 'ro:');grid on; hold on;
end   

tie_slope = ones(1, length(tie_xc));
tie_slope(1)=60;
for i = 2:length(tie_xc)
  tie_slope(i) = (tie_xb(i)-tie_yb(i))/(tie_xc(i)- tie_yc(i));
end

%fitting to polynomials
 p1 = polyfit(xc, xb, 1);
 f1 = polyval(p1, xc);
 p2 = polyfit(yc, yb, 1);
 f2 = polyval(p2, yc);

 %given data
 S =  2100;
 F =  2000 ;
 xcf = 0.26;
 xbf = 0;
 ycs = 0.005;
 ybs = 1-ycs;
 xcr_n = 0.015;

 %cross current process
 i=1;
 xcr=1;%any values so that loop can run for first step
 Sum_solvent = 0;
 solute_leaving = 0;

 while (xcr > xcr_n)
    %plotting the M(i) points
    M = F+S;
    Mx = (F*xcf+S*ycs)/M;
    My = (F*xbf+S*ybs)/M;
    plot(Mx, My, 'ko')
    text(Mx, My, ['M' num2str(i) ''])

    % tie line interpolation
    for k =2:length(tie_xc)
        if Mx<=((tie_yc(k)+tie_xc(k))/2) && Mx>=((tie_yc(k-1)+tie_xc(k-1))/2)
            break;
        end
    end

    if (k==2)
        slope = tie_slope(k-1) + ((Mx-tie_yc(k-1))/(tie_yc(k)-tie_yc(k-1))) * (tie_slope(k) - tie_slope(k-1));
    else
        slope = tie_slope(k-1) + (Mx-tie_yc(k-1))/(tie_yc(k)-tie_yc(k-1));
    end

    % calculating the concentrations after each stage
    syms x y
    x_vals = linspace(0, 1, 10);
    y_vals = slope * (x_vals - Mx ) + My;
    intercept = slope * (0 - Mx ) + My; %y at x_vals=0
    [xcr, xbr] = vpasolve([y ==poly2sym(p1,x), y == slope * x + intercept], [x,y], [0 1;0 1]);
    [yce, ybe] = vpasolve([y == poly2sym(p2,x), y == slope * x + intercept], [x,y], [0 1;0 1]);
 
    plot([xcr yce], [xbr ybe],'bo-');grid on; hold on;
    plot([xcr ycs], [xbr ybs], 'c-');grid on; hold on;
    %calculating the amount of underflow and overflow after each stage
 
    [r, e] = vpasolve([y == M-x, y ==x*((Mx-xcr)/(yce-Mx))], [x,y], [0, Inf; 0, Inf]);
    %printing the results
     fprintf('Amount of underflow after stage_%d =', i);
     disp(r);
     fprintf('Amount of overflow after stage_%d =', i);
     disp(e);
     fprintf('Underflow concentration after stage_%d = ', i);
     disp(xcr);
     fprintf('Overflow concentration after stage_%d = ', i);
     disp(yce);
    %calculating fractional recovery after each stage
    fprintf('Fraction of solute removed in stage_%d =', i);
    disp((yce*e-S*ycs)/(F*xcf)); 

    Sum_solvent = Sum_solvent+S;
    solute_leaving = solute_leaving + (e*yce);

    xbf = xbr;
    xcf = xcr;
    F= r;
    text( xcf, xbf,['L' num2str(i) '']);
    text( yce, ybe,['V' num2str(i) '']);
    i=i+1;
 end
