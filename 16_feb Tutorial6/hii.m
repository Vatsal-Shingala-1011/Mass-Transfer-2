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

for i = 1:length(tiexc)
plot([tiexb(i) tieyb(i)], [tiexc(i) tieyc(i)],'mo:','linewidth',1.25);
end

tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc)
tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i));
end

F = [0, 0.45];
S = [0.96, 0.04];
plot(S,F , 'g', 'linewidth', 0.35)
text(0,  0.45,'F');
text(0.96,  0.04,'S');

p1 = polyfit(B,C,7);
X_=0.04:0.02:1;
f1 = polyval(p1,X_);
plot(X_,f1)

S = 2500;
F = 2000;
M = F + S;

xcf = 0.45;

ep=0.01;
for i=0:0.01:0.5 
    if polyval(p1, i)<=.04+ep && polyval(f1, i)>=.04-ep
%         disp(i);
        Rnx=i;
        text(i,0.04 ,'Rn');
        break;
    end
end
Rny = 0.04; 
xbf = 0;
ycs = Rny;
ybs = 1-Rny;
My = (F*xcf + S*ycs)/M; 
Mx = (F*xbf + S*ybs)/M; 
text(Mx, My,'M');
plot([Rnx, Mx], [Rny, My], 'o'); 
hold on;
m = (My-Rny)/(Mx-Rnx);
b = Rny - (m * Rnx);
xx = linspace(Rnx, 1, 500); 
yy = m * xx + b;
plot(xx, yy, 'r');

syms x y;
[ybe1,yce1] = vpasolve([y == poly2sym(p1), y == m * x + b],[x,y],[.5 1;0 0.6]);
text(ybe1,  yce1,'E');

m2 = (yce1 - xcf) / (ybe1 - xbf);
m3 = (ycs - Rny) / (ybs - Rnx);
b2 = xcf - m2 * xbf; 
b3 = Rny - m3 * Rnx;
xx2 = linspace(0, 2, 200); 
yy2 = m2 * xx2 + b2;
xx3 = linspace(0, 2, 200); 
yy3 = m3 * xx3 + b3;

plot(xx2, yy2, 'b');
plot(xx3, yy3, 'b'); 

syms x y;
eq1 = y == m2*x + b2;
eq2 = y == m3*x + b3;
[px, py] = solve([eq1, eq2], [x, y]);
text(px,  py,'P(delta)');


stages=5;
xbr = ones(1,stages);
xcr = ones(1,stages);
ybe = ones(1,stages);
yce = ones(1,stages);
yce(1)=yce1;
ybe(1)=ybe1;
syms x y

for i=1:stages

    if (0 < yce(i)) && (yce(i) <= 0.098)
    slope = 0 + (yce(i) - 0) * 0 / 0.098;
    elseif (0.1 < yce(i)) && (yce(i) <= 0.25)
    slope = tie_slope(1) + (yce(i) - 1) * (tie_slope(2) - tie_slope(1)) / (0.25 - 0.1);
    elseif (0.25 < yce(i)) && (yce(i) <= 0.41)
    slope = tie_slope(2) + (yce(i) - 0.25) * (tie_slope(3) - tie_slope(2)) / (0.41 - 0.25);
    end

    x_vals = linspace(0, ybe(i), 10);
    y_vals = slope * (x_vals - ybe(i) ) + yce(i);

    c = (yce(i) - slope * ybe(i)); 
    plot(x_vals, y_vals, 'r', 'linewidth', .05);
  
    syms x y;
    [xbr(i),xcr(i)] = vpasolve([y == poly2sym(p1), y == slope * x + c],[x,y],[0 .5;0 .7]); 
    text(xbr(i),xcr(i),['R', num2str(i)]);
    
    m4 = (py - xcr(i)) / (px - xbr(i));
    b4 = xcr(i) - m4 * xbr(i);
    xx4 = linspace(0, px, 200); 
    yy4 = m4 * xx4 + b4;
    plot(xx4, yy4, 'b'); 

    [ybe(i+1),yce(i+1)] = vpasolve([y == poly2sym(p1), y == m4 * x + b4],[x,y],[.5 1;0 .7]); 
    text(ybe(i+1),  yce(i+1),['E', num2str(i+1)]);
    
end
