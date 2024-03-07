clc;
clear;
wA = [0.002 0.001 0.0 0.0 0.0 0.0];
wB = [0.952 0.967 0.979 0.989 0.994 0.998];
wC = [0.046 0.032 0.021 0.011 0.006 0.002];
wA_dash = [0.432 0.417 0.401 0.398 0.397 0.396];
wB_dash = [0.524 0.564 0.586 0.5954 0.5994 0.6028];
wC_dash = [0.026 0.019 0.013 0.0066 0.0036 0.0012];

for i = 1:length(wA)
    xc(i) = (wC_dash(i))/(wC_dash(i)+wB_dash(i)); %under
    Z(i) = (wA(i))/(wB(i)+wC(i)); %over
    yc(i) = (wC(i))/(wC(i)+wB(i));%over
    z(i) = (wA_dash(i))/(wB_dash(i)+wC_dash(i));%
end     

F = 400;
F_dash = 400*(1-0.48);
S = 300;
S_dash = S;
zf = 48.1/(49+2.9);
xcf = 2.9/(2.9+49);
zs = 0;
ycs = 0;

figure(1)
plot(xc, z, 'bo-'); grid on; hold on;
plot(yc, Z, 'bo-'); grid on; hold on;
xlabel('xC,yC'); ylabel('Z');
title('Ponchon-Savarit System');

%plotting FS line
F1 = [xcf 0];
S1 = [zf 0];
plot(F1, S1, 'go-');grid on; hold on;
text(xcf,zf, 'F');
text(0, 0, 'S');

%DETERMINING MIXTURE POINT
M = F_dash+S_dash;
zM = (((F_dash*zf)+(S_dash*zs))/M);
xM = (((F_dash*xcf)+(S_dash*ycs))/M); 
text(xM,zM, 'M');
a4 = plot([xcf ycs],[zf, zs],'-or');
text(xcf, zf, 'Fdash');
text(ycs, zs, 'Sdash');
Wb = 0.6017;
Wc = 0.002;
Wa = 1 - (Wb+Wc);
xcn = Wc/(Wc+Wb);

p1 = polyfit(xc,z,2);
X_=0:0.001:0.06;
f1 = polyval(p1,X_);
plot(X_,f1,'black')

zn = polyval(p1, xcn);
text(xcn,zn, 'Ln');


p2 = polyfit(yc,Z,2);
X_=0:0.001:0.06;
f1 = polyval(p2,X_);
plot(X_,f1,'black')

%finding line form Ln to M 
slop_Ln_M = (zn-zM)/(xcn-xM);
x_values = linspace(0,0.055, 100);
y = zM + slop_Ln_M*(x_values-xM);
plot(x_values, y, 'r-', 'LineWidth', 1.5);  % Plot Ln M line

% finding v1 point by skoving p2 and Ln to M line
syms x y;
[v1,yc1] = vpasolve([y == poly2sym(p2), y == zM + slop_Ln_M*(x-xM)],[x,y],[0 1;0 1]); %range=[x1 x2;y1 y2]
text(v1,  yc1,'V1');

%fingding delta point
%Finding FE and RnS lines
%for FV1 line F(xcf,zf) V1(v1,yc1)
m2 = (zf-yc1)/(xcf-v1);
b2 = yc1 - m2 * v1; %intercept
x_val2 = linspace(-0.14, 0.06, 100); 
y_val2 = m2 * x_val2 + b2;
plot(x_val2, y_val2, 'b');
%for LnS line Ln(xcn,zn) S(0,0)
m3 = zn/xcn;
b3 = 0 - m3 * 0;
x_val3 = linspace(-0.14, 0.06, 100); 
y_val3 = m3 * x_val3 + b3;
plot(x_val3, y_val3, 'b'); 

%Finding delta point by solving FE and RnS lines
syms x y;
eq1 = y == m2*x + b2;
eq2 = y == m3*x + b3;
[Px, Py] = solve([eq1, eq2], [x, y]);
text(Px,  Py,'delta');

%finding L1(l1,xc1)
syms x y;
l1=v1;
l1=double(l1);
xc1=polyval(p1,l1);
text(l1,  xc1,'L1');


m2 = (zf-yc1)/(xcf-v1);
b2 = yc1 - m2 * v1; %intercept
x_val2 = linspace(-0.14, 0.06, 100); 
y_val2 = m2 * x_val2 + b2;
plot(x_val2, y_val2, 'b');
%finding the L1 delta line
%for L1 delta line L1(l1,xc1) delta(Px,Py)
m2 = (xc1-Py)/(l1-Px);
b2 = Py - m2 * Px; %intercept
x_val2 = linspace(-0.14, 0.06, 100); 
y_val2 = m2 * x_val2 + b2;
plot(x_val2, y_val2, 'b');

%finding the V2 point by soving L1 delta line and P2 curve






