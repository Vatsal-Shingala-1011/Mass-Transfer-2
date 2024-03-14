clc;
clear;

% wC = [0.046 0.032 0.021 0.011 0.006 0.002];
wC = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7];
wB_dash = [0.33 0.302 0.272 0.242 0.213 0.1817 0.1492 0.1148];
wC_dash = [0 0.0336 0.0682 0.1039 0.1419 0.1817 0.224 0.268];
for i =1:length(wC)
    wA_dash(i) = 1-wB_dash(i)-wC_dash(i);
    wA(i) = 0;
    wB(i)= 1-wC(i);
end

%SETTING UP THE UNDERFLOW AND OVERFLOW
for i = 1:length(wA)
xc(i) = (wC_dash(i))/(wC_dash(i)+wB_dash(i));
Z(i) = (wA(i))/(wB(i)+wC(i));
yc(i) = (wC(i))/(wC(i)+wB(i));
z(i) = (wA_dash(i))/(wB_dash(i)+wC_dash(i));
end

F = 2000;
F_dash = 2000*(0.26);
S = 800;
S_dash = S;
zf = (100-26)/26; %solid/non solid
xcf = 1;
zs = 0;%no solid in overflow
ycs = 0.005;


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
Wb = 0;
Wc = 0.015;
Wa = 1 - (Wb+Wc);
xcn = 0.039;

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
x_values = linspace(0,0.28, 100);
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
x_val2 = linspace(-0.2, 1.01, 100); 
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
i=1;
l1=xcn+1; %to run while once
while (xcn<l1)
    %finding L1(l1,xc1)
    l1=v1;
    l1=double(l1);
    xc1=polyval(p1,l1);
    text(l1,  xc1,['L' num2str(i)]);
    
    %finding the L1 delta line
    %for L1 delta line L1(l1,xc1) delta(Px,Py)
    m2 = (xc1-Py)/(l1-Px);
    b2 = Py - m2 * Px; %intercept
    x_val2 = linspace(-0.14, 0.06, 100); 
    y_val2 = m2 * x_val2 + b2;
    plot([l1,Px],[xc1,Py])
    
    %finding the V2 point by soving L1 delta line and P2 curve
    syms x y;
    [v2,yc2] = vpasolve([y == poly2sym(p2), y == m2 * x + b2],[x,y],[-0.1 1;-0.1 1]); %range=[x1 x2;y1 y2]
    text(v2,  yc2,['V' num2str(i+1)]);
      
    v1=v2;
    yc1=yc2;
    i=i+1;
end



