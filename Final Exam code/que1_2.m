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

xb = [0.33 0.302 0.272 0.242 0.213 0.1817 0.1492 0.1148];
 xc = [0 0.0336 0.0682 0.1039 0.1419 0.1817 0.224 0.268];
 xa = 1-xb-xc;
 yc = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7];
 for i = 1: length(yc);
    yb(i) = -yc(i) + 1;
 end
 ya = 1-yb-yc;
 
%SETTING UP THE UNDERFLOW AND OVERFLOW
for i = 1:length(wA)
Xc(i) = (wC_dash(i))/(wC_dash(i)+wB_dash(i));
Z(i) = (wA(i))/(wB(i)+wC(i));
Yc(i) = (wC(i))/(wC(i)+wB(i));
z(i) = (wA_dash(i))/(wB_dash(i)+wC_dash(i));
end

F = 2000;
F_dash = 2000*(0.26);
S = 2100;
S_dash = S;
zf = (100-26)/26; %solid/non solid
xcf = 1;
zs = 0;%no solid in overflow
ycs = 0.005;


figure(1)
plot(Xc, z, 'bo-'); grid on; hold on;
plot(Yc, Z, 'bo-'); grid on; hold on;
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
disp("xcn is ")
disp(xcn)
disp("above line is underflow and below line in graph is overflow")
