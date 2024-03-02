clc;
clear;
wA = [0.002 0.001 0.0 0.0 0.0 0.0];
wB = [0.952 0.967 0.979 0.989 0.994 0.998];
wC = [0.046 0.032 0.021 0.011 0.006 0.002];
wA_dash = [0.432 0.417 0.401 0.398 0.397 0.396];
wB_dash = [0.524 0.564 0.586 0.5954 0.5994 0.6028];
wC_dash = [0.026 0.019 0.013 0.0066 0.0036 0.0012];

%SETTING UP THE UNDERFLOW AND OVERFLOW
for i = 1:length(wA)
xc(i) = (wC_dash(i))/(wC_dash(i)+wB_dash(i));
Z(i) = (wA(i))/(wB(i)+wC(i));
yc(i) = (wC(i))/(wC(i)+wB(i));
z(i) = (wA_dash(i))/(wB_dash(i)+wC_dash(i));
end

F = 400;
Fdash = 400*(1-0.48);
S = 500;
Sdash = S;
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

tie_z = z;
tie_xc = xc;
tie_Z = Z;
tie_yc = yc;

%plotting the tie lines
for i = 1:length(tie_z)
plot([tie_xc(i) tie_yc(i)],[tie_z(i) tie_Z(i)], 'ro:');grid on; hold on;
end

%fitting to polynomials
p1 = polyfit(xc, z, 1);
f1 = polyval(p1, xc);

p2 = polyfit(yc, Z, 1);
f2 = polyval(p2, yc);

plot(xc, f1); hold on;
plot(yc, f2); hold on;

stages = 4;

for i=1:stages
   %plotting the M(i) points
    Mdash = Fdash + Sdash;
    Mx = (Fdash*xcf + Sdash*ycs)/Mdash;
    My = (Fdash*zf + Sdash*zs)/Mdash;
    plot(Mx, My, 'ko')
    text(Mx, My, ['M' num2str(i) ''])
    
    xcr = Mx;
    yce = Mx;
    zr = polyval(p1, xcr);
    ze = polyval(p2, yce);
    plot([xcr yce], [zr ze],'bo-');grid on; hold on;
    plot([xcr ycs], [zr zs], 'c-');grid on; hold on;
    
    %calculating the amount of underflow and overflow after each stage
    Ldash = (Mdash*(My-ze))/(zr-ze);
    Vdash = Mdash-Ldash;
    
    zf = zr;
    xcf = xcr;
    Fdash= Ldash;
    text( xcf, zf,['L' num2str(i) '']);

    text( yce, ze,['V' num2str(i) '']);
end