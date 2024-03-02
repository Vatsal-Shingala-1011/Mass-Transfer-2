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
    z(i) = (wA_dash(i))/(wB_dash(i)+wC_dash(i));
    Z(i) = (wA(i))/(wB(i)+wC(i));
    yc(i) = (wC(i))/(wC(i)+wB(i));
 end     

%plotting in right triangular system
 figure(1)
 plot(xc, z, 'bo-'); grid on; hold on;

 plot(yc, Z, 'bo-'); grid on; hold on;
 xlabel('xC,yC'); ylabel('z,Z');
 title('Ponchon-savarit ');

%plotting the tie lines
for i = 1:length(xc)
   plot([xc(i) yc(i)],[z(i) Z(i)], 'ro:');grid on; hold on;
end   

F=400;
S=500;
F_desh=400*(1-.48); %solid free
S_desh=S;  %solid free
M1_desh=F_desh+S_desh;
Zf_desh=48.1/(49+2.9);
Xcf=2.9/(2.9+49);
Zs=0;

%plotting FS line
 F1 = [Xcf,Zf_desh];
 S1 = [0 0];
%  plot(F1, S1, 'go-');grid on; hold on;
 plot([0 Xcf], [0 Zf_desh],'bo-');grid on; hold on;
 text(Xcf,Zf_desh ,'F');
 text(0, 0,'S');

 %DETERMINING THE UNDEFLOW AND OVERFLOW POINTS
    p1 = polyfit(xc, z, 1);
    Ly = polyval(p1,xc);
    p2 = polyfit(yc, Z, 1);
    Vy = polyval(p2,z);

    M1_desh=F_desh+S_desh;
    M1x=(F_desh*Zf_desh+S_desh*Zs)/M1_desh; %x cor of M1_desh
    slop=Zf_desh/Xcf;
    M1y=slop*M1x;
    text(M1x,M1y ,'M1_desh');
    
%     a6 = plot([…..],[…..], '--ok');
%     text(xM, Ly, ['L' num2str(i) '']); 
%     text(xM, Vy, ['V' num2str(i) '']);