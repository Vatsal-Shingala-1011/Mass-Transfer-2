%% Question 2:

B = [0 0.004 0.006 0.01 0.02 0.03 0.036 0.07 0.13 0.16 0.19 0.23 0.26 0.5 0.63 0.71 0.78 0.84 0.9 0.95 1];
C = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.4025 0.405 0.402 0.4 0.35 0.3 0.25 0.2 0.15 0.1 0.05 0];
A = 1 - B - C;

%==========================================================================
%------ extra data points are added to data to make the curve smoother
%==========================================================================

tiexc = [0.04 0.083 0.13 0.215 0.395];
tieyc = [0.035 0.068 0.09 0.145 0.31];

% figure 1
figure(1)
plot(B,C,'bo-');grid on;hold on;
plot([0 1 0 0],[0 0 1 0],'k-.','linewidth',1.25)
xlabel('xB');ylabel('xC,yC');grid on;
title('raw data')

p1 = polyfit(B,C,4);
f1= polyval(p1,B);

figure(2)
plot(B,C,'--go',...
    'linewidth',2,...
    'MarkerSize',6.5,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.6,0.6,0.6]);hold on;grid on;
xlabel('xB');ylabel('xC,yC');hold on; 
plot([0 1 0 0],[0 0 1 0],'k-.','linewidth',1.25)
plot(B,f1,'k');grid on;
xlabel('xB');ylabel('xC,yC');
title('fitted data');ylim([0 1]);

% interpolation
B1 = zeros(1,10001);
C1 = zeros(1,10001);

for i = 1:10001
    B1(i) = (i-1)/10000;
    C1(i) = round(interp1(B,C,(i-1)/10000),5);
end

%calculations
%stages = input('enter the number of stages');

solute_in_extract = zeros(1, length(600:20:840));
index = 1;

for p=600:20:840
stages = 3;
figure(p)
plot(B1,C1,'b','linewidth',1.00);grid on;hold on;
plot([0 1 0 0],[0 0 1 0],'k-.','linewidth',1.25)
xlabel('xB');ylabel('xC,yC');title('interpolated data - overall process schematic')

tiexb = zeros(1,length(tiexc));
tieyb = zeros(1,length(tiexc));

for i = 1:length(tiexc)
    index1 = find(abs(C1-tiexc(i)) < 0.001,1,'first');
    index2 = find(abs(C1-tieyc(i)) < 0.001,1,'last');    
    XX = [B1(index1) B1(index2)] ;
    YY = [tiexc(i) tieyc(i)];
    tiexb(i) = B1(index1);
    tieyb(i) = B1(index2);
    plot(XX',YY','mo:','linewidth',1.25)
end   

F1 = [0 0.98];
S1 = [0.35 0.02];
plot(F1,S1,'g^-','linewidth',0.35)
text(0,0.35,'F')
text(0.98,0.02,'- S1,S2,S3')

S = [p  p p];
F = 1000;
xbf = 0;
xcf = 0.35;
ybs = 0.98;
ycs = 0.02;
    
tie_slope = zeros(1,length(tiexc));
for i =  1:length(tiexc)
    tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i));
end

r = ones(1,stages);
e = ones(1,stages);
xbr = ones(1,stages);xcr = ones(1,stages);
ybe = ones(1,stages);yce = ones(1,stages);
M = ones(1,stages);
Mx = ones(1,stages);
My = ones(1,stages);

solute = 0;

syms R E
for i = 1:stages

    M(i) = F + S(i);
    My(i) = (F*xcf + S(i)*ycs)/M(i);
    Mx(i) = (F*xbf + S(i)*ybs)/M(i);
    
    if  ((0 < My(i)) && (My(i) <= 0.04))
        slope = 0 + (My(i) - 0)*tie_slope(1)/(0.04);
        
    elseif((0.04< My(i)) && (My(i) <= 0.083))
        slope = tie_slope(1) + (My(i) - 0.04)*tie_slope(2)/(0.083 - 0.04);
        
    elseif ((0.083 < My(i)) && (My(i) <= 0.13))
        slope = tie_slope(2) + (My(i) - 0.083)*tie_slope(3)/(0.13 - 0.083);    
        
    elseif ((0.13 < My(i)) && (My(i) <= 0.215))
        slope = tie_slope(3) + (My(i) - 0.13)*tie_slope(4)/(0.215 - 0.13);    
        
    elseif ((0.215 < My(i)) && (My(i) <= 0.395))
        slope = tie_slope(4) + (My(i) - 0.215)*tie_slope(5)/(0.395 - 0.215);    
          
    elseif((My(i) > 0.395))
        slope = tie_slope(5) + (My(i) - 0.395)*(-0.155)/(0.4 - 0.395);
    
    end
    
    syms x y
    [xbr(i),xcr(i)] = vpasolve([y == poly2sym(p1),y == My(i) + slope*(x - Mx(i))],[x,y],[0 0.2; 0 0.405]);
    [ybe(i),yce(i)] = vpasolve([y == poly2sym(p1),y == My(i) + slope*(x - Mx(i))],[x,y],[0.2 1; 0 0.405]);
 figure(p)
    plot([double(xbr(i)) double(ybe(i))],[double(xcr(i)) double(yce(i))],'o-','Color',[0,0.25,1],'linewidth',0.35)
    plot(double(Mx(i)),double(My(i)),'bo')
    
    text(double(xbr(i)-0.05),double(xcr(i)),['R',num2str(i),' - '])   
    text(double(ybe(i)),double(yce(i)),[' - ','E',num2str(i)])   
    text(double(Mx(i)),double(My(i)+0.02),['M',num2str(i)])  
    
    [r(i),e(i)] = solve([R + E - M(i),((R*xcr(i)) + (E*yce(i)) - (M(i)*My(i))),R>0,E>0],[R,E]);

    % Calculate the amount of solute in the extract phase
    sol = e(i) * yce(i);
    solute = solute + sol;

    F = r(i);
    xcf = xcr(i);
    xbf = xbr(i);
    
end

solute_in_extract(index) = solute;
index = index + 1;

end

% Plot the amount of solute in the extract phase as a function of solvent amount
figure(10)
plot(600:20:840, solute_in_extract, 'b-');
xlabel('Solvent amount (kg)');
ylabel('Solute in extract phase');
title('Removal of solute as a function of solvent amount');
grid on;

figure(11)
for i = 1:length(xcr)-1
    index1 = find(abs(C1-xcr(i)) < 0.001,1,'first');
    XX2 = [B1(index1) 0.98] ;
    YY2 = [xcr(i) 0.02];    
    plot(XX2',YY2','*-','Color',[1,0,0],'linewidth',0.35)
end