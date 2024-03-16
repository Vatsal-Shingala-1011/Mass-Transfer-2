
%Data provided ---------------------------
B = [0.0155 0.017 0.025 0.038 0.06 0.122 0.225 0.354 0.472 0.609 0.735 0.857 0.9533 0.9788];
C = [0 0.0285 0.117 0.205 0.262 0.328 0.346 0.336 0.308 0.246 0.173 0.089 0.0187 0];
A = 1 - B - C;
tiexc = [0 0.0285 0.117 0.205 0.262 0.328 0.346];
tieyc = [0 0.0187 0.089 0.173 0.246 0.308 0.336];
tiexb = [0.0155 0.017 0.025 0.038 0.06 0.122 0.225];
tieyb = [0.9788 0.9533 0.857 0.735 0.609 0.472 0.354];

%%----------------------------------------
figure(1)
plot(B,C,'bo-');grid on;
hold on; 
plot([0 1 0 0], [0 0 1 0])
xlabel('');ylabel('');
title('raw data');

p1 = polyfit(B(1:7), C(1:7), 5);
f1= polyval(p1, B(1:7));
plot(B(1:7),f1);grid on;

p2 = polyfit(B(8:14), C(8:14), 3);
f2= polyval(p2, B(8:14));
plot(B(8:14),f2);grid on;

%%---------------------------------------
B1 = zeros(1,10001);
C1 = zeros(1,10001);

for i = 1:10001
    B1(i) = (i-1)/10000;
    C1(i) = round(interp1(B,C,B1(i)),3);
end
plot(B1,C1);grid on;hold on;
plot([0 1 0 0],[0 0 1 0])
xlabel('xB');ylabel('xC,yC');title('interpolated data - overall process schematic')

%%Calculation for given tie lines----------------------


for i = 1:length(tiexc)   
    plot([tiexb(i) tieyb(i)],[tiexc(i) tieyc(i)]);
end

tie_slope=zeros(1,length(tiexc));
for i=1:length(tiexc)
    tie_slope(i) = (tieyc(i)-tiexc(i))/(tieyb(i)-tiexb(i));
end 

%%Given Values--------------------------
% stages = 3;
S = 1300;
F = 1000;
xbf = 0;
xcf = 0.35;
ybs = 1;
ycs = 0;

plot([xbf ybs],[xcf ycs]) % Line FS

%%To be calculated-----------------
r = ones(1,20);
e = ones(1,20);
xbr = ones(1,20);
xcr = ones(1,20);
ybe = ones(1,20);
yce = ones(1,20);
M = ones(1, 20);
Mx = ones(1,20);
My = ones(1, 20);
i=1;

M(i) = F + S;
Mx(i) = (F*xbf + S*ybs)/M(i);
My(i) = (F*xcf + S*ycs)/M(i);
plot(Mx(i),My(i),'o','MarkerSize',5,'MarkerFaceColor','yellow')
text(Mx(i),My(i),'M')

if ((0 < My(i)) && (My(i) <= tiexc(1)))
    slope = 0 + (My(i)-0)/(tiexc(1)-0)*tie_slope(1);
elseif((tiexc(1) < My(i)) && (My(i) <= tiexc(2)))
    slope= tie_slope(1) + (My(i)-tiexc(1))/(tiexc(2)-tiexc(1))*tie_slope(2);
elseif((tiexc(2) < My(i)) && (My(i) <= tiexc(3)))
    slope= tie_slope(2) + (My(i)-tiexc(2))/(tiexc(3)-tiexc(2))*tie_slope(3);
elseif((tiexc(3) < My(i)) && (My(i) <= tiexc(4)))
    slope= tie_slope(3) + (My(i)-tiexc(3))/(tiexc(4)-tiexc(3))*tie_slope(4);
elseif((tiexc(4) < My(i)) && (My(i) <= tiexc(5)))
    slope= tie_slope(4) + (My(i)-tiexc(4))/(tiexc(5)-tiexc(4))*tie_slope(5);
elseif((tiexc(5) < My(i)) && (My(i) <= tiexc(6)))
    slope= tie_slope(5) + (My(i)-tiexc(5))/(tiexc(6)-tiexc(5))*tie_slope(6);
elseif((tiexc(6) < My(i)) && (My(i) <= tiexc(7)))
    slope= tie_slope(6) + (My(i)-tiexc(6))/(tiexc(7)-tiexc(6))*tie_slope(7);
elseif((My(i) > tiexc(7)))
    slope= tie_slope(7) + (My(i)-tiexc(7))/(0.5-tiexc(7))*1;
end

%%tie line and the polynomial solving together
syms x y
[A,B] = vpasolve([y==poly2sym(p1), y==(My(i) +slope*(x - Mx(i)) )],[x,y],[0 0.24;0 0.346]);
xbr(i)=A(1);
xcr(i)=B(1);
plot(xbr(i),xcr(i),'o','MarkerSize',5,'MarkerFaceColor','yellow',Color='black');
text(xbr(i),xcr(i),'R');
syms x y
[ybe(i),yce(i)] = vpasolve([y == poly2sym(p2),y == My(i) + slope*(x - Mx(i))],[x,y],[0.24 1; 0 0.346]);
plot(ybe(i),yce(i),'o','MarkerSize',5,'MarkerFaceColor','yellow',Color='black')
text(ybe(i),yce(i),'E')

plot([xbr(i) ybe(i)],[xcr(i) yce(i)],Color='black',LineWidth=1.0) %plotting the the required tie line.

R = [xbr(i) xcr(i)];
E = [ybe(i) yce(i)];

syms E R
[e(i),r(i)] = vpasolve(E == M(i)-R, E == R*((My(i)-xcr(i))/(yce(i)-My(i))),[E,R]); %calculation of E and R using equation 8.3 and 8.7 BKD 
disp(r(i))
disp(e(i))

solute_fraction = (yce(i)*e - ycs*S)/(F*xcf); %fraction of solute saparated.

F = r(i);
xbf = xbr(i);
xcf = xcr(i);
plot([xbf ybs],[xcf ycs])

%%finding of the of required tie line-------------------
while xcr(i)>0.03
    i=i+1;
    M(i) = F + S;
    Mx(i) = (F*xbf + S*ybs)/M(i);
    My(i) = (F*xcf + S*ycs)/M(i);
    plot(Mx(i),My(i),'o','MarkerSize',5,'MarkerFaceColor','yellow')
    text(Mx(i),My(i),'M')

    if ((0 < My(i)) && (My(i) <= tiexc(1)))
        slope = 0 + (My(i)-0)/(tiexc(1)-0)*tie_slope(1);
    elseif((tiexc(1) < My(i)) && (My(i) <= tiexc(2)))
        slope= tie_slope(1) + (My(i)-tiexc(1))/(tiexc(2)-tiexc(1))*tie_slope(2);
    elseif((tiexc(2) < My(i)) && (My(i) <= tiexc(3)))
        slope= tie_slope(2) + (My(i)-tiexc(2))/(tiexc(3)-tiexc(2))*tie_slope(3);
    elseif((tiexc(3) < My(i)) && (My(i) <= tiexc(4)))
        slope= tie_slope(3) + (My(i)-tiexc(3))/(tiexc(4)-tiexc(3))*tie_slope(4);
    elseif((tiexc(4) < My(i)) && (My(i) <= tiexc(5)))
        slope= tie_slope(4) + (My(i)-tiexc(4))/(tiexc(5)-tiexc(4))*tie_slope(5);
    elseif((tiexc(5) < My(i)) && (My(i) <= tiexc(6)))
        slope= tie_slope(5) + (My(i)-tiexc(5))/(tiexc(6)-tiexc(5))*tie_slope(6);
    elseif((tiexc(6) < My(i)) && (My(i) <= tiexc(7)))
        slope= tie_slope(6) + (My(i)-tiexc(6))/(tiexc(7)-tiexc(6))*tie_slope(7);
    elseif((My(i) > tiexc(7)))
        slope= tie_slope(7) + (My(i)-tiexc(7))/(0.5-tiexc(7))*1;
    end

    %%tie line and the polynomial solving together
    syms x y
    [C,D] = vpasolve([y==poly2sym(p1), y==(My(i) +slope*(x - Mx(i)) )],[x,y],[0 0.24;0 0.346]);
    xbr(i)=C(1);
    xcr(i)=D(1);
    plot(xbr(i),xcr(i),'o','MarkerSize',5,'MarkerFaceColor','yellow',Color='black');
    text(xbr(i),xcr(i),'R');
    syms x y
    [ybe(i),yce(i)] = vpasolve([y == poly2sym(p2),y == My(i) + slope*(x - Mx(i))],[x,y],[0.24 1; 0 0.346]);
    plot(ybe(i),yce(i),'o','MarkerSize',5,'MarkerFaceColor','yellow',Color='black')
    text(ybe(i),yce(i),'E')

    plot([xbr(i) ybe(i)],[xcr(i) yce(i)],Color='black',LineWidth=1.0) %plotting the the required tie line.

    R = [xbr(i) xcr(i)];
    E = [ybe(i) yce(i)];

    syms E R
    [e(i),r(i)] = vpasolve(E == M(i)-R, E == R*((My(i)-xcr(i))/(yce(i)-My(i))),[E,R]); %calculation of E and R using equation 8.3 and 8.7 BKD 
    disp(r(i))
    disp(e(i))

    solute_fraction = (yce(i)*e - ycs*S)/(F*xcf); %fraction of solute saparated.

    F = r(i);
    xbf = xbr(i);
    xcf = xcr(i);
    plot([xbf ybs],[xcf ycs])
    
    %plot([double(xbr(i)) double(ybe(i))],[double(xcr(i)) double(yce(i))],'o-','Color',[0,0.25,1],'linewidth',0.35)
    %plot(double(Mx(i)),double(My(i)),'bo')

end
overall_frac = (1000*0.35 - xcr(3)*r(3))/(1000*0.35)
fprintf('The minimum number of stages is %d\n',i);

hold on
plot([0 1 0 0],[0 0  1 0])
xlabel('xB,yB');ylabel('xC,yC');title('interpolated data - overall process schematic')
text(0,0.35,'F')
text(0.98,0.02,'- S')