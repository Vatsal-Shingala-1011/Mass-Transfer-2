clc
clear
close all


%% Step 1

B = [0.04 0.05 0.07 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.993 1];
C = [0 0.11 0.26 0.375 0.474 0.487 0.468 0.423 0.356 0.274 0.185 0.09 0.001 0];

p1 = polyfit(B(1:7),C(1:7),5);
f1 = polyval(p1,B(1:7));


p2 = polyfit(B(7:14),C(7:14),3);
f2 = polyval(p2,B(7:14));
tiexc = [0.1 0.245 0.426];
tiexb = [0.048 0.065 0.133];
tieyc = [0.098 0.242 0.409];
tieyb = [0.891 0.736 0.523];
lt = length(tiexc);
slope_tie = zeros(1,lt);


for i = 1:lt
slope_tie(i) = (tieyc(i)-tiexc(i))/((tieyb(i)-tiexb(i)));

end

B1 = zeros(1,10001);
C1 = zeros(1,10001);
for i = 1:10001
B1(i) = (i-1)/10000;
C1(i) = round(interp1(B,C, (i-1)/10000),4);
end


tiexb = zeros(1,length(tiexc));
tieyb = zeros(1,length(tiexc));
for i = 1:length(tiexc)
index1 = find(abs(C1-tiexc(i)) < 0.001,1,'first');
index2 =  find(abs(C1-tieyc(i)) < 0.001,1,'last');
XX = [B1(index1) B1(index2)] ;
YY = [tiexc(i) tieyc(i)];

end

% DPH = [0.04,0.08,0.015];
% stage = [1,2,3];
% fr_dph_r = ones(1,3);
% fr_dph_e = ones(1,3);
% r = ones(1,3);
% e = ones(1,3);
dataP = 25;
data = ones(dataP,4);

%% Step 2 and 3
for o=1:1:dataP
  
delx = ones(1,dataP);
dely = ones(1,dataP);
stageInp = randi([2,10],1);
S = randi([500,2500],1);

F=2000;

xbf= 0;

xcfInp = rand(1);
if xcfInp==0.0||xcfInp==1.0
    xcfInp = 0.01;
end

xcf = xcfInp;

ybs=1;
ycs=0;
M = F + S;
Mx = (F*xbf + S*ybs)/M;
My = (F*xcf + S*ycs)/M;

Rnx = 0.03;
index1 = find(abs(C1-Rnx) < 0.001,1,'first');
XX = B1(index1);
YY = C1(index1);




%% step 4 and 5 - Plotting E1

syms x y
slope = ((My - YY)/(Mx-XX));
[E1x,E1y] = vpasolve([y == poly2sym(p2),y == My + slope*(x - Mx)],[x,y],[0 1; 0 1]);

%% step 6 - locating delta point

% F - xbf,xcf , S- xcs,ycs , R - XX,YY  
slope_FE1 = ((E1y-xcf)/(E1x-xbf));
slope_RS = ((YY-ycs)/(XX-ybs));

syms x y
[delx(o),dely(o)] = vpasolve([y==E1y + slope_FE1*(x-E1x),y==YY + slope_RS*(x-XX)],[x,y]);



%% step 7 - R1

Ex = ones(1,100);
Ey = ones(1,100);
Rx = ones(1,100);
Ry = ones(1,100);

Ex(1) = E1x;
Ey(1) = E1y;
Ex(2) = E1x;
Ey(2) = E1y;
count = 2;

for p = 1:1:stageInp
   
if ((0 < Ey(count)) && (Ey(count) <= 0.098))
slope = 0 + (Ey(count) - 0)*slope_tie(1)/(0.098-0);
elseif((0.098< Ey(count)) && (Ey(count) <= 0.242))
slope = slope_tie(1) + (Ey(count) - 0.098)*slope_tie(2)/(0.242 - 0.098);
elseif ((0.242 < Ey(count)) && (Ey(count) <= 0.409))
slope = slope_tie(2) + (Ey(count) - 0.242)*slope_tie(3)/(0.409 - 0.242);
else
   slope = slope_tie(3) + (Ey(count)-0.409)*(-0.515)/(0.5-0.409);
end

syms x y
[C,D] = vpasolve([y==poly2sym(p1),y==Ey(count) + slope*(x-Ex(count))],[x,y],[0 0.3; 0 0.6]);
if isempty(C)||isempty(D)
    continue;
end
Rx(count)=C(1);
Ry(count)=C(1);
                      



slope_Rdel = ((Ry(count) - dely(o))/(Rx(count) - delx(o)));



syms x y 
[H,G] = vpasolve([y==poly2sym(p2),y== dely(o) + slope_Rdel*(x-delx(o))],[x,y],[0.5 1; 0 0.5]);
if isempty(H)||isempty(G)
    continue;
end
Ex(count + 1)=H(1);
Ey(count + 1)=G(1);

count = count+1;
Ry(count) = Ry(count-1);

end


% syms E R
% [e(o),r(o)] = vpasolve([E==R+del,E==],[E,R]);


R = ones(1, stageInp);
E = ones(1, stageInp);
syms x y
eqn11 = x == F + S - y;
eqn12 = x == (1/Ry(stageInp))*((S*ycs) + (F*xcf) - (y*Ey(1)));
[x,y] = vpasolve([eqn11, eqn12], [x, y]);
R(stageInp) = x(1);
E(1) = y(1);
frac = (R(stageInp)*Ry(stageInp)) - (S*ycs);
c = R(stageInp) - S;
for k = 2: stageInp
syms x y
    eqn13 = (x*Ry(k-1)) - (y*Ey(k)) == frac;
    eqn14 = x - y == c;
    [Q,P] = vpasolve([eqn13, eqn14], [x, y]);
    if isempty(P)||isempty(Q)
        continue;
    end
    R(k-1) = Q(1);
     E(k) = P(1);
end
percentageRemoved = ones(1, stageInp);
% percentageRemoved(1) = ((F*xcf) - (R(1)*Xcr(1)))/(F*xcf);
for m = 1 : stageInp
        percentageRemoved(m) = ((F*xcf) - (R(m)*Ry(m)))/(F*xcf);
end

%  disp('for a fraction of DPH of ')
%  disp(xcfInp) 
%  disp('No. of Stages: ');
%         disp(count-2);
        disp(percentageRemoved(stageInp));
      

        data(o,:) = [S stageInp xcfInp percentageRemoved(stageInp)*100];

end


try
    numbers = xlsread("counter_data.xlsx");
    lastRow = size(numbers, 1);
    nextRow = lastRow + 1;
    cellReference = sprintf('A%d:D%d', nextRow, nextRow+100);
    xlswrite("counter_data.xlsx",data, 'Sheet1', cellReference);
    catch
    disp("Can't write excel")
end




%fr_dph_r(1) = Ry(1)*r(1)/(F*xcf);
%fr_dph_e(1) = Ey(1)*e(1)/(F*xcf);
