clc; 
clear;
close all;

%Fitting of Langmuir isotherm and determination of heat of adsorption

%At T = 310.9K
p1 = [15.6 31.74 59.6 99.97 293 479.2 679.1];       %kPa
q1 = [2.819 3.48 3.97 4.342 4.94 5.294 5.304];      %mmol/g
y1 = p1./q1;
x1 = p1;
lnp1 = log(p1);

%At T = 338.7K
p2 = [27.2 53.2 99.97 244.8 424 617.1 803.2];       %kPa
q2 = [2.469 3.078 3.635 4.188 4.475 4.71 5.289];    %mmol/g
y2 = p2./q2;
x2 = p2;
lnp2 = log(p2);

%At T = 366.5K
p3 = [21.07 45.33 93.34 279.2 461.9 634.3];         %kPa
q3 = [1.677 2.386 2.954 3.584 3.922  4.224];        %mmol/g
y3 = p3./q3;
x3 = p3;
lnp3 = log(p3);
%%
%plotiing p vs q
figure(1)
f1 = plot(p1, q1, 'bo-'); grid on; hold on;
f2 = plot(p2, q2, 'go-'); grid on; hold on;
f3 = plot(p3, q3, 'ro-'); grid on; hold on;
xlabel('Pressure of propane (kPa)');
ylabel('Amount of propane adsorbed at equilibrium (mmol/g carbon)');
title('Adsoption isotherm');
legend([f1, f2, f3],'at T1','at T2','at T3','Location','southeast');
%%
% fitting a straight line through the points
P1 = polyfit(x1, y1,1);
P2 = polyfit(x2, y2,1);
P3 = polyfit(x3, y3,1);

% slopes and intercepts
m1 = P1(1);   
c1 = P1(2); 
m2 = P2(1);   
c2 = P2(2); 
m3 = P3(1);   
c3 = P3(2); 

% langmuir isotherm
%------ q = (qm*k*p)/(1+k*p)
%------ p/q = 1/k*qm + (1/qm)*p
%------ y = c + m*x
%------ m =1/qm
%------ c = 1/k*qm

disp('parameters')
qm1 = 1/m1
k1  = 1/(qm1*c1)
qm2 = 1/m2
k2  = 1/(qm2*c2)
qm3 = 1/m3
k3  = 1/(qm3*c3)
%%
figure(2)
plot(x1,y1,'b*'); grid on; hold on;
plot(x2,y2,'g*'); grid on; hold on;
plot(x3,y3,'r*'); grid on; hold on;
xlabel('p'); ylabel('p/q');
title('Fitting the Langmuir Isotherm');

% plotting the fitted curve
syms x
b1 = fplot(m1*x+c1);
b2 = fplot(m2*x+c2);
b3 = fplot(m3*x+c3);
legend([b1, b2, b3],'at T1','at T2','at T3','Location','southeast');
%%

%calculating simulated q
qsim1 = (qm1.*k1.*p1)./(1+k1.*p1);
qsim2 = (qm2.*k2.*p2)./(1+k2.*p2);
qsim3 = (qm3.*k3.*p3)./(1+k3.*p3);


%plotting q vs qsim
figure(3)
a1 = plot(q1, qsim1, 'bo-'); grid on; hold on;
a2 = plot(q2, qsim2, 'go-'); grid on; hold on;
a3 = plot(q3, qsim3, 'ro-'); grid on; hold on;
xlabel('q experimental'); ylabel('q simulated');
title('q experimental vs q simulated');
legend([a1, a2, a3],'at T1','at T2','at T3','Location','southeast');
%%
% Root mean square error
RMSE1 = rms(qsim1-q1) 
RMSE2 = rms(qsim2-q2) 
RMSE3 = rms(qsim3-q3) 

% mean 
Xbar1 = mean(q1);
Xbar2 = mean(q2);
Xbar3 = mean(q3);

Ybar1 = mean(qsim1);
Ybar2 = mean(qsim2);
Ybar3 = mean(qsim3);

% numerator for calculation of corelation coefficient
numerator1 = sum((q1-Xbar1).*(qsim1-Ybar1));
numerator2 = sum((q2-Xbar2).*(qsim2-Ybar2));
numerator3 = sum((q3-Xbar3).*(qsim3-Ybar3));

% denomenator for calculation of corelation coefficient
denominator1 = ((sum((q1-Xbar1).^2)).*(sum((qsim1-Ybar1).^2))).^0.5;
denominator2 = ((sum((q2-Xbar2).^2)).*(sum((qsim2-Ybar2).^2))).^0.5;
denominator3 = ((sum((q3-Xbar3).^2)).*(sum((qsim3-Ybar3).^2))).^0.5;

disp('correlation coefficient=')
correlation_coefficient1 = numerator1/denominator1
correlation_coefficient2 = numerator2/denominator2
correlation_coefficient3 = numerator3/denominator3

%%
% determination of heat of adsorption
% from Clausius-Clayperon equation:
% delta_H_iso = -R*((T)^2)*d(lnp)/dT = R*d(lnp)/d(1/T)

% fitting a straight line through the points
P_1 = polyfit(lnp1, q1,1);
P_2 = polyfit(lnp2, q2,1);
P_3 = polyfit(lnp3, q3,1);

% slopes and intercepts
m_1 = P_1(1);   
c_1 = P_1(2); 
m_2 = P_2(1);   
c_2 = P_2(2); 
m_3 = P_3(1);   
c_3 = P_3(2);

q = [1.5 2.5 3.5];
lnp1_sim = polyval(P_1, q);
lnp2_sim = polyval(P_2, q);
lnp3_sim = polyval(P_3, q);

figure(4)
plot(lnp1, q1,'b*'); grid on; hold on;
plot(lnp2, q2,'g*'); grid on; hold on;
plot(lnp3, q3,'r*'); grid on; hold on;
xlabel('lnp'); ylabel('q');
title('plot q vs lnp');

%plotting the fitted curve
d1 = fplot(m_1*x+c_1);
d2 = fplot(m_2*x+c_2);
d3 = fplot(m_3*x+c_3);
legend([d1, d2, d3],'at q = 1.5','at q = 2.5','at q = 3.5','Location','southeast');
%%
T = [310.9 338.7 366.5];

% plot 1/T vs lnp for calculation of slope to calculate isosteric heat of adsorption
figure(5)
n1 = plot(1./T, lnp1_sim,'bo-'); grid on; hold on;
n2 = plot(1./T, lnp2_sim,'go-'); grid on; hold on;
n3 = plot(1./T, lnp3_sim,'ro-'); grid on; hold on;
xlabel('1/T'); ylabel('lnp obtained');
title('plot 1/T vs lnp obtained ');
legend([n1, n2, n3],'at q = 1.5','at q = 2.5','at q = 3.5','Location','southwest');
%%
P_iso1 = polyfit(1./T, lnp1_sim ,1);
P_iso2 = polyfit(1./T, lnp2_sim ,1);
P_iso3 = polyfit(1./T, lnp3_sim ,1);

slope1 = P_iso1(1);
slope2 = P_iso2(1);
slope3 = P_iso3(1);

slope = (slope1+slope2+slope3)/3;

% ------delta_H_iso = R * slope

disp('isosteric heat of adsorption in cal/gmol');
Delta_H = 1.987 * slope