%Determination of heat of adsorption from Clausius-Clayperon equation:
% delta_H_iso = -R*((T)^2)*d(lnp)/dT = -R*d(lnp)/d(1/T)
%Fitting of Langmuir isotherm and determination of heat of adsorption

%At T = 310.9K
p1 = [15.6 31.74 59.6 99.97 293 479.2 679.1]; %kPa
q1 = [2.819 3.48 3.97 4.342 4.94 5.294 5.304]; %mmol/g
y1 = p1./q1;
x1 = p1;

%At T = 338.7K
p2 = [27.2 53.2 99.97 244.8 424 617.1 803.2]; %kPa
q2 = [2.469 3.078 3.635 4.188 4.475 4.71 5.289]; %mmol/g
y2 = p2./q2;
x2 = p2;

%At T = 366.5K
p3 = [21.07 45.33 93.34 279.2 461.9 634.3]; %kPa
q3 = [1.677 2.386 2.954 3.584 3.922 4.224]; %mmol/g
y3 = p3./q3;
x3 = p3;

%plotiing p vs q
figure(1)
f1 = plot(p1, q1, 'bo-'); grid on; hold on;
f2 = plot(p2, q2, 'go-'); grid on; hold on;
f3 = plot(p3, q3, 'ro-'); grid on; hold on;
xlabel('Pressure of propane (kPa)');
ylabel('Amount of propane adsorbed at equilibrium (mmol/g carbon)');
title('Adsorption isotherm');
