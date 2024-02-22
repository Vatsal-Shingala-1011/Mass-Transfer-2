%Determination of heat of adsorption from Clausius-Clayperon equation:
% delta_H_iso = -R*((T)^2)*d(lnp)/dT = -R*d(lnp)/d(1/T)
%Fitting of Langmuir isotherm and determination of heat of adsorption

%At T = 310.9K
p1 = [15.6 31.74 59.6 99.97 293 479.2 679.1]; %kPa
q1 = [2.819 3.48 3.97 4.342 4.94 5.294 5.304]; %mmol/g
y1 = p1./q1;
x1 = p1;
lnp1 =log(p1);

%At T = 338.7K
p2 = [27.2 53.2 99.97 244.8 424 617.1 803.2]; %kPa
q2 = [2.469 3.078 3.635 4.188 4.475 4.71 5.289]; %mmol/g
y2 = p2./q2;
x2 = p2;
lnp2 =log(p2);
%At T = 366.5K
p3 = [21.07 45.33 93.34 279.2 461.9 634.3]; %kPa
q3 = [1.677 2.386 2.954 3.584 3.922 4.224]; %mmol/g
y3 = p3./q3;
x3 = p3;
lnp3 =log(p3);

pp1=polyfit(p1,p1./q1,1); %x,y
pp2=polyfit(p2,p2./q2,1);
pp3=polyfit(p3,p3./q3,1);

display(p1(1))

% Langmuir isotherm
%------ q = (qm*k*p)/(1+k*p)
%------ p/q = 1/k*qm + (1/qm)*p
%------ y = c + m*x
%------ m =1/qm
%------ c = 1/k*qm

m1=pp1(1);
m2=pp2(1);
m3=pp3(1);
c1=pp1(2);
c2=pp2(2);
c3=pp3(2);

qm1=1/m1;
qm2=1/m2;
qm3=1/m3;

k1=1/(c1*qm1);
k2=1/(c2*qm2);
k3=1/(c3*qm3);
qexp1=zeros(7);
for i=1:7
    qexp1(i)=(qm1*k1*p1(i))/(1+(k1*p1(i)));
end

figure(1)
f1 = plot(qexp1, q1, 'bo-'); grid on; hold on;
xlabel('q expw');
ylabel('q simual');