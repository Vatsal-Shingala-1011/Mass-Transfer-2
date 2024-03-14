clear;
clc;

MSI = [1.41 1.19 0.841 0.595 0.42 0.297];
MSF = [1.19 0.841 0.595 0.42 0.297 0];
M = [28.5 29.2 37.5 27 24.7 3.1];
sf = 0.42;
d = 1400;


L = (10^-1)*(MSI+MSF)/2;
L_diff = (MSI - MSF)*10^-1;

for i = 1:6
n(i) = M(i)*10^-2/ (d*sf*L(i)^3*L_diff(i));
end

scatter(L, log(n));
xlabel("L");
ylabel("ln(n)");
title("ln(n) vs avg particle size");
hold on;

p = polyfit(L, log(n), 1);
f1 = polyval(p, L);
plot(L, f1);
rmse = rmse(f1, log(n))
t = 200/250;
G = 1/p(1)*t
n0 = exp(p(2))