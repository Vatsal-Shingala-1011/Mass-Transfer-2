% Question 1
MSI = [1.41, 1.19, 0.841, 0.595, 0.42, 0.297]
MSF = [1.19, 0.841, 0.595, 0.42, 0.297, 0];  

L_avg = (MSI + MSF)*(0.1)/2
L_diff = (MSI - MSF)*(0.1)

W_del = [28.5, 29.2, 37.5, 27, 24.7, 3.1]

n = W_del./(1400*0.42.*L_avg.*L_avg.*L_avg.*L_diff)

log_n = log(n)

X = L_avg
Y = log_n


plot(L_avg, log_n, 'bo-')


coefficients = polyfit(X, Y, 1);
fit_line = polyval(coefficients, X);
scatter(X, Y, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
hold on;
plot(X, fit_line, 'r-', 'LineWidth', 2);

xlabel('Lavg');
ylabel('ln(n)');
slope = coefficients(1)
y_intercept = coefficients(2)
title('ln(n) VS Lavg');

m = slope
tau = 200/250

G = -1/m*tau



