B = [0 0.0011 0.0021 0.0031 0.0058 0.0136 0.0372 0.11 0.1259 0.1376 0.1508 0.16 0.3357 0.4751 0.608 0.7698 0.8872 0.9982 1];
C = [0 0 0.1 0.2 0.3 0.4 0.5 0.5821 0.6 0.6058 0.6107 0.6093 0.5819 0.4944 0.3748 0.2223 0.1079 0 0];
figure(1)
plot(B,C,'bo-');
grid on;
hold on;
plot([0 1],[1 0],'k-','linewidth',1.25)
plot([0 0],[1 0],'k-','linewidth',1.25)
plot([0 1],[0 0],'k-','linewidth',1.25)
xlabel("xB, yB");
ylabel("xC, yC");
title("raw data");

tiexc = [0.0 0.2 0.4 0.6];
tieyc = [0.0 0.2223 0.4944 0.6107];
tiexb = [0.0011 0.0031 0.0136 0.1259];
tieyb = [0.9982 0.7698 0.4751 0.1508];

for i = 1:length(tiexc)
plot([tiexb(i) tieyb(i)], [tiexc(i) tieyc(i)], 'mo:', 'linewidth', 1.25);
end

tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc)
    tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i));
end

My = 0.5;
Mx = 0.2;
if (0 < My) && (My <= 0.2)
slope = 0 + (My - 0) * tie_slope(1) / 0.2;
elseif (0.2 < My) && (My <= 0.4)
slope = tie_slope(2) + (My - 0.2) * (tie_slope(3) - tie_slope(2)) / (0.4 - 0.2);
elseif (0.4 < My) && (My <= 0.6)
slope = tie_slope(3) + (My - 0.4) * (tie_slope(4) - tie_slope(3)) / (0.6 - 0.4);
elseif (My > 0.6)
slope = tie_slope(4) + (My - 0.6) * 0.25 / (My - 0.6);
end

x_vals = linspace(0, 1, 10);
y_vals = slope * (x_vals - Mx) + My;
plot(x_vals, y_vals, 'r-', 'linewidth', 1.5);
plot(Mx, My, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'g');

p1 = polyfit(B(1:10), C(1:10), 5);
f1 = polyval(p1, B(1:10));

p2 = polyfit(B(11:19), C(11:19), 3);
f2 = polyval(p2, B(11:19));

figure(2)
plot(B(1:10), f1, 'g');
hold on;
plot(B(11:19), f2, 'r');
hold on;
scatter(B,C,'o');
legend("raffinate", "extract", "experimental data");
title("curve fitting");
xlabel("xB, yB");
ylabel("xC, yC");

figure(3)
f = [f1 f2];
plot(C, f, 'b');
xlabel("experimental");
ylabel("simulated");
title("goodness of fit");

sum1=0;
for i=1:10
    sum1=sum1+(f1(i) - C(i))^2;
end
rmse = sqrt((sum1)/10);
disp(rmse);

sum2=0;
for i=11:19
    sum2=sum2+(f2(i-10) - C(i))^2;
end
rmse1 = sqrt((sum2)/9);
disp(rmse1);

figure(4)
plot(B,C,'bo-');
grid on;
hold on;
plot(B(1:10), f1, 'g');
hold on;
plot(B(11:19), f2, 'r');
hold on;
legend("eqlb curve", "raffinate", "extract", "experimental data");
plot([0 1],[1 0],'k-','linewidth',1.25)
plot([0 0],[1 0],'k-','linewidth',1.25)
plot([0 1],[0 0],'k-','linewidth',1.25)
xlabel("xB, yB");
ylabel("xC, yC");
title("fitted data");