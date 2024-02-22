clc; clear;
%Problem1
%for figure 1
B = [0 0.004 0.006 0.01 0.02 0.03 0.036 0.07 0.13 0.16 0.19 0.23 0.26 0.5 0.63 0.71 0.78 0.84 0.9 0.95 1];
C = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.4025 0.405 0.402 0.4 0.35 0.3 0.25 0.2 0.15 0.1 0.05 0];
% figure 1
figure(1)
plot(B,C,'bo-');grid on;hold on;
plot([0 1 0 0],[0 0 1 0],'k-.','linewidth',1.25)
xlabel('xB');ylabel('xC,yC');grid on;
title('raw data')
p1 = polyfit(B,C,4);

X_=0:0.02:1;
f1 = polyval(p1,X_);
plot(X_,f1,'black')
disp("coefficient of polynomials are");
for i=1:5
    disp(f1(i));
end

%for figure 2
 %% Example tie lines 
tiexb =[ 0.0032    0.0053    0.0084    0.0229    0.1229]
tieyb = [ 0.9659    0.9329    0.9109    0.8471    0.6065]
tiexc= polyval(p1,tiexb);
tieyc = polyval(p1,tieyb);
%%slop of tie lines
for i=1:5
    disp("slopes of tie line are");
    disp((tieyc-tiexc/tieyb-tiexb));
end


 %% Plotting the tie lines
for i = 1:length(tiexc)
plot([tiexb(i) tieyb(i)], [tiexc(i) tieyc(i)],'mo:','linewidth',1.25);
end

%% Slopes of the tie lines
tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc)
tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i));
end

