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
