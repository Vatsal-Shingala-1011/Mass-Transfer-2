% Step 1: Generate equilibrium data
x = 0:0.01:0.2;  
y = 3.7* x;       

F=1000;
x0=0.15;
y0=0;
R=850;
X0=x0/(1-x0);
Y0=y0/(1-y0);

x_init=X0;
y_init=0;
E=100:10:500;
slope=-(R./E);

solue_remove_pre=zeros(1,length(slope));
%  figure
for i=1:length(slope)
    x_init=X0;
    y_init=0;
     figure
    for j=1:3
        % Step 2: Plot concentration of solute in extract vs. concentration of solute in raffinate
         plot(x, y, 'LineWidth', 1.5); hold on;
        xlabel('Concentration of Solute C in Raffinate');
        ylabel('Concentration of Solute C in Extract');
        title('Equilibrium Data: Concentration of Solute in Extract vs. Raffinate');
        grid on;
        hold on;
        plot(x_init, y_init, 'ro', 'MarkerSize', 8);  % Plot the point X0,Y0  X1',Y1' .....
        hold on;
    
        syms xn_d yn_d;
        eq1 = eq(yn_d, 3.75 * xn_d);
        eq2 = eq(yn_d, y_init + slope(i) * (xn_d - x_init));
    
        [xn_d,yn_d] = solve([eq1, eq2], [xn_d, yn_d]);
    
        xn_d = double(xn_d);
        yn_d = double(yn_d);
         plot([x_init,xn_d],[y_init,yn_d],'m-');
        x_init=xn_d;
        y_init=0;
         plot([x_init,xn_d],[y_init,yn_d],'m-');
    end
    solue_remove_pre(i) = (1-(xn_d*R/150))*100;
end
solvent_amount=E;
figure
plot(solvent_amount,solue_remove_pre);
xlabel('Solvent amount/stage');
ylabel('Percentage removal of solute');
title('Solvent amount/stage vs Percentage removal of solute');
grid on;

desired_y=95;
%estimate function so can calculate relation B at 95% c removed
figure
f1=polyfit(solvent_amount,solue_remove_pre,3);
p1=polyval(f1,desired_y);
x_values = linspace(100,500, 100);
y_estimated = polyval(f1, x_values);
plot(x_values, y_estimated, 'r-', 'LineWidth', 1.5);  % Plot estimated polynomial curve
xlabel('Solvent amount/stage');
ylabel('Estimated Percentage removal of solute');
title('Estimation of function between Solvent amount/stage and Percentage removal of solute');
hold off;

ep=0.01;
disp("Amount of B at 95% C removed is ")
for i=100:500
    if polyval(f1, i)<=95+ep && polyval(f1, i)>=95-ep
        disp(i);
        break;
    end
end
