B = [0 0.0011 0.0021 0.0031 0.0058 0.0136 0.0372 0.11 0.1259 0.1376 0.1508 0.16 0.3357 0.4751 0.608 0.7698 0.8872 0.9982 1];
C = [0 0 0.1 0.2 0.3 0.4 0.5 0.5821 0.6 0.6058 0.6107 0.6093 0.5819 0.4944 0.3748 0.2223 0.1079 0 0];
figure(1)
plot(B,C,'bo-');
grid on;
hold on;
plot([0 1 0 0],[0 0 1 0], 'k-.', 'linewidth',1.25);
xlabel("xB, yB");
ylabel("xC, yC");
title('raw data');

tiexc = [0, 0.2, 0.4, 0.6];
tieyc = [0.0, 0.2223, 0.4944, 0.6107];
tiexb = [0.0011, 0.0031, 0.0136, 0.1259];
tieyb = [0.9982, 0.7698, 0.4751, 0.1508];

for i = 1:length(tiexc)
plot([tiexb(i) tieyb(i)], [tiexc(i) tieyc(i)],'mo:','linewidth',1.25);
end

tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc)
tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i));
end

F = [0, 0.5];
S = [0.995, 0.005];
plot(S,F,'g','linewidth',0.35)
text(0, 0.5,'F');
text(0.995, 0.005,'S');

p1 = polyfit(B(1:10),C(1:10),5);
p2 = polyfit(B(10:19),C(10:19),3);

stages = 10;
S = repmat(800, 1, 10);
F = 1000;
xbf = 0;
xcf = 0.5;
ybs = 0.995;

ycs = 0.005;
r = ones(1,stages);
e = ones(1,stages);
xbr = ones(1,stages);xcr = ones(1,stages);
ybe = ones(1,stages);yce = ones(1,stages);
M = ones(1,stages);
Mx = ones(1,stages);
My = ones(1,stages);
syms R E
syms x y

for i = 1:10
    M(i) = F + S(i);
    My(i) = (F*xcf + S(i)*ycs)/M(i); 
    Mx(i) = (F*xbf + S(i)*ybs)/M(i); % x coordinate of the mixture point
    
    if (0.1 < My(i)) && (My(i) <= 0.3)
    slope = tie_slope(1) + ( My(i)-0.1) * (tie_slope(2) - tie_slope(1)) / (0.3 - 0.1 );
    elseif (0.3 < My(i)) && (My(i) <= 0.5)
    slope = tie_slope(2) + (My(i) - 0.3) * (tie_slope(3) - tie_slope(2)) / (0.5 - 0.3);
    elseif (0.5 < My(i)) && (My(i)<=0.6)
    slope = tie_slope(3) + (My(i) - 0.3) * (tie_slope(4) - tie_slope(3)) / (0.5 - 0.3);
    elseif (My(i)>0.6)
    slope = tie_slope(4) + (My(i) - 0.3) * 0.12 / (0.6107 - 0.5);
    end
    
    [xbr(i),xcr(i)] = vpasolve([y == poly2sym(p1),y == My(i) + slope*(x - Mx(i))],[x,y],[0 0.1508; 0 0.6107]); % for the raffinate part
    
    [ybe(i),yce(i)] = vpasolve([y == poly2sym(p2), y == My(i) + slope*(x-Mx(i))],[x,y],[0.16 1; 0.6093 0]); % for the extract part
   

    [r(i),e(i)] = solve([R + E - M(i),((R*xcr(i)) + (E*yce(i)) - (M(i)*My(i))),R>0,E>0],[R,E]);
    F = r(i); 
    xcf = xcr(i); 

end

figure(3)
plot(yce, '-bo');
xlabel("Number of stages");
ylabel("Solute Fraction");
title("Solute Fraction vs no of stages");
display("From the graph we can say 4 stages required for 96% removal of the solute for the second solvent scheme.")

% From the graph we can say 4 stages required for 96% removal of the solute
% for the second solvent scheme.