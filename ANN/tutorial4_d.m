B = [0 0.0011 0.0021 0.0031 0.0058 0.0136 0.0372 0.11 0.1259 0.1376 0.1508 0.16 0.3357 0.4751 0.608 0.7698 0.8872 0.9982 1];
C = [0 0 0.1 0.2 0.3 0.4 0.5 0.5821 0.6 0.6058 0.6107 0.6093 0.5819 0.4944 0.3748 0.2223 0.1079 0 0];

tiexc = [0, 0.2, 0.4, 0.6];
tieyc = [0.0, 0.2223, 0.4944, 0.6107];
tiexb = [0.0011, 0.0031, 0.0136, 0.1259];
tieyb = [0.9982, 0.7698, 0.4751, 0.1508];

tie_slope = zeros(1,length(tiexc));
for i = 1:length(tiexc)
tie_slope(i) = (tieyc(i) - tiexc(i))/(tieyb(i) - tiexb(i));
end

F = [0, 0.5];
S = [0.995, 0.005];

p1 = polyfit(B(1:10),C(1:10),5);
p2 = polyfit(B(10:19),C(10:19),3);

solute_removed = zeros(1,80);
k = 1;
for j = 400:10:1200
    stages = 3;
    S = [j j j];
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
    for i = 1:3
        M(i) = F + S(i);
        My(i) = (F*xcf + S(i)*ycs)/M(i); % y coordinate of the mixture point
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

        [xbr(i),xcr(i)] = vpasolve([y == poly2sym(p1),y == My(i) + slope*(x - Mx(i))],[x,y], [0 0.1508; 0 0.6107]); % for the raffinate part

        [ybe(i),yce(i)] = vpasolve([y == poly2sym(p2), y == My(i) + slope*(x-Mx(i))],[x,y],[0.16 1; 0.6093 0]); % for the extract part

        [r(i),e(i)] = solve([R + E - M(i),((R*xcr(i)) + (E*yce(i)) - (M(i)*My(i))),R>0,E>0],[R,E]);
        F = r(i); % This is the feed in the subsequent stage
        xcf = xcr(i); % Respective concentration
        solute_removed(k) = yce(3)*e(3);
    end
    k = k+1;
end

figure(4)
x = [400:10:1200];
plot(x,(500-solute_removed)/500, '-bo');
xlabel("amount of solvent");
ylabel("Solute level");
title("Solvent level vs no of stages");