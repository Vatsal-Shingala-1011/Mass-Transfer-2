p=[10 20 40 100 150 200 250 300 350 400 500]
q=[71.3 142.2 287.7 679.4 1025 1053 1175 1316 1996 3451 5283]
pv=760
chi=p/pv;

for i=1:11
    kappa(i)= p(i)/(q(i)*(pv-p(i)));
end

p=polyfit(chi,kappa,1);
f1=polyval(p,chi);

%paramerer esstimation


c=1+(p(1)/p(2));
qm=1/(p(2)*c);

figure(3)
plot(chi,q,'o');
%hold on;
%plot(q,q,'r');
%xlabel("q experimental");
%ylabel("q theoretical");
%hold off;