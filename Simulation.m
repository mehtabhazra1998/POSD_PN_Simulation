clc;clear;close all;

T=300; % Temperature in Kelvin
k=8.617e-5; % Boltzmann constant (eV/K)
e0=8.85e-14; % permittivity of free space (F/cm)
q=1.602e-19; % charge on an electron (coul)
KS=11.8; % Dielectric constant of Si
ni=1.0e10; % intrinsic carrier conc. in Silicon at 300K
EG=1.12; % Silicon band gap (eV)

xleft = -3.5e-4; % Leftmost x position
xright = -xleft; % Rightmost x position
NA=10E15;
ND=10E15;
VA=0.2;


%Input Parameters
disp("The simulation is starting. Pls enter the input parameters:")
ND = input("Donor concentration per cc: ");
NA = input("Acceptor concentration per cc: ");
ni = input("intrinsic carrier concentration per cc: ");
T = input("Simulation Temperature in Kelvin: ");
VA = input("Biasing Voltage in Volts: ");
EG = input("Energy Bandgap in ev:");

%er = input("Relative permittivity of the semiconductor: ");
%length = input("Length of the PN junction in metres: ");

Vbi = k*T*log((NA *ND)/ni^2);
xN=sqrt(2*KS*e0/q*NA*(Vbi-VA)/(ND*(NA+ND))); % Depletion width n-side
xP=sqrt(2*KS*e0/q*ND*(Vbi-VA)/(NA*(NA+ND))); % Depletion width p-side
x = linspace(xleft, xright, 200);

divD = q*ND/(2*KS*e0);
divA = q*NA/(2*KS*e0);
V_prime = Vbi - VA;
Vx1 = (V_prime - divD*((xN-x).^2)).*(x>=0 & x<=xN);
Vx2 = divA*((xP+x).^2).*(x>=-xP & x<0);
Vx3 = V_prime.*(x>=xN);

Vx=Vx1 + Vx2 +Vx3; 
axis ([xleft xright -1 1]);
title("voltage as a function of x")
plot(x,Vx);

figure();
title("Electric field")
plot(x(1:end-1),-1*diff(Vx));
figure()
q1 = +q*ND.*(x>=0 & x<=xN);
q2 = -q*NA.*(x>=-xP & x<0);
q=q1+q2;
axis ([xP-0.1 xN+0.1 -5 5]);
plot(x,q);
figure();

VMAX = 3; 
if(VA~=0)
    EF=Vx(1)+VMAX/2-k*T*log(NA/ni)+VA; % Fermi level
    EF2=Vx(1)+VMAX/2-k*T*log(ND/ni); 
end    
if(VA==0)
    EF=Vx(1)+VMAX/2-k*T*log(NA/ni); % Fermi level
    EF2=Vx(1)+VMAX/2-k*T*log(NA/ni); 
end
    
%Plot Diagram
title("Energy Band Diagram")
plot (x, -Vx+EG/2+VMAX/2);


axis ([xleft xright -VMAX VMAX]);
axis ('off');
hold on
plot (x, -Vx-EG/2+VMAX/2);
plot (x, -Vx+VMAX/2,'w:');
plot (x,(-Vx+EG/2+VMAX/2+-Vx-EG/2+VMAX/2)/2,'--');
plot ([xP+xP/4 xright], [EF EF]);
plot ([xleft xN-3e-05], [EF2 EF2]);
plot ([0 0], [0.15 VMAX-0.5], 'w--');
%plot(x,EF);
xstr = (- Vx(1)+ EG/2 + VMAX/2);
xstr2=(- Vx(1)-EG/2+ VMAX/2);
xstr3=(- Vx(1)+ VMAX/2);
xstr4=EF;
xstr5=(- Vx(200)+ EG/2 + VMAX/2);
xstr6=(- Vx(200)- EG/2 + VMAX/2);

text(xleft* 1.17,(- Vx(1)+ EG/2 + VMAX/2-.1),'Ec=');
text(xleft* 1.08,(- Vx(1)+ EG/2 + VMAX/2-.1),num2str(xstr));

text(xright* 1.02,(- Vx(200)+ EG/2 + VMAX/2-.05),'Ec=');
text(xright* 1.11,(- Vx(200)+ EG/2 + VMAX/2-.05),num2str(xstr5));

text(xleft* 1.17,(- Vx(1)-EG/2+ VMAX/2-.1),'Ev=');
text(xleft* 1.08,(- Vx(1)-EG/2+ VMAX/2-.1),num2str(xstr2));

text(xright* 1.02,(- Vx(200)-EG/2+ VMAX/2-.05),'Ev=');
text(xright* 1.11,(- Vx(200)-EG/2+ VMAX/2-.05),num2str(xstr6));

text(xleft* 1.17,(- Vx(1)+ VMAX/2-.1),'Ei =');
text(xleft* 1.08,(- Vx(1)+ VMAX/2-.1),num2str(xstr3));

text(xright* 1.02,(- Vx(200)+ VMAX/2-.1),'Ei=');
text(xright* 1.11,(- Vx(200)+ VMAX/2-.1),num2str((- Vx(200)+ VMAX/2-.1)));

text(xright* 1.02, EF-0.15,'EFn');


text(xleft*1.17,EF2+0.15,'EFp');


set(gca, 'DefaultTextUnits' , 'normalized ')
text(.18, 0,'p-side');
text(.47, 0, 'x=0');
text(.75, 0,'n-side'); 
set(gca,'DefaultTextUnits','data')
hold off
