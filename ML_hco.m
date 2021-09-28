
% Morris-Lecar half-center oscillator as in Skinner et al, (1996)
% Escape vs Release

% mex Morris_lecar_Skinner_HCO.cpp

TT=30000; dt=0.02; N=TT/dt;

gbarK=0.02; gbarCa=0.015; gl=0.005; gbarh=0; gsyn=0.01; % mS/cm^2
Iapp=0.8; % uA/cm^2

addpath('C:\Users\moroz\Documents\ML_Marder')

% synaptic escape
Vthe=-30; % mV
% control
gmi=0;
[V1e,m1e,V2e,m2e] = Morris_lecar_Skinner_HCO(TT,gmi,gbarK,gbarCa,gbarh,gl,gsyn,Vthe,Iapp);
% with gmi
gmi=0.005; % mS/cm^2
[V1e_mi,m1e_mi,V2e_mi,m2e_mi] = Morris_lecar_Skinner_HCO(TT,gmi,gbarK,gbarCa,gbarh,gl,gsyn,Vthe,Iapp);

% synaptic release
Vthr=20; % mV
%control
gmi=0;
[V1r,m1r,V2r,m2r] = Morris_lecar_Skinner_HCO(TT,gmi,gbarK,gbarCa,gbarh,gl,gsyn,Vthr,Iapp);
%with gmi
gmi=0.005; % mS/cm^2
[V1r_mi,m1r_mi,V2r_mi,m2r_mi] = Morris_lecar_Skinner_HCO(TT,gmi,gbarK,gbarCa,gbarh,gl,gsyn,Vthr,Iapp);

Ve = [V1e;V2e;V1e_mi;V2e_mi]; Vr = [V1r;V2r;V1r_mi;V2r_mi];


% phase planes ans traces

C=1; ECa=100; EK=-80; Eleak=-50; Eh=-20; Emi=50; Esyn=-80;
v1=0; v2=15; v3=0; v4=15; v5=78.3; v6=10.5; v7=-42.2; v8=87.3;
vhalf_mi=-10; k_mi=-10; phi = 0.0005;
V_slope=2;

gmi=0;
myfun1 = @(Vn,wn) (1/C)*(gbarCa*(0.5*(1+tanh((Vn-v1)/v2)))*(ECa-Vn)+gbarK*wn*(EK-Vn)+...
    gl*(Eleak-Vn)+ gsyn*(Esyn-Vn)+gmi*(1/(1+exp((Vn-vhalf_mi)/k_mi)))*(Emi-Vn)+Iapp);

myfun2 = @(Vn,wn) phi*cosh((Vn-v3)/(2*v4))*(0.5*(1+tanh((Vn-v3)/v4))-wn);

myfun3 = @(Vn,wn) (1/C)*(gbarCa*(0.5*(1+tanh((Vn-v1)/v2)))*(ECa-Vn)+gbarK*wn*(EK-Vn)+...
    gl*(Eleak-Vn)+gmi*(1./(1+exp((Vn-vhalf_mi)/k_mi)))*(Emi-Vn)+Iapp);

gmi=0.005;
myfun4 = @(Vn,wn) (1/C)*(gbarCa*(0.5*(1+tanh((Vn-v1)/v2)))*(ECa-Vn)+gbarK*wn*(EK-Vn)+...
    gl*(Eleak-Vn)+gsyn*(Esyn-Vn)+gmi*(1/(1+exp((Vn-vhalf_mi)/k_mi)))*(Emi-Vn)+Iapp);

myfun5 = @(Vn,wn) (1/C)*(gbarCa*(0.5*(1+tanh((Vn-v1)/v2)))*(ECa-Vn)+gbarK*wn*(EK-Vn)+...
    gl*(Eleak-Vn)+gmi*(1./(1+exp((Vn-vhalf_mi)/k_mi)))*(Emi-Vn)+Iapp);

clf
% escape
subplot(2,2,1)
plot(V1e(1,5000/dt:end),m1e(1,5000/dt:end),'color','k','linewidth',1.5), hold on
plot(V1e_mi(1,5000/dt:end),m1e_mi(1,5000/dt:end),'color','r','linewidth',1.5)
ylim([0 1.5]); xlim([-50 70]);
%
a1 = fimplicit(@(Vn,wn) myfun1(Vn,wn), [-80 70 -0.2 2])
set(a1,'color','k','linestyle','--', 'linewidth',1.5); hold on
h=display.plotvertline(Vthe); set(h,'linewidth',1,'color','g')

a2=fimplicit(@(Vn,wn) myfun2(Vn,wn), [-80 70 -0.2 2])
set(a2,'color','k','linewidth',1.5);

a3 = fimplicit(@(Vn,wn) myfun3(Vn,wn), [-80 70 -0.2 2])
set(a3,'color','k','linestyle','--','linewidth',2)

a4 = fimplicit(@(Vn,wn) myfun4(Vn,wn), [-80 70 -0.2 2])
set(a4,'color','r','linestyle','--','linewidth',2)

a5 = fimplicit(@(Vn,wn) myfun5(Vn,wn), [-80 70 -0.2 2])
set(a5,'color','r','linestyle','--','linewidth',2)

h=display.plotvertline(Vthe); set(h,'linewidth',1,'color','g')

ylim([0 1.5]); xlim([-50 70]);
xlabel('Voltage,mV'); ylabel('IK gating variable, w')
title('Synaptic escape')
set(gca,'Fontsize',14)

% plot traces for the escape case
for i=[1,2]
    subplot(4,2,3+i*2)
    plot(dt/1000:dt/1000:TT/1000,Ve(i,:),'k','linewidth',1.5), hold on
    plot(dt/1000:dt/1000:TT/1000,Ve(i+2,:),'r','linewidth',1.5)
    %hold on,  plot(dt/1000:dt/1000:TT/1000,m1e*20,'color','b','linewidth',1)
    display.plothorzline(Vthe);set(h,'linewidth',1,'color','g')
    xlabel('Time,s'); ylabel('Voltage,mV'); 
    xlim([10 30]); ylim([-50 70])
    set(gca,'Fontsize',14)
end

% release
subplot(2,2,2)
plot(V1r(1,5000/dt:end),m1r(1,5000/dt:end),'color','k','linewidth',1.5), hold on
plot(V1r_mi(1,5000/dt:end),m1r_mi(1,5000/dt:end),'color','r','linewidth',1.5), hold on

a1 = fimplicit(@(Vn,wn) myfun1(Vn,wn), [-80 70 -0.2 2])
set(a1,'color','k','linestyle','--', 'linewidth',1.5); hold on

a2=fimplicit(@(Vn,wn) myfun2(Vn,wn), [-80 70 -0.2 2])
set(a2,'color','k','linewidth',1.5);

a3 = fimplicit(@(Vn,wn) myfun3(Vn,wn), [-80 70 -0.2 2])
set(a3,'color','r','linestyle','--','linewidth',1.5)

a4 = fimplicit(@(Vn,wn) myfun4(Vn,wn), [-80 70 -0.2 2])
set(a4,'color','r','linestyle','--','linewidth',2)

a5 = fimplicit(@(Vn,wn) myfun5(Vn,wn), [-80 70 -0.2 2])
set(a5,'color','r','linestyle','--','linewidth',2)

h=display.plotvertline(Vthr); set(h,'linewidth',1,'color','g')
xlabel('Voltage,mV'); ylabel('IK gating variable, w')
title('Synaptic release')
set(gca,'Fontsize',14)
ylim([0 1.5]); xlim([-50 70]);

% plot traces for the release case
for i=[1,2]
    subplot(4,2,4+i*2)
    plot(dt/1000:dt/1000:TT/1000,Vr(i,:),'k','linewidth',1.5), hold on
    plot(dt/1000:dt/1000:TT/1000,Vr(i+2,:),'r','linewidth',1.5)
    display.plothorzline(Vthr)
    xlabel('Time,s'); ylabel('Voltage,mV'); 
    if i==1; title('Release'); end
    xlim([10 30]); ylim([-50 70])
    set(gca,'Fontsize',14)
end
