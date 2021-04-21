x = xolotl.examples.networks.LeechHeartbeat;

A=1;

% add Ca mechanism
f=940; % 0.94uM/1000 (because all the conductances are 1000 times smaller than in stg... (asume small A?))

x.HN3R.add('prinz/CalciumMech','f',f);
x.HN3L.add('prinz/CalciumMech','f',f);

% add KCa
gKCa=0;
x.HN3R.add('prinz\KCa','gbar',gKCa/(A*1000),'E',-70);
x.HN3L.add('prinz\KCa','gbar',gKCa/(A*1000),'E',-70);

% add Kir
gKir=0;
x.HN3R.add('amarillo\Kir','gbar',gKir/(A*1000),'E',-70);
x.HN3L.add('amarillo\Kir','gbar',gKir/(A*1000),'E',-70);

Vth = -50; % synaptic threshold (mV)
x.set('*.Vth',Vth);

x.t_end = 50e3; % msec
x.output_type = 1;

x.closed_loop = false;
data = x.integrate;

% set conductances to the values similar to other K currents in the model
gKCa=80/1000; % uS/mm^2
x.set('HN3R.KCa.gbar',gKCa); x.set('HN3L.KCa.gbar',gKCa);

gKir=80/1000;
x.set('HN3R.Kir.gbar',gKir); x.set('HN3L.Kir.gbar',gKir);

x.closed_loop = false;
data1 = x.integrate;

%figure('outerposition',[300 300 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
clf
subplot(3,1,1), 
ax1=plot(x.time,data.HN3L.V,'k','linewidth',1)
hold on, plot(x.time,data1.HN3L.V,'b','linewidth',1)
yyaxis right
plot(x.time,data1.HN3L.Ca,'r','linewidth',2);
ylabel('Ca,uM')
yyaxis left
line([0 x.t_end/10^3],[Vth Vth],'color','k','linestyle','--')
ylabel('V_m (mV)');
ylim([-70 5]);

subplot(3,1,2), plot(x.time,data.HN3R.V,'k','linewidth',1)
hold on, plot(x.time,data1.HN3R.V,'b','linewidth',1)
line([0 x.t_end/10^3],[Vth Vth],'color','k','linestyle','--')
ylabel('V_m (mV)')
ylim([-70 5]);
legend('control','+KCa +Kir','location','best')

subplot(3,1,3), plot(x.time,data1.HN3R.KCa.I,'linewidth',1.5)
hold on, plot(x.time,data1.HN3R.KACurrent.I,'linewidth',1.5)
hold on, plot(x.time,data1.HN3R.Kir.I,'linewidth',1.5)
ylabel('I,nA'); xlabel('Time (s)');
legend('KCa','KA','Kir','location','best')

figlib.pretty()