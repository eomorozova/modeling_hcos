x = HalfCenterOscillator;

A=0.005; % mm^2

Vth = -50; % synaptic threshold (mV)
x.set('*.Vth',Vth);

x.t_end = 50e3; % msec
x.output_type = 1;

x.closed_loop = false;
data = x.integrate;

% set gKCa and gKir to 0 for comparison 
gKCa=0/(A*1000); % uS/mm^2
x.set('cell1.KCa.gbar',gKCa); x.set('cell2.KCa.gbar',gKCa);

gKir=0/(A*1000);
x.set('cell1.Kir.gbar',gKir); x.set('cell2.Kir.gbar',gKir);

x.closed_loop = false;
data1 = x.integrate;

%figure('outerposition',[300 300 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
clf
subplot(3,1,1), 
plot(x.time,data1.cell1.V,'k','linewidth',1)
hold on, plot(x.time,data.cell1.V,'b','linewidth',1)
yyaxis right
plot(x.time,data.cell1.Ca,'r','linewidth',2);
ylabel('Ca,uM')
yyaxis left
line([0 x.t_end/10^3],[Vth Vth],'color','k','linestyle','--')
ylabel('V_m (mV)');
ylim([-70 5]);

subplot(3,1,2), 
plot(x.time,data1.cell2.V,'k','linewidth',1)
hold on, plot(x.time,data.cell2.V,'b','linewidth',1)
line([0 x.t_end/10^3],[Vth Vth],'color','k','linestyle','--')
ylabel('V_m (mV)')
ylim([-70 5]);
legend('canonical','+KCa +Kir','location','best')

subplot(3,1,3), plot(x.time,data.cell1.KCa.I,'linewidth',1.5)
hold on, plot(x.time,data.cell1.KACurrent.I,'linewidth',1.5)
hold on, plot(x.time,data.cell1.Kir.I,'linewidth',1.5)
ylabel('I,nA'); xlabel('Time (s)');
legend('KCa','KA','Kir','location','best')

figlib.pretty()