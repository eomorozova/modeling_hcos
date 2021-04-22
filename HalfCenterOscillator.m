

%{
              _       _   _
   __  _____ | | ___ | |_| |
   \ \/ / _ \| |/ _ \| __| |
    >  < (_) | | (_) | |_| |
   /_/\_\___/|_|\___/ \__|_|

## half-center oscillator


**Description**

makes a half-center oscillator model based on leech heartbeat interneurons
from Hill et al. 2001, Journal of Computational Neuroscience, with the addition of 
calcuim dynamics, calcium-dependent potassium current and inward-rectifier potassium current 

Synapse is a graded synapse as in Sharp et al. 1996

(In the Hill et al. 2001 model there are two types of inhibitory synaptic transmission
between the model elemental oscillator neurons:
graded transmission, which is dependent on the influx
of presynaptic Ca2+ through low-threshold Ca2+
channels (Angstadt and Calabrese, 1991); and spike mediated
transmission, which is dependent on influx
of presynaptic Ca2+ through high-threshold Ca2+ channels
during a spike (Lu et al., 1997). Spike-mediated
inhibition was modeled as a postsynaptic conductance
that is triggered by presynaptic spikes.)

%}

function x = HalfCenterOscillator()

%channels = {'hill/CaF','hill/CaS','hill/HCurrent','hill/KACurrent','prinz/KCa',...
%    'hill/KCurrent1','hill/KCurrent2','amarillo/Kir','Leak','hill/NaP','hill/NaV'};

channels = {'hill/CaF','hill/CaS','hill/HCurrent','hill/KACurrent','prinz/KCa',...
    'hill/KCurrent1','hill/KCurrent2','amarillo/Kir','Leak','hill/NaP','hill/NaV'};

x = xolotl;

A = 0.005; % area(mm^2) (assuming that the cell diameter is 40 uM (Tolbert LP, Calabrese RL (1985)))

x.add('compartment','cell1','A',A,'Cm',0.5/A) % right heart interneuron
x.add('compartment','cell2','A',A,'Cm',0.5/A) % left heart interneuron

% add calcium mechanism
f=20; % based on thin shell calculation (f=tau/(2*F*A*d)), d=0.1 um, tau=200 msec
x.cell1.add('prinz/CalciumMech','f',f);
x.cell2.add('prinz/CalciumMech','f',f);

for j = 1:length(channels)
    x.cell1.add(channels{j});
    x.cell2.add(channels{j});        
end
    
% configure gbars
gCaF=5; gCaS=3.2; gH=4; gKA=80; gK1=100; gK2=80; gKCa=100; gKir=100; gNaP=5; gNaV=200; gL=8; % nS
gbar=[gCaF, gCaS, gH, gKA, gKCa, gK1, gK2, gKir, gL, gNaP, gNaV];
gbar = gbar/(A*1000); % convert to uS/mm^2

x.cell1.set('*gbar', gbar);
x.cell2.set('*gbar', gbar);

% configure reversal potentials
ECaF=135; ECaS=135; EH=-21; EKA=-70; EKCa=-70; EK1=-70; EK2=-70; EKir=-70; ENaP=45; ENaV=45; EL=-60; % mV
E=[ECaF, ECaS, EH, EKA, EKCa, EK1, EK2, EKir, EL, ENaP, ENaV];

x.cell1.set('*E',E);
x.cell2.set('*E',E);


% add reciprocal inhibitory synapses

% Vth=-40; % synaptic threshold (mV) for release mechanism
Vth=-55; % synaptic threshold (mV) for escape mechanism

x.connect('cell1','cell2','Graded','gmax',5,'Vth',Vth,'Delta',2); 
x.connect('cell2','cell1','Graded','gmax',5,'Vth',Vth,'Delta',2); 

% set initial conditions to be slightly different for two cells
x.cell1.V = -50.1; 
x.cell2.V = -50; 







