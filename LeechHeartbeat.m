

%{
              _       _   _
   __  _____ | | ___ | |_| |
   \ \/ / _ \| |/ _ \| __| |
    >  < (_) | | (_) | |_| |
   /_/\_\___/|_|\___/ \__|_|

## half-center oscillator

**Syntax**

```matlab
x = xolotl.examples.networks.LeechHeartbeat();
```

**Description**

makes a half-center oscillator model based on leech heartbeat interneurons
from Hill et al. 2001, Journal of Computational Neuroscience. 

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

function x = LeechHeartbeat()

channels = {'hill/CaF','hill/CaS','hill/HCurrent','hill/KACurrent','hill/KCurrent1','hill/KCurrent2','Leak','hill/NaP','hill/NaV'};

x = xolotl;

A = 0.005; % mm^2 (assuming that the cell diameter is 40 uM (Tolbert LP, Calabrese RL (1985)))
x.add('compartment','HN3R','A',A,'Cm',0.5/A) % right heart interneuron
x.add('compartment','HN3L','A',A,'Cm',0.5/A) % left heart interneuron

for j = 1:length(channels)
    x.HN3R.add(channels{j});
    x.HN3L.add(channels{j});        
end
    
% configure gbars
gCaF=5; gCaS=3.2; gH=4; gKA=80; gK1=100; gK2=80; gNaP=5; gNaV=200; gL=8; % nS
gbar=[gCaF, gCaS, gH, gKA, gK1, gK2, gL, gNaP, gNaV];
gbar = gbar/(A*1000); % convert into uS/mm^2
%gCaF=5; gCaS=3.2; gH=4; gKA=100; gK1=80; gK2=200; gNaP=5; gNaV=200; gL=8/(A*1000); % nS
x.HN3R.set('*gbar', gbar);
x.HN3L.set('*gbar', gbar);

% configure reversal potentials
ECaF=135; ECaS=135; EH=-21; EKA=-70; EK1=-70; EK2=-70; ENaP=45; ENaV=45; EL=-60; % mV
E=[ECaF, ECaS, EH, EKA, EK1, EK2, EL, ENaP, ENaV];
x.HN3R.set('*E',E);
x.HN3L.set('*E',E);


% add reciprocal inhibitory synapses

% Vth=-40; % synaptic threshold (mV) for release mechanism
Vth=-55; % synaptic threshold (mV) for escape mechanism


x.connect('HN3R','HN3L','Graded','gmax',5,'Vth',Vth,'Delta',2); 
x.connect('HN3L','HN3R','Graded','gmax',5,'Vth',Vth,'Delta',2); 


x.HN3R.V = -50.1; 
x.HN3L.V = -50; % set initial conditions to be slightly different for two cells







