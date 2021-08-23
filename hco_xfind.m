% use xfind to sample random neual models and log them to disk

% Searching for half-center oscillators with
% Escape mechanism (Vth=-55 mV), one value of synaptic conductance gsyn=5nS
% and neurons with conducnatces randomly sampled from  a uniform distribution

% set up a neuron/ network
p = xfind;
p.x = HalfCenterOscillator;


% define parameters to randomly sample and bounds

p.ParameterNames = p.x.find('*gbar');

A=0.005; % area of the neuron (mm^2)

%         CaF  CaS     H     KA      KCa    K1      K2       Kir  Leak   NaP   NaV     
% neuron #1
Upper1 = [10,   6.4,   8,   200,   200,   200,    200,     200,   16,   15,   400]./(A*1000);
% neuron #2
Upper2 = [10,   6.4,   8,   200,   200,   200,    200,     200,   16,   15,   400]./(A*1000);

p.Upper=[Upper1, Upper2];
p.Lower = zeros(22,1);


% defining a smaple function: uniform random sampler

p.SampleFcn = @p.uniformRandom;


% defininf the simulation function: measureFiringRate
% method that calculates mean firind rate of a network

p.SimFcn = @p.measureFiringRate;


% we only want to keep non silent networks
p.DiscardFcn = @(data) data <= 0;

p.simulate
