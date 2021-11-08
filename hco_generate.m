% use xfind to sample random neual models and log them to disk

% Searching for half-center oscillators with
% Escape mechanism (Vth=-55 mV),
% and neurons with conducnatces randomly sampled from  a uniform distribution

% set up a neuron/ network
p = xfind_hco;
p.x = HalfCenterOscillator;
p.x.t_end = 90e3; % msec

% define parameters to randomly sample and bounds

p.ParameterNames = [p.x.find('*gbar'); p.x.find('*gmax')];

A=0.005; % area of the neuron (mm^2)

%         CaF  CaS     H     KA      KCa    K1      K2       Kir  Leak   NaP   NaV  Synapse   
% neuron #1
Upper = [10,   6.4,   8,   200,   200,   200,    200,     200,   16,   15,   400]./(A*1000);


p.Upper=Upper;
p.Lower = zeros(11,1);


% defining a smaple function: uniform random sampler

p.SampleFcn = @p.uniformRandom;


% defininf the simulation function: measureFiringRate
% method that calculates mean firind rate of a network

p.SimFcn = @p.measureFiringRate;


% keep all the networks
p.DiscardFcn = @(data) data <= 0;

p.simulate

%p.parallelSearch

% wait for a bit
%pause(1200)

%cancel(p.workers)
