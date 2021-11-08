
classdef xfind_hco < handle



properties (Access = private)

end % private props

properties

	x (1,1) xolotl

	% names of parameters
	ParameterNames (:,1) cell

	% bounds
	Upper (:,1) double
	Lower (:,1) double 


	% stores futures of executing parallel tasks
	workers (:,1) parallel.FevalFuture	

	% allow user-defined custom sample function
	SampleFcn function_handle

	% allow user to run some function on the model
	% after we sample it
	% should return a structure data with data
	SimFcn function_handle 

	% this function determines if the sampled point
	% is logged or not. 
	% data is the data returned by SimFcn
	DiscardFcn (1,1) function_handle = @(data) false

end % props



methods 

	function self = xfind_hco()
		self.SampleFcn = @self.uniformRandom;
		self.SimFcn = @self.measureFiringRate;
	end % constructor



	% helper method that writes to disk
	function writeToDisk(self,f, data)
		params = self.x.get(self.ParameterNames);
        V = self.x.integrate;
        
        %clf, subplot(2,1,1),plot(V(:,1)); ylim([-65 0])
        %subplot(2,1,2),plot(V(:,2)); ylim([-65 0])
        
        %self.x.get('*gmax')
        
        st1 = veclib.computeOnsOffs(V(:,1) > -30);
        st2 = veclib.computeOnsOffs(V(:,2) > -30);
        %length(st1)
        %length(st2)
		% stream data to disk
        fwrite(f,params,'double');
		fwrite(f,data,'double');
        fwrite(f,[-1,st1',-2],'double'); % save spike times for neuron 1
        fwrite(f,[-3,st1',-4],'double'); % save spike times for neuron 1
	end


	% returns a hash of the neuron model
	% that is used to save data
	function H = hash(self)
		H = self.x.hash(1:5);
	end


	% picks a point from a uniform random distribution
	% within bounds
	function uniformRandom(self)
		params1 = (self.Upper-self.Lower).*rand(length(self.Upper),1) + self.Lower;
        
        Uppersyn = 20; % 10;
        psyn1 = (Uppersyn).*rand(length(Uppersyn),1);
        %psyn2 = (Uppersyn).*rand(length(Uppersyn),1);
        
        % to have close to identical neurons
        params = [params1;params1+params1.*(rand(length(self.Upper),1)-0.5);...
            psyn1; psyn1+psyn1.*(rand(length(Uppersyn),1)-0.5)];
       % params = [params1;params1]; % to have identical neurons
        self.x.set(self.ParameterNames,params);
	end % uniformRandom


	% this example function simulates the model and
	% measures the firing rate in the first compartment
	function data = measureFiringRate(self)
		V = self.x.integrate;
		fr1 = xtools.findNSpikes(V(5e3/self.x.dt:end,1))/(self.x.t_end*1e-3-5);
        fr2 = xtools.findNSpikes(V(5e3/self.x.dt:end,2))/(self.x.t_end*1e-3-5);
        
        data = [fr1,fr2];
	end





end % methods



methods (Static)



end % static methods






end % classdef

