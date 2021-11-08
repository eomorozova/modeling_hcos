
%% calculate phases and percentile of the phases

function [phasesprtile,phases1,phases2] = PhasesPrtiles(spiketimes)

st1 = spiketimes(1,:);
st2 = spiketimes(2,:);

prevspk=[]; nextspk=[];

% phase for the 1st neuron
for i = 1:length(st1)
    if ~isempty(st1{i}) && ~isempty(st2{i})...
            && max(st2{i})>min(st1{i}) && max(st1{i})>min(st2{i})...
            && length(st1{i})>10 && length(st2{i})>10
        
        for j = 1:length(st1{i})
            prev_other_spike = st2{i}(find(st2{i}<st1{i}(j),1,'last'));
            
            if isempty(prev_other_spike)
                continue
            end
            
            next_other_spike = st2{i}(find(st2{i}>st1{i}(j),1,'first'));
            
            if isempty(next_other_spike)
                continue
            end
            
            prevspk(j) = prev_other_spike;
            nextspk(j) = next_other_spike;
            
            phases1{i}(j) = (st1{i}(j)-prev_other_spike)/(next_other_spike-prev_other_spike);
            
        end
    else
        phases1{i}=[];
    end
    
    if isempty(prevspk)==1 && isempty(nextspk)==1
        phases1{i}=[];
    end
    prevspk=[]; nextspk=[];
    
    if  ~isempty(phases1{i})
        phasesprtile1(i,:) = prctile(phases1{i},1:10:100);
    else
        phasesprtile1(i,:)=zeros(1,10);
    end
end

% phase for the 2nd neuron
for i = 1:length(st2)
    
    if ~isempty(st1{i}) && ~isempty(st2{i})...
            && max(st1{i})>min(st2{i}) && max(st2{i})>min(st1{i})...
            && length(st1{i})>10 && length(st2{i})>10   
        
        for j = 1:length(st2{i})
            prev_other_spike = st1{i}(find(st1{i}<st2{i}(j),1,'last'));
            
            if isempty(prev_other_spike)
                continue
            end
            
            next_other_spike = st1{i}(find(st1{i}>st2{i}(j),1,'first'));
            
            if isempty(next_other_spike)
                continue
            end
            
            prevspk(j) = prev_other_spike;
            nextspk(j) = next_other_spike;
            
            phases2{i}(j) = (st2{i}(j)-prev_other_spike)/(next_other_spike-prev_other_spike);
            
        end
    else
        phases2{i}=[];
    end
    
    if isempty(prevspk)==1 && isempty(nextspk)==1
        phases2{i}=[];
    end
    prevspk=[]; nextspk=[];
    
    if  ~isempty(phases2{i})
        phasesprtile2(i,:) = prctile(phases2{i},1:10:100);
    else
        phasesprtile2(i,:)=zeros(1,10);
    end
end

phasesprtile = [phasesprtile1, phasesprtile2];