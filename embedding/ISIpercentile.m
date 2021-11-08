
% calculate ISIs and percentiles of ISIs

function [isiprtile, ISI] = ISIpercentile(spiketimes)

for i=1:size(spiketimes,1)
    for j=1:size(spiketimes,2)
        
        ISI{i,j}=diff(spiketimes{i,j});
        
        isiprtile{i}(j,:)=prctile(ISI{i,j},1:10:100);
        
    end
    
    isiprtile{i}(isnan(isiprtile{i}))=0;
end

