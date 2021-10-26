%% get the data

file = dir('C:\Users\moroz\Desktop\HCO_git\HCO-modeling\hco_modeling\c4d20*.xfind')
data=[];
for i=1:length(file)
fileid=fopen(['C:\Users\moroz\Desktop\HCO_git\HCO-modeling\hco_modeling\',file(i).name]);
data =[data; fread(fileid,'double')];
end

%% get parameters

st1=[]; st2=[]; fr1=[]; params=[];
istart1=[]; istart2=[]; iend1=[]; iend2=[];
istart1 = find(data==-1); iend1 = find(data==-2);
istart2 = find(data==-3); iend2 = find(data==-4);

params(1,:) = data(1:istart1(1)-3);
fr1(1) = data(istart1(1)-1);
fr2(1) = data(istart1(1)-2);
for i=1:length(istart1)
   try; params(i+1,:) = data(iend2(i)+1:istart1(i+1)-3); end;
   try; fr1(i+1) = data(istart1(i+1)-1); end;
   try; fr2(i+1) = data(istart1(i+1)-2); end;
   st1{i} =  data(istart1(i)+1:iend1(i)-1)/10^4;
   st2{i} =  data(istart2(i)+1:iend2(i)-1)/10^4;
end
%% reintegrate with the found parameters and plot voltage traces + IMI

%V=[]; stmi1=[]; stmi2=[]; 

ie=4747; % escape only
ir1=17714; ir2=21512;
gMI = 10/(A*1000);

%x = HalfCenterOscillator; % control
x1 = HalfCenterOscillator_IMI; % with IMI

%x.t_end = 30e3;
x1.t_end = 35e3;

for i = ir1:ir2
    if mod(i,100)==0
        i
    end
    %x.set('*gbar',params(i,:));
    x1.set('*gbar',[params(i,1:9),gMI,params(i,10:20),gMI,params(i,21:end)]);

    %p.x.I_ext=0;
    %V{i} = x.integrate;
    VMI{i} = x1.integrate;
    
    stmi1_1 = analysis.crossings(VMI{i}(5e3/x1.dt:end,1)',-30); % remove transient
    stmi2_1 = analysis.crossings(VMI{i}(5e3/x1.dt:end,2)',-30);
    
    
    stmi1{i} = stmi1_1{1}./10^4;  stmi2{i} = stmi2_1{1}./10^4;

    VMI{i}=[];
end


%% plot traces

i=19;
gMI = 10/(A*1000);

x = HalfCenterOscillator; % control
x1 = HalfCenterOscillator_IMI; % with IMI

x.t_end = 30e3; x1.t_end = 30e3;

x.set('*gbar',params(i,:));
x1.set('*gbar',[params(i,1:9),gMI,params(i,10:20),gMI,params(i,21:end)]);

%p.x.I_ext=0;
V{i} = x.integrate;
VMI{i} = x1.integrate;
    
figure(2); clf
time = (1:length(V{i}))*1e-3*x.dt;
subplot(2,1,1), plot(time,V{i}(:,1),'k')
hold on, plot(time,VMI{i}(:,1),'r')
subplot(2,1,2), plot(time,V{i}(:,2),'k')
hold on, plot(time,VMI{i}(:,2),'r')
figlib.pretty('LineWidth', 1, 'PlotLineWidth', 1, 'PlotBuffer', 0)
% pause

%% plot spike raster in control and with Imi

i=19;
figure(1),clf
%analysis.plot_raster([st1(1:50),st2(1:50)],1,st1{1}(end),'k',1)
%analysis.plot_raster([st1(i),st2(i)],1,st1{1}(end),'k',1)
analysis.plot_raster([st1(i),st2(i),stmi1(i),stmi2(i)],1,st1{1}(end),'k',1)
%hold on, display.plothorzline(50,'r') 
xlabel('Time, sec'); %xlim([1 10])
set(gca,'Fontsize',16)

%% make it so that the first neurons has always more spikes
% because it doesn't matter which neuron is which

st_1 = st1; st_2 = st2;

for i=1:length(st1)
    if length(st1{i})<length(st2{i})
        st_1(i)=st2(i);
        st_2(i)=st1(i);
    end
end

st1=st_1; st2=st_2;


stmi_1 = stmi1; stmi_2 = stmi2;

for i=1:length(stmi1)
    if length(stmi1{i})<length(stmi2{i})
        stmi_1(i)=stmi2(i);
        stmi_2(i)=stmi1(i);
    end
end

stmi1=stmi_1; stmi2=stmi_2;
%% calculate ISIs, phases and percentiles of ISIs and phases

%st=[st1(1:4747),stmi1;st2(1:4747),stmi2];
st=[st1(1:24066),stmi1;st2(1:24066),stmi2];

[isiprtile, ISIs] = ISIpercentile(st);

phasesprtile = PhasesPrtiles(st);

allprtile=[isiprtile{1},isiprtile{2},phasesprtile(:,1:10),phasesprtile(:,11:20)];

%% with firing rates
% firing rate
for i=1:length(st)
    FR1(i) = mean(diff(st{1,i}));
    FR2(i) = mean(diff(st{2,i}));
    if ~isempty(ISIs{1,i})==1 && length(ISIs{1,i})>1
        smax1(i) = max(diff(sort(ISIs{1,i})));
    else
        smax1(i) = 0;
    end
    if ~isempty(ISIs{2,i})==1 && length(ISIs{2,i})>1
        smax2(i) = max(diff(sort(ISIs{2,i})));
    else
        smax2(i) = 0;
    end
end

FR1(isnan(FR1))=0; FR2(isnan(FR2))=0;

allprtile1=[isiprtile{1},isiprtile{2},phasesprtile(:,1:10),phasesprtile(:,11:20),FR1',FR2',smax1',smax2'];

%%
i=37;
%clf, subplot(2,2,1), bar([isiprtile{1}(2,:),isiprtile{2}(2,:), phasesprtile(2,1:10), phasesprtile(2,11:20)])
clf, subplot(2,2,1), bar(allprtile(i,:))
subplot(2,2,2), bar(zscore(allprtile(i,:)))
subplot(2,2,3), analysis.plot_raster(st(:,i),1,st{1}(end),'k',1)

subplot(2,2,4),plot((diff(sort(ISIs{1,i}))))

%% do tsne

%addpath(genpath('C:\Users\moroz\Desktop\HCO_git'))

opts.perplexity=30;
%opts.late_exag_coeff = 2.5;
%opts.start_late_exag_iter = 800;

%[mappedX, costs] = fast_tsne(allprtile, opts);
%[mappedZ, costs] = fast_tsne((zscore(allprtile)), opts);
[mappedFR, costs] = fast_tsne(((allprtile1)), opts);
%[mappedFRZ, costs] = fast_tsne((zscore(allprtile1)), opts);

%% plot tsne map and spike rasters

i_i= [5,3,35,49,67,102,9,55,23,69,93,89,108,129,126,107,275,329,2,537,560];

colors = colormaps.dcol(length(i_i));

clf, subplot(1,2,1), scatter(mappedX(:,1),mappedX(:,2),5,'k','filled'), hold on,

for i=1:length(i_i);
    scatter(mappedX(i_i(i),1),mappedX(i_i(i),2),50,colors(i,:),'filled');
end

ylabel('t-SNE 2'); xlabel('t-SNE 1')
title('Based on ISIs and phases')
axis square

% plot rasters
for i=1:length(i_i)
    subplot(length(i_i),4,[i*4-1,i*4]); analysis.plot_raster([st1(i_i(i));st2(i_i(i))],1,20,colors(i,:),1);
    xlim([1 20]); axis off
    display.plothorzline(0)
end

%% colorcode based on control vs IMI
%ictrl = 4747;
ictrl = length(st)/2;

clf, subplot(1,2,1), scatter(mappedX(:,1),mappedX(:,2),5,'k','filled'), hold on,
h1=scatter(mappedX(1:ictrl,1),mappedX(1:ictrl,2),5,'k','filled'); hold on,
h2=scatter(mappedX(ictrl:end,1),mappedX(ictrl:end,2),5,'r','filled');
legend([h1,h2],'control','+IMI')

%% colorcode based on the firing patten

%% find bursts

%[hcostat] = analysis.hco_stat(data.V{i}{j}{jj}(ii,:),data.Fs{i});

for i = 1:length(st)
    for n = 1:2 % neuron
        %burststat{i,n}=analysis.detect_bursts(st{n,i},mean(diff(st{n,i})),mean(diff(st{n,i}))+0.5,1);
        burststat{i,n}=analysis.detect_bursts(st{n,i},mean(diff(st{n,i})),mean(diff(st{n,i}))+0.01,1);
    end
end

%% classify the state (silent, spiking, HCO)
addpath('C:\Users\moroz\Documents\Katya_DynamicClamp\hco_related_codes')

networkstate=[]; state=[]; percentSingleSpikeBursts=[];

for j = 1:length(st)
    
    %if ~isempty(st(:,j))==1
        
        percentSingleSpikeBursts(j) = calcPercentSingleSpikeBursts(st{1,j}',st{2,j}');
        
        if length(st{1,j})<=5 && length(st{2,j})<=5 % less than 5 spike in 20 sec in both neurons
            
            networkstate(j) = 1; % silent
            
        elseif (length(st{1,j})<5  && length(st{2,j})>5)...
                || (length(st{1,j})>5  && length(st{2,j})<5) % asymmetric
            
            if (length(burststat{j,1}.N)>1 || length(burststat{j,2}.N)>1)
                
                networkstate(j) = 7; % asymmetric bursting     
            else
                networkstate(j) = 2; % asymmetric (only one neuron spikes)  
            end
            
        elseif (length(st{1,j})>5  && length(st{2,j})>5 ...
                && length(burststat{j,1}.N)==1 && length(burststat{j,2}.N)==1) % spiking
            
            if percentSingleSpikeBursts(j)>0.8 % more than 80% alternating
                networkstate(j) = 4; % antiphase spiking
            else
                networkstate(j) = 3; % irregular spiking
            end
         
        elseif (length(st{1,j})>5  && length(st{2,j})>5 && ...
                length(burststat{j,1}.N)>1 && length(burststat{j,2}.N)>1) % HCO
            
            networkstate(j) = 5; % HCO
            
      %  elseif (length(hcostat{j}{k}{1}{i}.st(hcostat{j}{k}{1}{i}.st<60*Fs(j)))>5  && length(hcostat{j}{k}{2}{i}.st(hcostat{j}{k}{2}{i}.st<60*Fs(j)))>5 && ...
      %          hcostat{j}{k}{1}{i}.T1_mean>0 && hcostat{j}{k}{2}{i}.T1_mean>0 && hcostat{j}{k}{1}{i}.nSpks_mean==0 && hcostat{j}{k}{2}{i}.nSpks_mean==0) % single spike HCO
            
      %      networkstate{j}{k}(i) = 4; % antiphase spiking
            
        else
            
            networkstate(j) = 6; % something is not right
            
        end
    %end
end

%%
i=iasym(8);
n=1;
burststat{i,n}=analysis.detect_bursts(st{n,i},mean(diff(st{n,i})),mean(diff(st{n,i}))+0.1,1);

burststat{i,1}
%clf,analysis.plot_raster(st(:,i),1,st{1}(end),'k',1)


%% subdivide assymatric spiking patterns into tonic and bursting

%% colorcode t-SNE map based on the activity pattern
isilent = find(networkstate==1);
iasym = find(networkstate==2);
iasymburst = find(networkstate==7);
iirreg = find(networkstate==3);
iantiphase = find(networkstate==4);
ihco = find(networkstate==5);
iother = find(networkstate==6);

colors = display.linspecer(9);

clf, subplot(3,2,1)
%scatter(mappedX(:,1),mappedX(:,2),5,'k','filled'); hold on,
hold on, scatter(mappedX(isilent,1),mappedX(isilent,2),5,'k','filled')
hold on, scatter(mappedX(iasym,1),mappedX(iasym,2),5,colors(3,:),'filled')
hold on, scatter(mappedX(iasymburst,1),mappedX(iasymburst,2),5,colors(2,:),'filled')
hold on, scatter(mappedX(iirreg,1),mappedX(iirreg,2),5,colors(7,:),'filled')
hold on, scatter(mappedX(iantiphase,1),mappedX(iantiphase,2),5,colors(4,:),'filled')
hold on, scatter(mappedX(ihco,1),mappedX(ihco,2),5,colors(1,:),'filled')
set(gca,'xtick',[]); set(gca,'ytick',[]); axis square
ylabel('t-SNE 2'); xlabel('t-SNE 1')
%legend('silent','asymmetric (tonic)', 'asymmetric (bursting)','irregular','antiphase','half-center')

for i=1:length(st)
    CV_isi(i) = std(diff(st{1,i}))./mean(diff(st{1,i}));
end

CV_isi(find(CV_isi>5))=0;

[~,icv] = sort(CV_isi);
mappedX_sorted = mappedX(icv,:);
subplot(3,2,3)
col=display.linspecer(numel(icv));
scatter(mappedX_sorted(:,1),mappedX_sorted(:,2),5,1:length(icv),'filled'), hold on
ylabel('t-SNE 2'); xlabel('t-SNE 1')
axis square

axes('Position',[.3 .6 .1 .05])
[f]=ksdensity(CV_isi(CV_isi>0.05)),
patch(linspace(1,max(CV_isi),100),f,1:100)
xlabel('CV')

% firing rate
for i=1:length(st)
    FR(i) = mean(diff(st{1,i}));
end
T = 1./FR;
[~,ifr] = sort(T);
mappedX_sortedT = mappedX(ifr,:);

%%
clf
subplot(1,2,1)
scatter(mappedX_sortedT(:,1),mappedX_sortedT(:,2),5,1:length(ifr),'filled'), hold on
ylabel('t-SNE 2'); xlabel('t-SNE 1')
axis square


%i_i= [5,3,35,49,67,102,9,55,23,69,93,89,108,129,126,107,275,329,2,537,560];
%i_i = [18,22, 49,102, 89,108, 23, 168, 55,93, 129,126,107]
i_i = [18,22, 49,102, 89,108, 23, 168, 55,93, 129,126,107];
colors = colormaps.dcol(length(i_i));

for i=1:length(i_i);
    scatter(mappedX(i_i(i),1),mappedX(i_i(i),2),50,colors(i,:),'filled');
end

%axes('Position',[.3 .3 .1 .05])
%[f]=ksdensity(T),
%patch(linspace(1,max(T),100),f,1:100)
%xlabel('Period, sec')


% plot rasters
for i=1:length(i_i)
    subplot(length(i_i),4,[i*4-1,i*4]); analysis.plot_raster([st1(i_i(i));st2(i_i(i))],1,20,colors(i,:),1);
    xlim([1 20]); axis off
    display.plothorzline(0)
end


%% plot rasters for bursting patterns only
clf
subplot(2,2,1)
scatter(mappedX_sortedT(:,1),mappedX_sortedT(:,2),5,1:length(ifr),'filled'), hold on
ylabel('t-SNE 2'); xlabel('t-SNE 1')
axis square


%i_i= [5,3,35,49,67,102,9,55,23,69,93,89,108,129,126,107,275,329,2,537,560];

iirreg1 = [29,126,107,52,64,107,362,290,452];
iweird = [54,56,137,194];

%i_i = [205,117,90,11,19,25, 77,537, 290,452, 2,560]; %123 % half-centers

%i_i = [34,74,46,126,64,81,279,52,117,233,19,111, 2,77,537]; % irregular

i_i = [34,74,46,126,64,81,279,52,117,233,19,111, 2,77,537]; % irregular/hco

%colors = colormaps.dcol(length(i_i));
colors = display.linspecer(length(i_i));

for i=1:length(i_i);
    scatter(mappedX(i_i(i),1),mappedX(i_i(i),2),50,colors(i,:),'filled');
end

%axes('Position',[.3 .3 .1 .05])
%[f]=ksdensity(T),
%patch(linspace(1,max(T),100),f,1:100)
%xlabel('Period, sec')

%% --------------------------------------------------------------------
% colorcode t-SNE map based on the activity pattern
isilent = find(networkstate==1);
iasym = find(networkstate==2);
iasymburst = find(networkstate==7);
iirreg = find(networkstate==3);
iantiphase = find(networkstate==4);
ihco = find(networkstate==5);
iother = find(networkstate==6);

colors = display.linspecer(9);

clf, g=0.07; l=0.07;
subtightplot(2,2,1,g,l,l)
%scatter(mappedX(:,1),mappedX(:,2),5,'k','filled'); hold on,
hold on, scatter(mappedFR(isilent,1),mappedFR(isilent,2),5,'k','filled')
hold on, scatter(mappedFR(iasym,1),mappedFR(iasym,2),5,colors(3,:),'filled')
hold on, scatter(mappedFR(iasymburst,1),mappedFR(iasymburst,2),5,colors(2,:),'filled')
hold on, scatter(mappedFR(iirreg,1),mappedFR(iirreg,2),5,colors(7,:),'filled')
hold on, scatter(mappedFR(iantiphase,1),mappedFR(iantiphase,2),5,colors(4,:),'filled')
hold on, scatter(mappedFR(ihco,1),mappedFR(ihco,2),5,colors(1,:),'filled')
set(gca,'xtick',[]); set(gca,'ytick',[]); axis square
ylabel('t-SNE 2'); xlabel('t-SNE 1')
%legend('silent','asymmetric (tonic)', 'asymmetric (bursting)','irregular','antiphase','half-center')


%% colorcode based on frequency
clf
%i_i= [5,3,35,49,67,102,9,55,23,69,93,89,108,129,126,107,275,329,2,537,560];

%ifind = find(mappedFR(:,1)>-20 & mappedFR(:,1)<0 & mappedFR(:,2)>-30 & mappedFR(:,2)<-10)

%i_i = [164,5,3,208] % tonic (asymmetric)
%i_i = [43,67,49,11,78,23,121,192, 85,6,29, 21,98]; % one cell is bursting

i_i = [34,46,126,81,117,233,19,111, 77,537,2, 279]; % irregular/hco

%i_i = [164,5,3,208, 43,67,49,11,78,23,121,192, 85,6,29, 21,98, 34,46,126,81,117,233,19,111, 77,537,2, 279]; % all

%colors = colormaps.dcol(length(i_i));
colors = display.linspecer(length(i_i));

for i=1:length(st)
    FR(i) = mean(diff(st{1,i}));
end
T = 1./FR;
[~,ifr] = sort(T);
mappedFR_sortedT = mappedFR(ifr,:);

colormap(colormaps.magma)

g=0.07; l=0.07;
%subtightplot(2,2,3,g,l,l)
subtightplot(1,2,1,g,l,l)
%scatter(mappedZ(:,1),mappedZ(:,2),5,'k','filled'), hold on
%scatter(mappedFR(:,1),mappedFR(:,2),5,'k','filled'), hold on
scatter(mappedFR_sortedT(:,1),mappedFR_sortedT(:,2),5,1:length(ifr),'filled'), hold on
%scatter(mappedFRZ(:,1),mappedFRZ(:,2),5,'k','filled'), hold on
ylabel('t-SNE 2'); xlabel('t-SNE 1')
axis square
for i=1:length(i_i);
    %scatter(mappedZ(i_i(i),1),mappedZ(i_i(i),2),50,colors(i,:),'filled');
    scatter(mappedFR(i_i(i),1),mappedFR(i_i(i),2),100,colors(i,:),'filled','MarkerEdgeColor',[0,0,0]); 
    %scatter(mappedFR_sortedT(i_i(i),1),mappedFR_sortedT(i_i(i),2),50,colors(i,:),'filled');     
    %scatter(mappedFRZ(i_i(i),1),mappedFRZ(i_i(i),2),50,colors(i,:),'filled');      
end



%ifind = find(mappedX(:,1)>34 & mappedX(:,1)<44 & mappedX(:,2)>8 & mappedX(:,2)<16)
%ifind = find(mappedX(:,1)>25 & mappedX(:,1)<45 & mappedX(:,2)>-6 & mappedX(:,2)<8)

% plot rasters
for i=1:length(i_i)
    subplot(length(i_i),4,[i*4-1,i*4]); analysis.plot_raster([st1(i_i(i));st2(i_i(i))],1,20,colors(i,:),1);
    xlim([1 20]); axis off
    display.plothorzline(0);
end


%% colorcode based on ctrl vs Imi

isilent_ctrl = isilent(isilent<=ictrl);
iasym_ctrl = iasym(iasym<=ictrl);
ihco_ctrl = ihco(ihco<=ictrl);
iother_ctrl = iother(iother<=ictrl);
izero_ctrl = izero(izero<=ictrl);

isilent_mi = isilent(isilent>ictrl);
iasym_mi = iasym(iasym>ictrl);
ihco_mi = ihco(ihco>ictrl);
iother_mi = iother(iother>ictrl);
izero_mi = izero(izero>ictrl);

states = [length(isilent)/length(st), length(iasym)/length(st),...
    length(ihco)/length(st),length(iother)/length(st), length(izero)/length(st)];

states_ctrl = [length(isilent_ctrl)/ictrl, length(iasym_ctrl)/ictrl,...
    length(ihco_ctrl)/ictrl,length(iother_ctrl)/ictrl, length(izero_ctrl)/ictrl];

states_mi = [length(isilent_mi)/ictrl, length(iasym_mi)/ictrl,...
    length(ihco_mi)/ictrl,length(iother_mi)/ictrl, length(izero_mi)/ictrl];


colors = display.linspecer(numel(states));

clf, subplot(2,2,1), scatter(mappedX(:,1),mappedX(:,2),5,'k','filled');

h1=scatter(mappedX(1:ictrl,1),mappedX(1:ictrl,2),5,'k','filled'); hold on,
h2=scatter(mappedX(ictrl:end,1),mappedX(ictrl:end,2),5,'r','filled');
legend([h1,h2],'control','+IMI')

%hold on, scatter(mappedX(iasym,1),mappedX(iasym,2),5,colors(2,:),'filled')
%hold on, scatter(mappedX(ihco,1),mappedX(ihco,2),5,colors(3,:),'filled')

set(gca,'xtick',[]); set(gca,'ytick',[]); axis square
ylabel('t-SNE 2'); xlabel('t-SNE 1')

% a mondrian plot (a treemap)
subplot(2,4,5), mondrian(states_ctrl, [0,1,2,3,4])
axis square
title('Control')

subplot(2,4,6), mondrian(states_mi, [0,1,2,3,4])
axis square
title('+IMI')
