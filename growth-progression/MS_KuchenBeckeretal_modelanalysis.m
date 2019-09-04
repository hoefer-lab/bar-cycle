% analysis and presentation of the growth-progression model results

% Erika Kuchen 2018


%% Figure 3 - Figure Supplement 3A
% best fit model parameters and confidence bounds
% loads in the simulation runs with all tested parameter sets

ABC_parameterBounds_weighted_alldatasets



%% Figure 3 - Figure Supplement 3B
% cell age and size distribution over time

% load in the population data frame generated in R by the growth-progression model.
% analyse cell age and cell size over time
% parameters below specifically for rep3.

% change the input of mat.
pops = double.empty;

for j=1:10  % load in several population data frames 
    load(['rep3exponential15gen',num2str(j),'reps_stationary.mat']);
    % when loading data frame from R: 
    pops = [pops; [pop.Tdiv,pop.absoluteTimeofDivision, pop.initialSize, pop.generation]]; % loading in the data frame from R.
    % when loading as array: 
    %pops = [pops;pop(:,[1,34,32,33])]; % only need information on cycle length (1), absolute division time (2), initial cell size (3) and generation (4)
end

% ADJUST PARAMETERS depending on the dataset and model growth parameters
k = 0.0405; % for rep3, % use the correct growth rate of the fitted dataset
growthTerm = 'exponential'; % 'logistic'; % growth model
snorm =1;
smax = 5;

fsize =8; % font size

% end of user input

pop = pops;
clear pops

alive = pop(:,2)-pop(:,1); % time when cell was born
gen = max(pop(:,4)); % maximum cell generation

tend = round(min(alive(pop(:,4)==gen)));
int = 0:(tend+3); % bins
intstats = NaN(length(int),8);

% calculate current cell age and size distribution over time in 1h bins
for l=1:length(int)
    use = logical((alive<=int(l)).*(pop(:,2)>=int(l))); % find cells alive during the time interval, i.e. cells already born and that have not yet divided
    allcells = int(l)-alive(use);  % determine the cell age at the current bin (time)
    switch growthTerm % calculate the current cell size, dependent on the growth model and growth rate
        case 'exponential'
            snow = pop(use,3).*exp(k*allcells);
        case 'logistic'
            snow = smax/(1+ ((smax-pop(use,3)) / pop(use,3))*exp(-k*allcells));
    end
    sage = sort(allcells);
    siage = length(allcells);
    sagex = sort(snow);
    intstats(l,:) = [mean(allcells),median(allcells),sage(max(1,round(siage*0.25))), sage(round(siage*0.75)), ...
        mean(snow),median(snow),sagex(max(1,round(siage*0.25))), sagex(round(siage*0.75))];
end

% display the current cell age and cell size (median and interquartile range)
figure
subplot(2,1,1) % cell age
sub = 0:tend;
subs = intstats(1:(tend+1),:);
order = sort(1:length(sub),'descend');
hp = patch([sub,sub(order)],[subs(:,3);subs(order,4)],'red'); hold on
set(hp,'FaceColor',[0.8,0.8,0.8],'EdgeColor','none')
hold on

plot(sub,intstats(1:tend+1,2),'k');
ylabel('age (h)','Fontsize',fsize);
box off
set(gca,'TickDir','out','Fontsize',fsize)
xlim([0 tend])

% size of cells
subplot(2,1,2)
hp = patch([sub,sub(order)],[subs(:,7);subs(order,8)],'red'); hold on
set(hp,'FaceColor',[0.8,0.8,0.8],'EdgeColor','none')
plot(sub,intstats(1:tend+1,6),'k');
ylabel('cell size (a.u.)','Fontsize',fsize);
xlabel('Time (h)','Fontsize',fsize)
box off
set(gca,'TickDir','out','Fontsize',fsize,'YTick',0.4:0.2:1.4,'XColor','k','YColor','k')
xlim([0 tend])
ylim([0.5 1.4])