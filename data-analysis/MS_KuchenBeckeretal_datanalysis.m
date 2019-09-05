%% Data analysis to generate the experimental results

% Erika Kuchen 2018


% figure plotting parameters
fs = 8; % font size
cd('your project path here') % path to working directory of files


%% Correct lineage trees for finite observation time

% read in each dataset in turn
% this section can be adapted to correct the perturbation lineage trees
for m= 1:3 % correct the control replicate experiments
    fates = xlsread('LineageTrees_controlReplicates_Neuroblastoma_KuchenBeckeretal.xlsx',m);
    
    f = [[fates(:,8)+fates(:,4),fates(:,16)];[fates(:,8)+fates(:,7),fates(:,16)]];
    maxgen = max(fates(:,16)); % maximum generation number within dataset
    maxtime = max(f(:,1));% observation duration (cells still alive at the end)
    numalive = sum(f(:,1)>=(maxtime-0.01)); % number of cells alive at observation end (slight shift to account for numerical inaccuracies)
    genalive = zeros(1,maxgen-1);
    allowedsurvival = ceil(numalive*0.05); % number of cells alive at end to get 5% survival
    
    for j = 2:maxgen % for each generation find out how many cells hit the end of the obseration time
        genalive(j-1) = sum((f(:,2)==j).*(f(:,1)>=(maxtime-0.01)));
    end
    
    cumalive = cumsum(genalive);
    geninclude = find(cumalive<=allowedsurvival,1,'last');
    fates = fates(fates(:,16)<=(geninclude+1),:);
    
    save(['CorrectedLineageTree_rep',num2str(m),'.mat'],'fates');
    
end


%% Figure 1 and Figure Supplement 1 and 2


for m=1:3 % iterate over the three replicate control experiments
    load(['CorrectedLineageTree_rep',num2str(m),'.mat']); % load in the observation-time corrected lineage tree data
    
    
    %% Cycle length distribution - Figure 1B and Figure 1 - Figure Supplement 1B
    
    x = 0:1:50;
    d = [fates(:,2);fates(:,5)]; % These are the columns of the cycle lengths of sister 1 and 2, respectively.
    d = d(~isnan(d));
    ch = histc([fates(:,2);fates(:,5)],x);
    figure('Units','centimeters','Position',[5,5,6,2.5])
    bar(x,ch,'k')
    xlim([10,50])
    box off
    title(sprintf('%1.2f (%1.2f,%1.2f)',median(d),quantile(d,0.25),quantile(d,0.75))) % shows the median cycle length and quartiles
    set(gca,'FontSize',fs,'TickDir','out','Xtick',0:2:50)
    
    
    % for model comparison - get confidence bounds on distribution
    % statistics using bootstrap on the level of family trees
    
    rep = 10000; % number of bootstrap repeats
    
    linID = fates(:,1); % lineage ID numbers
    fam = unique(linID);
    n=length(fam);
    pos = randi(n,n,rep); % generate radom indices for sampling
    sumstats = zeros(rep,4);
    distributionStats = zeros(3,4);
    d = [fates(:,2);fates(:,5)];
    distributionStats(1,:) = [nanmean(d), nanmedian(d),quantile(d,0.25),quantile(d,0.75)];
    
    for j=1:rep % draw lineage trees with replacement rep number of times.
        
        d = double.empty;
        
        for k=1:size(pos,1)
            use = linID==fam(pos(k,j));
            d = [d;[fates(use,2);fates(use,5)]];
        end
        sumstats(j,:) = [nanmean(d), nanmedian(d),quantile(d,0.25),quantile(d,0.75)];
    end
    for j=1:4
        distributionStats(2,j) = quantile(sumstats(:,j),0.025);
        distributionStats(3,j) = quantile(sumstats(:,j),0.975);
    end
    
        
    %% Cycle length distribution over time (moving average) - Figure 1C and Figure 1 - Figure Supplement 1C
    % different cutoff points for the different datasets:
    msize = 14; % marker size
    
    tick = 0:130; % find moving average over this time
    btime = [fates(:,[2,8]);fates(:,[5,8])]; % cycle length and birth time of each cell
    use = ~isnan(btime(:,1));
    btime = btime(use,:);
    n = 0;
    av = NaN(length(tick),1);
    int = 10; % bin cells according to their birthtime in 10 hour bins
    for j=tick+int
        n = 1+n;
        use = logical((btime(:,2)>=(j-int)).*(btime(:,2)<=(j+int)));
        av(j) = median(btime(use,1));
    end
    
    med = median(btime(:,1)); % overall median cycle length
    figure('Units','centimeters','Position',[20,20,7.3,3])
    plot([0,110],[med,med],'Color',[0.4,0.4,0.4]); hold on
    plot(btime(:,2),btime(:,1),'.k','MarkerSize',5);
    plot(tick,av(1:end-(int-1)),'r') % median cycle length per bin
    ylim([10,50]);
    xlim([0,115])
    box off
    set(gca,'tickdir','out','FontSize',fs,'YTick',10:10:50)
    xlabel('cell birth time (h)','FontSize',fs)
    ylabel('cell-cycle duration (h)','FontSize',fs)
    
    
    %% Cycle length correlation coefficients  - Figure 1E
    
    % for the replicate experiments
    PlotCorrelations(['rep',num2str(m)])
    
    %% Partial Correlations and randomisation - Figure 1 - Figure Supplement 2C
    
    TemporalDrift_PartialCorrelations(fates,['rep',num2str(m)])
    
    
end


%% Spatial Drift - Figure 1 - Figure Supplement 2D
ms = 7; % marker size
fullinfo =1;

% measurements are in um
scale = 0.4; %(0.4um per pixel for 20x and 0.8um for 10x with the Nikon Ti-Eclipse microscope)
position = 'division'; % plot cell position at 'division' or 'birth'
bs0 = [375,450,375]; %grid sizes for each experiment (rep1,rep2,rep3) 

% load in all three MYCN overexpression replicate experiments and plot side-by-side
dtypes = [{'rep1'},{'rep2'},{'rep3'}];
figure('Units','centimeters','Position',[5,5,14,8])
hold on

for expp = 1:3  
    load(['CorrectedLineageTree_rep',num2str(expp),'.mat']); % load in observation-time corrected trees  
   
    confa = [fates(:,[2,9,11,12,13,1]);fates(:,[5,10,11,14,15,1])]; % columns: cell cycle length, cell ID, mother ID, birth x and y position sibling 1, birth x and y position sibling 2
    mo = repmat(confa(:,3)',size(confa,1),1); % mother ID
    da = repmat(confa(:,2),1,size(confa,1)); % daughter ID
    use = da==mo;
    use = triu(use,1);
    use2 = use';
    [e1,e2,e3] = find(use2);
    % only use the information once:
    [pos, eind] = unique(e2);
    xy = confa(pos,1); % cycle length
    x = confa(e1(eind),4)*scale; % x position
    y = confa(e1(eind),5)*scale; % y position
    famnum = confa(pos,6); % lineage ID number
    
    % plot the cell position at division coloured by lineage ID
    lx = unique(fates(:,1));
    cmap = colormap(jet(max(lx))); % set the colours
    subplot(2,3,expp)
    for j = lx'
        use = famnum==j;
        plot(x(use),y(use),'.','MarkerSize',ms,'Color',cmap(j,:)); hold on
    end
    axis([0,bs0(expp)*0.4*4,0,bs0(expp)*0.4*4])
   
    set(gca,'TickDir','out','XTick',0:100:700,'YTick',0:100:700,'Fontsize',fs);
    
    % plot the cell position and colour by cycle length
    % divide field of view into a grid and check for significantly different cycle distributions within each region.
    switch position
        case 'birth' % position at birth
            x = [fates(:,12);fates(:,14)]*scale; % x-position at birth
            y = [fates(:,13);fates(:,15)]*scale; % y-position at birth
            d = [fates(:,3);fates(:,6)]; % cell death time, where applicable
            xy = [fates(:,2);fates(:,5)]; % cycle length
            famnum = [fates(:,1);fates(:,1)]; % family ID number
        case 'division' % identify the position of division (birth position of daughter cells)
            confa = [fates(:,[2,9,11,12,13,1]);fates(:,[5,10,11,14,15,1])]; % columns: cell cycle length, cell ID, mother ID, birth x and y position sibling 1, birth x and y position sibling 2
            mo = repmat(confa(:,3)',size(confa,1),1);
            da = repmat(confa(:,2),1,size(confa,1));
            use = da==mo;
            use = triu(use,1);
            use2 = use';
            [e1,e2,e3] = find(use2);
            [pos, eind] = unique(e2); % only use the information once
            xy = confa(pos,1);
            x = confa(e1(eind),4)*scale;
            y = confa(e1(eind),5)*scale;
            famnum = confa(pos,6);
    end
    
    % divide this into a grid of 200 pixels.
    bs = bs0(expp)*scale;
    divSpace = 0:bs:(1800*scale);   
    n = length(divSpace)-1;
    meanDiv = zeros(n,n);
    range = 15:22; % cycle length
    col = colormap(jet(length(range)));
    DivTimes = cell(n,n);
    minimumcellnum = 5; %minimum number of cells within a grid point to do stats
    
    count = 0;
    statscomp  = double.empty;
    loopover = 1:n;
    inds = [0,1];
    
    subplot(2,3,expp+3)
    hold on
    for j = loopover % retrieve the cycle lengths within each region and compare the distribution to the full distribution.
        for l = loopover
            count = count+1;
            use1= (x>=divSpace(l+inds(1))).*(x<divSpace(l+inds(2)));
            use2 =  (y>=divSpace(j+inds(1))).*(y<divSpace(j+inds(2)));
            use = logical(use1.*use2);
            ndiv = xy(use);
            ndiv = ndiv(~isnan(ndiv));
            meanDiv(j,l) = mean(ndiv);
            DivTimes{j,l} = ndiv;
            
            if sum(~isnan(xy(use)))>=minimumcellnum
                statscomp = [statscomp; [ndiv,ones(length(ndiv),1)*count]];
            end
            [val,pos] = min(abs(range-repmat(meanDiv(j,l),1,length(range))));
            if isnan(meanDiv(j,l)) || sum(~isnan(xy(use)))<minimumcellnum
                patch([divSpace(l),divSpace(l+1),divSpace(l+1),divSpace(l)],[divSpace(j),divSpace(j),divSpace(j+1),divSpace(j+1)],[1,1,1]);
            else
                patch([divSpace(l),divSpace(l+1),divSpace(l+1),divSpace(l)],[divSpace(j),divSpace(j),divSpace(j+1),divSpace(j+1)],col(pos,:),'FaceAlpha',0.15);
                if fullinfo ==1
                    text((divSpace(l)+divSpace(l+1))/2, (divSpace(j)+divSpace(j+1))/2, num2str(round(mean(ndiv)*10)/10),'Color',[0,0,0],'FontWeight','bold','Fontsize',8)
                end
            end
        end
    end
    %title(num2str(round(range)))
    
    lx = [10:3:23,30,40,80]; % cycle length bins
    cmap = colormap(jet(length(lx)));
    for j = 1:length(xy) % plot each cell as a coloured dot
        switch position
            case 'birth'
                if ~isnan(xy(j))
                    plot(x(j),y(j),'.','MarkerSize',ms,'Color',cmap(logical(histc(xy(j),lx)),:)); hold on
                elseif ~isnan(d(j))
                    plot(x(j),y(j),'.k','MarkerSize',ms); hold on
                end
            case 'division'
                if ~isnan(xy(j))
                    plot(x(j),y(j),'.','MarkerSize',ms,'Color',cmap(logical(histc(xy(j),lx)),:)); hold on
                end
        end
    end
    
    % plotting the colorbar
    if fullinfo ==1
        lxs = cell(length(lx),1);
        for mkstr = 1:length(lx)
            lxs{mkstr} = num2str(lx(mkstr));
        end
        nt = 1/length(lx);
        colorbar('Ticks',nt:nt:1,'TickLabels',lxs);
    end
    
    set(gca,'Fontsize',8,'TickDir','out','XTick',0:100:700,'YTick',0:100:700)
    axis([0 4*bs 0 4*bs])
    
    subplot(2,3,4)
    xlabel('X-position (um)','Fontsize',fs);
    ylabel('Y-position (um)','Fontsize',fs);
    
    
    % convert the regional distributions into an array that can be saved
    % for stats analysis in R (for multiple testing)
    nl = zeros(1,numel(DivTimes));
    for j=1:numel(DivTimes)
        nl(j) = length(DivTimes{j});
    end 
    DivArray = NaN(max(nl),numel(DivTimes));
    for j=1:numel(DivTimes)
      DivArray(1:nl(j),j) = DivTimes{j};  
    end
    csvwrite(['RegionalDistributions_spatial_rep',num2str(expp),'.csv'],DivArray);
end


%% Figure 1 - Figure Supplement 2E
% distance between first cousins and their cycle length correlation

scale = 0.4; %(0.4um per pixel for 20x and 0.8um for 10x with the Nikon Ti-Eclipse microscope)
figure('Units','centimeters','Position',[5,5,14,4])
hold on
for expp = 1:3
    load(['CorrectedLineageTree_rep',num2str(m),'.mat']); % load in observation-time corrected trees
    
    % find first cousins - and then their position at birth
    grandids = repmat(fates(:,19)',size(fates,1),1); %grandmother ID
    cousids = repmat(fates(:,19),1,size(fates,1));
    use4 = cousids==grandids;
    use4 = triu(use4,1);
    use2 = use4';
    [e1,e2,e3] = find(use2);
    fa = [[fates(e1,[2,12,13,1]);fates(e1,[2,12,13,1]);fates(e1,[5,14,15,1]);fates(e1,[5,14,15,1])],...
        [fates(e2,[2,14,15]);fates(e2,[5,14,15]);fates(e2,[2,12,13]);fates(e2,[5,14,15])]];
    fa(:,[2,3,6,7]) = fa(:,[2,3,6,7])*scale;
    % calculate the euclidean distance between the cousin cells
    dist = sqrt((fa(:,2)-fa(:,6)).^2+(fa(:,3)-fa(:,7)).^2);
    % determine difference in cycle length
    diff = fa(:,1)-fa(:,5);
    use = ~isnan(diff);
    dist = dist(use); diff=diff(use); famid = fa(use,4);
    % calculate correlation, to bootstrapping on levels of trees again.
    outd = calcCorrSpearman([dist,diff,famid,famid],10000,'sampleTrees');
    
    subplot(1,3,expp)
    plot(dist,diff,'.k','Markersize',5);
    title(sprintf('\\rho = %1.2f (%1.2f,%1.2f) N=%1.0f',outd),'Fontsize',8);
    set(gca,'TickDir','out','Fontsize',8); box off
    ylim([-20,22])
end
subplot(1,3,1)
xlabel('Euclidean distance first cousins (um)','Fontsize',8)
ylabel('\Delta cycle length (h)','Fontsize',8);



%% Figure 4 and Supplement


%% Cycle length correlation coefficients of the perturbation conditions - Figure 4D and Figure 4 - Figure Supplement 1D

for j = 1:2 % for the two replicate experiments
    
    PlotCorrelations(['-myc',num2str(m)]) % -myc experiments
    PlotCorrelations(['rap',num2str(m)]) % rapamycin experiments
end

%% Figure 4 - Figure Supplement 1A
% RNA-seq of Neuroblastoma TET21N cell populations
% published in Ryl et al. Cell Systems 2017 and available at GEO (GSE98274).

% hard coded, mTOR expression data extracted from source above.
MTOR = [2348	2246	2334	2819]; % column1: MYCN-overexpression experiment B, col2: MYCN-inhibited experiment B, col3: MYCN overexpression experiment C, col4: MYCN-inhibited experiment C

d = [mean(MTOR([1,3])),mean(MTOR([2,4]))]; % calculate mean expression per condition
SEM(1) = std(MTOR([1,3]))/sqrt(length(MTOR([1,3]))); % standard error of the mean per condition
SEM(2) = std(MTOR([2,4]))/sqrt(length(MTOR([2,4])));

figure
bar(d); hold on
errorbar(1,d(1),SEM(1),'k')
errorbar(2,d(2),SEM(2),'k')
set(gca, 'Fontsize',8)
xticklabels({'control','MYCN-inhibited'})
ylabel('mTOR mRNA')
box off

% a standard paired t-test was performed in R (code not provided).



%% Figure 5 and Supplement

%% Cycle length correlation coefficients of the embryonic stem cells  - Figure 5B and Figure 5 - Figure Supplement 1B

for j = 1:3 % for the three replicate experiments
    PlotCorrelations(['esc',num2str(m)])
end