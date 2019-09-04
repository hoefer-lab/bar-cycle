% Approximate Bayesian Computing to determine best fit parameters and
% confidence bounds

% Erika Kuchen 2018

% uniform prior parameter distributions were generated and the model
% simulated.
% here the simulation runs are loaded one at a time and the summary statistics are
% compared to those of the data.

fs = 8; % font size
epsilons = 1:4; % tolerance for ABC.

pars0 = [{'k'},{'sig_g'},{'alpha'},{'gamma'},{'mu'},{'sig_p'}]; % parameter names
intd = [0.0005,0.01,0.02,0.02,0.02,0.02]; % bin size of parameter range tested 
type0 = [{'rep1'},{'rep2'},{'rep3'},{'-myc1'},{'-myc2'},{'rap1'},{'rap2'},{'esc1'},{'esc2'},{'esc3'}]; % dataset tested

% initialse arrays
next = 1;
parBest = NaN(3,6,length(epsilons));
parSet = zeros(10^3,18);
c = zeros(10,13);
lb = c;
ub = c;


% hard coded data summary statistics for the three ESC experiments
% c = cycle length correlation coefficients and distribution statistics
% lb = lower confidence bound
% ub = upper confidence bound
% order: % correlation of the reference cells with the sister, mother,
% grandmother, greatgrandmother, aunt, first cousins, greataunt, cousins
% once removed, second cousins, mean cycle length, median cycle length, first quartile, 3rd quartile

%'rep1'
c(1,:) = [0.7,	0.43,	0.05,	-0.19,	0.27,	0.37,	0.05,	0.26,	0.33,	18.93,	17.83,	16.1667,	20.3333];
lb(1,:) = [0.56	0.28	-0.12	-0.44	0.11	0.15	-0.13	0.01	0.1	18.0445	16.9167	15.625	18.8333];
ub(1,:) = [0.78	0.51	0.21	0.06	0.35	0.57	0.2	0.39	0.48	20.0765	18.8333	17	21.8333];

%'rep2'
c(2,:) = [0.55	0.08	-0.05	0.1	-0.03	0.16	-0.06	0.02	0.09	15.55	15.1	14	16.4];
lb(2,:) = [0.44	0	-0.2	-0.05	-0.15	0.03	-0.16	-0.11	0.01	15.16	14.8	13.8667	15.9667];
ub(2,:) = [0.63	0.16	0.08	0.21	0.07	0.24	0.02	0.12	0.13	16.05	15.47	14.2667	16.8];

%'rep3'
c(3,:) = [0.54	0.17	-0.09	-0.06	0.02	0.25	-0.07	0.02	0.16	17.5314	16.9	15.725	18.475];
lb(3,:) = [0.42	0.06	-0.23	-0.27	-0.08	0.08	-0.24	-0.19	-0.02	17.0928	16.6	15.4	18];
ub(3,:) = [0.65	0.28	0.02	0.13	0.13	0.4	0.09	0.19	0.29	18.0448	17.3	16.1	19.3];

%'-myc1'
c(4,:) = [0.59	0.33	0.19	0.11	0.07	0.18	0.13	0.07	-0.13	24.2	23	20.5	25.8333];
lb(4,:)= [0.48	0.05	-0.01	-0.14	-0.16	0	-0.2	-0.19	-0.51	22.96	22	19.7083	24.6667];
ub(4,:) = [0.72	0.56	0.37	0.37	0.35	0.39	0.34	0.23	0.11	25.67	24	21.6667	28.875];

%'-myc2'
c(5,:) = [0.58	0.26	0.07	-0.03	0.12	0.15	0.06	0	-0.05	24.16	22.8	20	26.6667];
lb(5,:) = [0.46	0.13	-0.12	-0.25	0	-0.04	-0.11	-0.13	-0.21	23.27	22.27	19.4667	25.3333];
ub(5,:) = [0.67	0.35	0.22	0.19	0.22	0.33	0.21	0.11	0.06	25.1	23.73	20.9	28.1333];

%'rap1'
c(6,:) = [0.7	0.43	0.43	0.08	0.22	0.42	0.29	0.25	0.4	24.4847	24	22	26.475];
lb(6,:) = [0.57	0.22	0.21	-0.18	0.04	0.18	0.14	0.02	0.08	23.3292	22.9	21	25.25];
ub(6,:) = [0.8	0.58	0.58	0.31	0.48	0.61	0.54	0.43	0.62	25.9109	25.9	23.1	28];

%'rap2'
c(7,:) = [0.78	0.48	0.51	0.19	0.41	0.6	0.51	0.42	0.53	19.6	18.6	17	20.8];
lb(7,:) = [0.67	0.27	0.26	-0.17	0.19	0.36	0.24	0.18	0.27	18.4	17.8	16	19.8];
ub(7,:) = [0.84	0.58	0.67	0.45	0.51	0.73	0.69	0.53	0.69	20.9	20	18	22.8];

%'esc1'
c(8,:) = [0.74,0.45,0.22,0.14,0.37,0.55,0.29,0.36,0.29,12.0957  , 11.6865,    10.5006 ,  13.5006];
lb(8,:) = [0.53,0.33,0.12,-0.1,0.06,0.35,0.03,0,-0.28, 11.6640  , 11.0442  ,  10.0006   ,13.0006];
ub(8,:) = [0.91,0.53,0.31,0.45,0.52,0.75,0.43,0.5,0.55,12.6289 ,  12.5007   , 11.0006  , 14.0008];

%'esc2'
c(9,:) = [0.79,0.39,0.31,0.12,0.38,0.54,0.25,0.23,0.39,10.9473  , 10.5144   ,9.5035 ,  12.0008];
lb(9,:) = [0.73,0.2,0.03,-0.19,0.21,0.37,0.04,0.01,0.13, 10.3944 ,  10.0143  ,9.2506  , 11.5004];
ub(9,:)  = [0.84,0.58,0.5,0.41,0.56,0.68,0.47,0.48,0.61, 11.6175   ,11.2644 ,10.0006  , 13.0008];

%'esc3'
c(10,:) = [0.60 , 0.34  , 0.2 ,0.26   , 0.37 ,   0.36 ,   0.22 ,0.2, 0.17 ,10.3636 ,  10.0007, 9.0008  , 11.0010];
lb(10,:) = [0.52 ,0.19 ,0.01 ,0.05 ,0.25 ,0.26 ,0.03 ,0.07 ,0.02,   10.0602,    9.9999,     9.0004  , 10.5011];
ub(10,:) = [0.66 ,0.41 ,0.28 ,0.33 ,0.44 ,0.44 ,0.30 ,0.30 ,0.28 , 10.7908,   10.5006   ,9.5008  , 11.5008];

stdC = (ub-lb)/4; % standard deviation for chi2 calculation


all_folders = dir; % find all the simulation runs with different parameter combinations

for j=1:size(all_folders,1) % load the results of each run and calculare the summary statistics
    if ne(all_folders(j).bytes,0) && strcmp(all_folders(j).name(end-3:end),'.mat')
        
        corrStruc = load(all_folders(j).name); % when loading in from R.
        
        p = corrStruc.pars; % parameters
        cm = corrStruc.corr; % correlations
        m = corrStruc.pop; % lineage tree data frame
        
        % generate an array with all the tested parameter combinations (avoid rounding issues)
        l = length(p.sig_p)*length(p.mu);
        k0 = floor(repmat(p.k,l,1)*10000)/10000;
        sig_g0 = floor(repmat(p.sig_g,l,1)*100)/100;
        alpha0 = floor(repmat(p.alpha,l,1)*100)/100;
        gamma0 = floor(repmat(p.gamma,l,1)*100)/100;
        mu7 = repmat(p.mu,1,length(p.sig_p))';
        mu0 = reshape(mu7,numel(mu7),1);
        
        sig_p0 = floor(repmat(p.sig_p,length(p.mu),1)*100)/100;
        sig_p0 = reshape(sig_p0,numel(sig_p0),1);
        
        parSet(next:(next+l-1),1:6) = [k0,sig_g0,alpha0,gamma0,mu0,sig_p0]; % array with all parameter combinations
        
        % calculate chi-squared value between simulation with each parameter set and each dataset
        % all calculated correlation coefficients, distribution of cycle lengths: mean, median, 1st quartile, 3rd quartile
        for dset = 1:length(type0)
            
            chi2all =  ((c(dset,1)-cm.S)./stdC(dset,1)).^2 + ((c(dset,2)-cm.M)./stdC(dset,2)).^2 + ((c(dset,6)-cm.C)./stdC(dset,6)).^2 + ((c(dset,3)-cm.G)./stdC(dset,3)).^2 ...
                + ((c(dset,4)-cm.GG)./stdC(dset,4)).^2  + ((c(dset,9)-cm.C2)./stdC(dset,9)).^2 + ((c(dset,5)-cm.A)./stdC(dset,5)).^2 ...
                + ((c(dset,7)-cm.GA)./stdC(dset,7)).^2  + ((c(dset,8)-cm.CR)./stdC(dset,8)).^2 ...
                +((c(dset,10)-m.mean)./stdC(dset,10)).^2 +((c(dset,11)-m.median)./stdC(dset,11)).^2 + ((c(dset,12)-m.q1)./stdC(dset,12)).^2 +((c(dset,13)-m.q3)./stdC(dset,13)).^2;
            
            parSet(next:(next+l-1),6+dset) = reshape(chi2all,numel(chi2all),1);
            
        end
        next = next+l; 
    end
end

% cut the array if too many rows were initialised
parSet = parSet(sum(parSet,2)>0,:);

%% plot distribution of prior marginals of each parameter

lims = zeros(6,2);

figure
for j =1:6
    lims(j,:) = [min(parSet(:,j)),max(parSet(:,j))]; % find parameter range of simulation 
    xbins = lims(j,1):intd(j):lims(j,2);
    count = histc(parSet(:,j),xbins); %

    subplot(2,3,j)
    bar(xbins,count,'k'); hold on
    xlabel(pars0{j},'Fontsize',fs)
end

%% weighting
% identify local densities of parameter combinations.
% divide the parameter space into a grid and weigh each parameter set by
% the inverse of the number of parameter sets within its grid.

% fine
% coarse
% divide each parameter range into bins, both coarse and fine grids were
% tested.
intk = lims(1,1):0.005:(lims(1,2)+0.005); %0.001 (fine), 0.005 (coarse) % parameter k
intsig_g = lims(2,1):0.01:(lims(2,2)+0.01); %0.01 % parameter sig_g
intalpha = lims(3,1):0.1:(lims(3,2)+0.1); %0.05 (fine), 0.1 (coarse) % parameter alpha
intgamma = lims(4,1):0.1:(lims(4,2)+0.1); %0.05 (fine), 0.1 (coarse) % parameter gamma
intmu = lims(5,1):0.1:(lims(5,2)+0.1); %0.05 (fine), 0.1 (coarse) % parameter mu
intsig_p = lims(6,1):0.1:(lims(6,2)+0.1); %0.05 (fine), 0.1 (coarse) % parameter sig_p

for kp = 2:length(intk)
    use1 = parSet(:,1)<=intk(kp);
    for sgp = 2:length(intsig_g)
        use2 = parSet(:,2)<=intsig_g(sgp);
        for mdp = 2:length(intalpha)
            use3 = parSet(:,3)<=intalpha(mdp);
            for ssp = 2:length(intgamma)
                use4 = parSet(:,4)<=intgamma(ssp);
                for mp = 2:length(intmu)
                    use5 = parSet(:,5)<=intmu(mp);
                    for sp = 2:length(intsig_p)
                        use6 = parSet(:,6)<=intsig_p(sp);
                        use7 = parSet(:,17)==0; % only use those combinations not counted yet.
                        
                        % now count the number of runs within the grid - and give a weight to these cells inverse to the density.
                        posW = use1.*use2.*use3.*use4.*use5.*use6.*use7;
                        weight = sum(posW);
                        if ne(weight,0)
                            parSet(logical(posW),17) = 1/weight;
                        end
                    end
                end
            end
        end
    end
end



%% plot the weighted ABC distributions

parScaling = cell(6,1);

for dj=1:10 % for each dataset
    
    [bf,po] = min(parSet(:,dj+6)); % find the overall minimun chi-squared calculated out of all of the parameter sets
    
    for tol = 1:length(epsilons)  % test for all of the tolerances
        
        parSet(:,18) = parSet(:,dj+6)<=(epsilons(tol)+bf); % find all parameter sets that achieved a chi-squared <= than the overall minimum chi-squared+the tolerance
        
        figure
        for j=1:6 % loop over all of the parameters
            
            parScaling{j}(:,1) = lims(j,1):intd(j):lims(j,2); % x tick marks - bins
            parScaling{j}(:,2)= zeros(size(parScaling{j},1),1);
            
            [count,bin] = histc(parSet(logical(parSet(:,18)),j),parScaling{j}(:,1)); %count how often a specifc parameter value gives good agreement with the data.
            for b = unique(bin)'
                parScaling{j}(b,2) = sum(parSet(bin==b,17)); % sum the weights
            end
            pos =find((cumsum(parScaling{j}(:,2))./sum(parScaling{j}(:,2)))>=0.5,1,'first'); % find the median
            medv = parScaling{j}(pos,1);
            parBest(dj,j,tol) = medv; % median of posterior used as best fit parameter
            
            subplot(2,3,j) % display the output - one figure per tolerance and experiment showing 6 panels for the parameters
            plot(parScaling{j}(:,1),parScaling{j}(:,2),'k'); hold on
            plot([medv,medv],[0,max(parScaling{j}(:,2))],'r') % red line at position of median of posterior distribution
            if j==1
                title(sprintf('dataset %s, tolerance=%1.0f',type0{dj},tol),'Fontsize',fs);
            end
            plim = max(find(ne(parScaling{j}(:,2),0),1,'first')-1,1);
            if ne(parScaling{j}(end,2),0)
              xlim([parScaling{j}(plim,1),parScaling{j}(end,1)]) % confidence bounds have no upper bound
            else
            xlim([parScaling{j}(plim,1),parScaling{j}(find(ne(parScaling{j}(:,2),0),1,'last')+1,1)])
            end
            xlabel(pars0{j},'Fontsize',fs)
            box off
            set(gca,'Fontsize',fs,'TickDir','out');
        end
    end
end


%% Display the parameter confidence bounds as lines

epsilonchosen = 2;

xlimmat = [0.027,0.07; 0,0.1; 0.15,1; 0.4,1; 1.7,3.18; 0.1, 0.7]; % limits of x-axis for each paramater
partick{1} = 0.026:0.004:0.046; % tick marks for parameter k
partick{2} = 0:0.02:0.08; % sig_g
partick{3} = 0.2:0.2:1; % alpha
partick{4} = 0.4:0.1:1; % gamma
partick{5} = 2.6:0.1:3.2; % mu
partick{6} = 0.2:0.1:0.6; % sig_p

orderexp = sort(1:10,'descend'); % reverse the order of the experiments to have 'rep1' at the top.

figure('Units','centimeters','Position',[20,20,14.5,3])
for dj=1:10
    
    [bf,po] = min(parSet(:,dj+6)); % find the overall minimun chi-squared calculated out of all of the parameter sets
    parSet(:,18) = parSet(:,dj+6)<=(bf+epsilonchosen);
    
    for j = 1:6
        subplot(1,6,j)
        
        count = histc(parSet(logical(parSet(:,18)),j),parScaling{j}(:,1)); %count how often a specifc parameter value gives good agreement with the data.
        use = parScaling{j}(count>0,1);
        vs = [min(use),max(use)]; % get boundaries of minimum and maximum of parameter value confidence bounds
        
        % plot confidence bounds as lines and median parameter value (best
        % fit) as a red dot.
        plot(vs,[orderexp(dj), orderexp(dj)],'k','LineWidth',1); hold on
        plot(parBest(dj,j,epsilonchosen),orderexp(dj),'.r','MarkerSize',8)
    end
end

for j=1:6
    subplot(1,6,j)
    xlim(xlimmat(j,:))
    ylim([0.5,10.5])
    box off
    set(gca,'Fontsize',fs,'TickDir','out','Xtick',partick{j},'YTick',1:10,'YTickLabel','');
    xlabel(pars0{j},'Fontsize',fs)
    if j ==1
        set(gca,'YTickLabel',type0(orderexp));
    end
end
