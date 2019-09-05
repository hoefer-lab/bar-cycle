function TemporalDrift_PartialCorrelations(fates,expname)

% calculate partial correlations and bootstrap confidence bounds

% Erika Kuchen 2017


rep = 10000; % resampling repeats for bootstrap

% correlate family members
compAll = [{'sibs'},{'mother'},{'granny'},{'greatgranny'},{'aunt'},{'firstcousins'},{'greataunt'},{'counsinsremoved'},{'secondcousins'}];

% calculate partical correlation coefficients with respect to: the last common ancestor
% ('common'), the average birth time of the two relatives ('av'), the birth
% time of the reference cell ('sibs'), the birthtime of the other relative
% ('mother')
birthcorrectionAll = [{'common'},{'av'},{'sibs'},{'mother'}];

partialcorrs = NaN(4,4,9); % initialise array to save correlation coefficients

for rel = 1:9 % go through each pair of relations  
    comp = compAll{rel};
    switch comp
        case 'sibs'
            fa = fates(:,[2,5,8,1,17]); % reference cell cycle length, sister cell cycle length, birth time, lineage ID, mother cell cycle length
            bta = fa(:,3); % sibling birthtime -> birthcorrection 'sibs'
            btan = fa(:,5); % common ancestor -> birthcorrection 'common'
        case 'mother'
            fa = [fates(:,[2,17,8,1,18]);fates(:,[5,17,8,1,18])]; % reference cell cycle length, mother cell cycle length, birth time reference cell, lineage ID, grandmother cell cycle length
            btm = fa(:,3)-fa(:,2); % mother birthtime -> birthcorrection 'mother'
            bts = fa(:,3);         % daughter birthtime -> birthcorrection 'sibs'
            bta = mean([btm,fa(:,3)],2); % average birthtime -> birthcorrection 'av'
            btan = fa(:,5); % common ancestor -> birthcorrection 'common'
        case 'granny'
            fa = [fates(:,[2,18,8,1,17]);fates(:,[5,18,8,1,17])]; % reference cell cycle length, grandmother cell cycle length, birth time reference cell, lineage ID, mother cell cycle length
            btm = fa(:,3)-fa(:,2)-fa(:,5); % grandmother birthtime -> birthcorrection 'mother'
            bts = fa(:,3); % sibling birthtime -> birthcorrection 'sibs'
            bta = mean([btm,fa(:,3)],2); % average birthtime -> birthcorrection 'av'
            btan = NaN(size(fates,1),1);
            for j = 1:size(fates,1) % greatgrandmother division time
                g = fates(fates(:,11)==fates(j,19),18);
                if ~isnan(g)
                    btan(j) = g;
                end
            end
            btan = [btan;btan];
        case 'greatgranny'
            bgran = NaN(size(fates,1),1);
            btan = NaN(size(fates,1),1);
            for j = 1:size(fates,1) % greatgrandmother division time
                g = fates(:,11)==fates(j,19);
                if ne(sum(g),0)
                    bgran(j) = fates(g,18);
                    gid = fates(g,19);
                    gg = fates(fates(:,11)==gid,18);
                    if ~isnan(gg)
                        btan(j) = gg;  % to also get common ancestor with greatgrandmother
                    end
                end
            end
            fa = [fates(:,[2,1,8,1,17,18]);fates(:,[5,1,8,1,17,18])]; % second entry is placeholder
            fa(:,2) = [bgran; bgran];
            btm = fa(:,3)-fa(:,2)-fa(:,5)-fa(:,6); % greatgrandmother birthtime -> birthcorrection 'mother'
            bts = fa(:,3); % sibling birthtime -> birthcorrection 'sibs'
            bta = mean([btm,fa(:,3)],2); % average birthtime -> birthcorrection 'av'
            btan = [btan;btan];
        case 'firstcousins'
            grandids = repmat(fates(:,19)',size(fates,1),1);
            cousids = repmat(fates(:,19),1,size(fates,1));
            use4 = cousids==grandids;
            use4 = triu(use4,1);
            use2 = use4';
            [e1,e2,e3] = find(use2);
            fa0 = [[fates(e1,[2,8,1]);fates(e1,[2,8,1]);fates(e1,[5,8,1]);fates(e1,[5,8,1])],...
                [fates(e2,[2,8,18]);fates(e2,[5,8,18]);fates(e2,[2,8,18]);fates(e2,[5,8,18])]];
            fa = [fa0(:,[1,4]),mean(fa0(:,[2,5]),2),fa0(:,[3,6])]; % average birthtime -> birthcorrection 'av' (only thing that makes sense)
            bta = fa(:,3);
            btan = fa(:,5); % common ancestor -> birthcorrection 'common'
        case 'aunt'
            grandids = repmat(fates(:,19)',size(fates,1),1);
            cousids = repmat(fates(:,19),1,size(fates,1));
            use = cousids==grandids;
            use = triu(use,1);
            use2 = use';
            [e1,e2,e3] = find(use2);
            c1 = [[fates(e1,2);fates(e1,5)],[fates(e2,17);fates(e2,17)]];
            c2 = [[fates(e2,2);fates(e2,5)],[fates(e1,17);fates(e1,17)]];
            c3 = [fates(e1,8);fates(e1,8);fates(e2,8);fates(e2,8)];
            btm = [fates(e1,8)-fates(e1,17);fates(e1,8)-fates(e1,17);fates(e2,8)-fates(e2,17);fates(e2,8)-fates(e2,17)]; % aunt birthtime -> birthcorrection 'mother'
            c4 = [fates(e2,1);fates(e1,1);fates(e2,1);fates(e1,1)];
            fa = [[c1;c2],c3,c4];
            bts = c3;% sibling birthtime -> birthcorrection 'sibs'
            bta = mean([btm,bts],2);  % average birthtime -> birthcorrection 'av'
            btan = repmat(fates(e1,18),4,1); % common ancestor -> birthcorrection 'common'
        case 'greataunt'
            grandids = repmat(fates(:,19)',size(fates,1),1);
            cousids = repmat(fates(:,9),1,size(fates,1));
            use = cousids==grandids;
            use = triu(use,1);
            use2 = use';
            [e1,e2,e3] = find(use2);
            fa = [[fates(e1,2);fates(e1,5)],[fates(e2,5);fates(e2,5)],[fates(e1,8);fates(e1,8)],[fates(e1,1);fates(e1,1)]];
            btm = [fates(e2,8);fates(e2,8)];
            btan = [fates(e2,17);fates(e2,17)];
            
            grandids = repmat(fates(:,19)',size(fates,1),1);
            cousids = repmat(fates(:,10),1,size(fates,1));
            use = cousids==grandids;
            use = triu(use,1);
            use2 = use';
            [e1,e2,e3] = find(use2);
            c2 = [[fates(e1,2);fates(e1,5)],[fates(e2,2);fates(e2,2)],[fates(e1,8);fates(e1,8)],[fates(e1,1);fates(e1,1)]];
            fa = [fa;c2];
            bts = fa(:,3);
            btm = [btm;fates(e2,8);fates(e2,8)];
            bta = mean([btm,bts],2);
            btan = [btan;fates(e2,17);fates(e2,17)];
        case 'counsinsremoved'
            % find the cells that have greataunt as mother
            grandids = repmat(fates(:,19)',size(fates,1),1);
            cousids = repmat(fates(:,9),1,size(fates,1));
            use = cousids==grandids;
            use = triu(use,1);
            use2 = use';
            [e1,e2,e3] = find(use2);
            fa = NaN(1,6);
            for j=1:length(e2)
                %  find the daughter cells, which are the cousins once removed
                use = fates(:,11)==fates(e2(j),10);
                if ne(sum(use),0)
                    addto = [[fates(e1(j),2),fates(use,2),fates(e1(j),8),fates(use,1),fates(use,8),fates(e2(j),17)];...
                        [fates(e1(j),2),fates(use,5),fates(e1(j),8),fates(use,1),fates(use,8),fates(e2(j),17)];...
                        [fates(e1(j),5),fates(use,2),fates(e1(j),8),fates(use,1),fates(use,8),fates(e2(j),17)];...
                        [fates(e1(j),5),fates(use,5),fates(e1(j),8),fates(use,1),fates(use,8),fates(e2(j),17)]];
                    fa = [fa;addto];
                end
            end
            grandids = repmat(fates(:,19)',size(fates,1),1);
            cousids = repmat(fates(:,10),1,size(fates,1));
            use = cousids==grandids;
            use = triu(use,1);
            use2 = use';
            [e1,e2,e3] = find(use2);
            for j=1:length(e2)
                %  find the daughter cells, which are the cousins once removed
                use = fates(:,11)==fates(e2(j),9);
                if ne(sum(use),0)
                    addto = [[fates(e1(j),2),fates(use,2),fates(e1(j),8),fates(use,1),fates(use,8),fates(e2(j),17)];...
                        [fates(e1(j),2),fates(use,5),fates(e1(j),8),fates(use,1),fates(use,8),fates(e2(j),17)];...
                        [fates(e1(j),5),fates(use,2),fates(e1(j),8),fates(use,1),fates(use,8),fates(e2(j),17)];...
                        [fates(e1(j),5),fates(use,5),fates(e1(j),8),fates(use,1),fates(use,8),fates(e2(j),17)]];
                    fa = [fa;addto];
                end
            end
            fa = fa(2:end,:);
            bts = fa(:,3);
            btm = fa(:,5);
            bta = mean([btm,bts],2);
            btan = fa(:,6);
        case 'secondcousins'
            fa = NaN(1,6);
            for j = 1:size(fates,1)
                use1 = find(fates(:,19)==fates(j,9));
                use2 = find(fates(:,19)==fates(j,10));
                if ne(sum(use1),0) && ne(sum(use2),0)
                    for l=use1'
                        for m=use2'
                            addto = [[fates(l,2);fates(l,2);fates(l,5);fates(l,5)],[fates(m,2);fates(m,5);fates(m,2);fates(m,5)],...
                                [fates(l,8);fates(l,8);fates(l,8);fates(l,8)],[fates(l,1);fates(l,1);fates(l,1);fates(l,1)],[fates(m,8);fates(m,8);fates(m,8);fates(m,8)],...
                                [fates(j,17);fates(j,17);fates(j,17);fates(j,17)]];
                            fa = [fa;addto];
                        end
                    end
                end
            end
            fa = fa(2:end,:);
            bta = mean([fa(:,3),fa(:,5)],2);
            btan = fa(:,6);
    end
    
    % now calculate the correlation between the vectors with family member cycle and birth times generated above.
    use = ~isnan(fa(:,1).*fa(:,2)); % cut out NaNs.
    fa = fa(use,:);
    
    for cpc = 1:length(birthcorrectionAll)
        birthcorrection = birthcorrectionAll{cpc};
        if cpc>2 && any([1,6,9]==rel)
            break
        end
        fa2 = fa;
        switch birthcorrection
            case 'mother' % older relative birthtime
                fa2(:,3) = btm(use);
            case 'sibs' % sibling birthtime
                fa2(:,3) = bts(use);
            case 'av' % average birthtime
                fa2(:,3) = bta(use);
            case 'common' % common ancestor
                fa2(:,3) = btan(use);
                use2 = ~isnan(fa2(:,3));
                fa2 = fa2(use2,:);
        end
        partialcorrs(cpc,:,rel) = calcCorrSpearman(fa2,rep,'sampleTreesPartial');
        
    end
end

%% randomise cells to test for correlation

reps = 10000; % number of randomisation repeats

out = zeros(1,reps);

for j = 1:reps
    sub = fates(:,[2,5]);
    s = reshape(sub,numel(sub),1);
    s = s(~isnan(s));
    s1 = s(randperm(length(s)));
    if rem(length(s1),2)
        s1 = s1(2:end);
    end
    pairs = reshape(s1,length(s1)/2,2);
    out(j) = corr(pairs(:,1),pairs(:,2),'type','Spearman');
end
totalrand = [nanmedian(out),quantile(out,0.025),quantile(out,0.975)]; % get the median Spearman rank correlation coefficient and the 95% confidence bounds.



%% Plot the data

% plot the normal Spearman rank corelations, the partial correlations and
% the total randomisation in one figure.

% plotting parameters
ms = [8,2,2,2]; % marker size
ofs = 0.15; % space between data points
col = [1,0,0; 0,1,0; 0,0,1; 0.6,0.6,0]; % colours
marker = ['.','d','^','v']; % marker types

% plot the previously calculated overall Spearman rank correlations
PlotCorrelations(expname) 


for j=1:9 % now add the Spearman rank partial correlations
    if j==1 % plot the randomisation
        plot([j,j],totalrand(2:3),'Color',[0.4,0.4,0.4],'LineWidth',1);
        plot(j,totalrand(1),'x','MarkerSize',4,'Color',[0.4,0.4,0.4]); % everything is randomised
        plot([j,j]+ofs,temprand(2:3),'Color',[0.7,0.7,0.7]);
        plot(j+ofs,temprand(1),'.','MarkerSize',4,'Color',[0.4,0.4,0.4]);
    end
    for l = 1:4
        plot([j,j]+ofs*l,partialcorrs(l,2:3,j),'Color',col(l,:),'LineWidth',1);
        plot(j+ofs*l,partialcorrs(l,1,j),marker(l),'MarkerSize',ms(l),'Color',col(l,:),'MarkerFaceColor',col(l,:)); %
    end
end
ylim([-0.55,0.8])

set(gca,'Fontsize',8,'YTick',-0.4:0.2:0.8);

