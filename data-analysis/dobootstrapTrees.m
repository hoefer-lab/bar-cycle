function [cor,weight] = dobootstrapTrees(v,rep,samp) % ERIKA - change to make runnable by statistics toolbox for Spearman correlations!

% performs booststrap sampling for correlation calculations in lineage
% trees 
% Erika Kuchen 2017

corn = zeros(3,rep);
v1 = v(:,1); v2 = v(:,2);

switch samp
    case 'sampleTrees' % sample on the level of lineage trees
        v4 = v(:,4); % lineage ID numbers
        fam = unique(v4);
        n=length(fam);
        pos = randi(n,n,rep); % generate radom indecies for sampling 
        
        for j=1:rep % draw lineage trees with replacement rep number of times. 
            
            sall = double.empty;
            
            for k=1:size(pos,1)
                use = v4==fam(pos(k,j));
                sall = [sall;[v1(use),v2(use)]]; 
            end
            if size(sall,1)>=2
                [coef,pvalc] = corr(sall(:,1),sall(:,2),'type','Spearman'); 
                % if statistics toolbox is missing, use code below for
                % Pearson correlations 
                %coP = corrcoef(sall(:,1),sall(:,2));
                %coef = coP(1,2); pvalc = coP(1,1); %pval not correct assignment here, just to fill the gap
                corn(:,j) = [coef, pvalc,size(sall,1)];
            else
                corn(:,j) = [NaN,NaN,NaN];
            end
        end
        
    case 'sampleAll' % sample on the level of individual cells
        n = length(v1);
        pos = randi(n,n,rep); % generate radom indecies for sampling
        
        for j=1:rep
            s1 = v1(pos(:,j));
            s2 = v2(pos(:,j));
            [coef,pvalc] = corr(s1,s2,'type','Spearman');
            corn(:,j) = [coef, pvalc, size(s1,1)];
            % if statistics toolbox is missing, use code below for Pearson correlations 
            %[coef,pvalc] = corrcoef(s1,s2);
            %corn(:,j) = [coef(2,1), pvalc(2,1),size(s1,1)];
        end
        
    case 'sampleTreesPartial' % sample on the level of lineage trees and calculate partial correlations
        v4 = v(:,4);
        v3 = v(:,3);
        fam = unique(v4);
        n=length(fam);
        pos = randi(n,n,rep);
        
        for j=1:rep
            sall = double.empty;
            for k=1:size(pos,1)
                use = v4==fam(pos(k,j));
                sall = [sall;[v1(use),v2(use),v3(use)]];
            end
            if size(sall,1)>=2
                cout = partialcorr(sall,'Type','Spearman'); 
                corn(:,j) = [cout(2,1), NaN,size(sall,1)];
            else
                corn(:,j) = [NaN,NaN,NaN];
            end
        end
end

% with MATLAB statistics toolbox:
weight = nanvar(corn(1,:)); 
a1 = quantile(corn(1,:),0.025);
a2 = quantile(corn(1,:),0.975);
% without MATLAB statistics toolbox:
%weight = 1; 
%[scor,ind] = sort(corn(1,:));
%nrep = ~isnan(scor);
%scor = scor(nrep);
%ind = ind(nrep);
%pcor = corn(2,ind);
%minmax = [floor(sum(nrep)*0.025), ceil(sum(nrep)*0.975)];
%a2 = [scor(minmax(2)), pcor(minmax(2))];
%a1 = [scor(minmax(1)), pcor(minmax(1))];
 
cor = [a1; a2];