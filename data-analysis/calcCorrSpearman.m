function [out,w] = calcCorrSpearman(d,rep,samp)

% calculate the correlations within lineage trees and calculate confidence
% boubds using bootstrap.

% input: [d=[cycle time cell 1, cycle time cell 2, family ID, family ID],
% rep = number of bootstrap repeats, samp=partial or full correlations]]
% samp options: sample on the level of family trees full correlations ('sampleTrees') or partial correlations ('sampleTreesPartial'), sample on the level of individual cells ('sampleAll')

% output: correlations, lower and upper 95% confidence bounds, number of
% cell pairs (all in variable out).

% Erika Kuchen 2017

a = size(d);
%if a(1)>a(2)
d1 = d(:,1);
d2 = d(:,2);
n = a(1);
%else
% d1 = d(1,:);
% d2 = d(2,:);
% n = a(2);
%end

if ~isempty(d1)
    % calculate the correlations for the full dataset
    if strcmp(samp,'sampleTreesPartial') % partial correlations
        coefcalc = partialcorr(d(:,1:3),'Type','Spearman'); %!!!CHECK OUTPUT ERIKA, is this a matlab function?
        coefcalc = coefcalc(2,1);
    else
        coefcalc = corr(d1,d2,'type','Spearman'); % all other cases
    end
    
    [rangecor,w] = dobootstrapTrees(d,rep,samp); % confidence bounds via bootstrap
    out = [coefcalc,rangecor(1,1),rangecor(2,1),n]; % output correlation, lower and upper confidence bounds, number of cell pairs

else
    out = NaN(1,4);
    w = NaN;
end

end