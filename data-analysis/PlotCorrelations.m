function PlotCorrelations(expname)

% Make a figure of the correlation coefficients within lineage trees.

% Erika Kuchen 2017


% hard coded previously calculated correlation coefficients and population statistics (vector c) and confidence
% bounds (lower bound vector 'lb' and upper bound vector 'ub') for the
% different datasets specified by 'expname'.

% 'c' contains information on the correlations of the reference cells with
% their sister, mother, grandmother, greatgrandmother, aunt, first cousins,
% greataunt, cousins once removed, second cousins, mean cycle length,
% median cycle length, ... 

% plotting parameters 
l = 1; % line width 
fsize = 8; % font size 

switch expname  
    case 'rep1'
        c = [0.7,	0.43,	0.05,	-0.19,	0.27,	0.37,	0.05,	0.26,	0.33,	18.93,	17.83,	16.1667,	20.3333];
        lb = [0.56	0.28	-0.12	-0.44	0.11	0.15	-0.13	0.01	0.1	18.0445	16.9167	15.625	18.8333];
        ub = [0.78	0.51	0.21	0.06	0.35	0.57	0.2	0.39	0.48	20.0765	18.8333	17	21.8333];
        
    case 'rep2'
        c = [0.55	0.08	-0.05	0.1	-0.03	0.16	-0.06	0.02	0.09	15.55	15.1	14	16.4];
        lb = [0.44	0	-0.2	-0.05	-0.15	0.03	-0.16	-0.11	0.01	15.16	14.8	13.8667	15.9667];
        ub = [0.63	0.16	0.08	0.21	0.07	0.24	0.02	0.12	0.13	16.05	15.47	14.2667	16.8];
        
    case 'rep3'
        c = [0.54	0.17	-0.09	-0.06	0.02	0.25	-0.07	0.02	0.16	17.5314	16.9	15.725	18.475];
        lb = [0.42	0.06	-0.23	-0.27	-0.08	0.08	-0.24	-0.19	-0.02	17.0928	16.6	15.4	18];
        ub = [0.65	0.28	0.02	0.13	0.13	0.4	0.09	0.19	0.29	18.0448	17.3	16.1	19.3];
        
    case '-myc1'
        c = [0.59	0.33	0.19	0.11	0.07	0.18	0.13	0.07	-0.13	24.2	23	20.5	25.8333];
        lb= [0.48	0.05	-0.01	-0.14	-0.16	0	-0.2	-0.19	-0.51	22.96	22	19.7083	24.6667];
        ub = [0.72	0.56	0.37	0.37	0.35	0.39	0.34	0.23	0.11	25.67	24	21.6667	28.875];
        N = [76   160   125    72   137   102   104   170   131]; % number of cell pairs
        
    case '-myc2'
        c = [0.58	0.26	0.07	-0.03	0.12	0.15	0.06	0	-0.05	24.16	22.8	20	26.6667];
        lb = [0.46	0.13	-0.12	-0.25	0	-0.04	-0.11	-0.13	-0.21	23.27	22.27	19.4667	25.3333];
        ub = [0.67	0.35	0.22	0.19	0.22	0.33	0.21	0.11	0.06	25.1	23.73	20.9	28.1333];
        N = [165   337   264   156   296   248   244   401   325];
        
    case 'rap1'
        c = [0.7	0.43	0.43	0.08	0.22	0.42	0.29	0.25	0.4	24.4847	24	22	26.475];
        lb = [0.57	0.22	0.21	-0.18	0.04	0.18	0.14	0.02	0.08	23.3292	22.9	21	25.25];
        ub = [0.8	0.58	0.58	0.31	0.48	0.61	0.54	0.43	0.62	25.9109	25.9	23.1	28];
        
    case 'rap2'
        c = [0.78	0.48	0.51	0.19	0.41	0.6	0.51	0.42	0.53	19.6	18.6	17	20.8];
        lb = [0.67	0.27	0.26	-0.17	0.19	0.36	0.24	0.18	0.27	18.4	17.8	16	19.8];
        ub = [0.84	0.58	0.67	0.45	0.51	0.73	0.69	0.53	0.69	20.9	20	18	22.8];
        
    case 'esc1' 
        c = [0.74,0.45,0.22,0.14,0.37,0.55,0.29,0.36,0.29,12.0957  , 11.6865,    10.5006 ,  13.5006];
        lb = [0.53,0.33,0.12,-0.1,0.06,0.35,0.03,0,-0.28, 11.6640  , 11.0442  ,  10.0006   ,13.0006];
        ub = [0.91,0.53,0.31,0.45,0.52,0.75,0.43,0.5,0.55,12.6289 ,  12.5007   , 11.0006  , 14.0008];
        
    case 'esc2' 
        c = [0.79,0.39,0.31,0.12,0.38,0.54,0.25,0.23,0.39,10.9473  , 10.5144   ,9.5035 ,  12.0008];
        lb = [0.73,0.2,0.03,-0.19,0.21,0.37,0.04,0.01,0.13, 10.3944 ,  10.0143  ,9.2506  , 11.5004];
        ub  = [0.84,0.58,0.5,0.41,0.56,0.68,0.47,0.48,0.61, 11.6175   ,11.2644 ,10.0006  , 13.0008];
        
    case 'esc3'
        c = [0.600785471683363 , 0.340635863092725  , 0.2 ,0.26   , 0.37 ,   0.36 ,   0.22 ,0.2, 0.17 ,10.3636 ,  10.0007, 9.0008  , 11.0010];
        lb = [0.52 ,0.19 ,0.01 ,0.05 ,0.25 ,0.26 ,0.03 ,0.07 ,0.02,   10.0602,    9.9999,     9.0004  , 10.5011];
        ub = [0.66 ,0.41 ,0.28 ,0.33 ,0.44 ,0.44 ,0.30 ,0.30 ,0.28 , 10.7908,   10.5006   ,9.5008  , 11.5008];
end


figure('Units','centimeters','Position',[20,20,13,10])
hold on
plot([0,10],[0,0],'Color',[0.7,0.7,0.7]) % plot zero line 
for j=1:9
    plot(j,c(j),'sk','MarkerSize',3,'MarkerFaceColor','k'); % confidence bounds of correlation coefficient 
    plot([j,j],[lb(j),ub(j)],'k','LineWidth',l); % correlation coefficient 
end
ylabel('Spearman rank correlation coefficient','FontSize',fsize)
set(gca,'XTick',1:length(c),'TickDir','out','XTickLabel',[{'Sibs'},{'MoDa'},{'GranDa'},{'GreatGran'},{'Aunt'},{'Cous'},{'Gaunt'},{'CousRem'},{'2nd Cous'}],'FontSize',8);
box off
title(sprintf('%s',expname))
xlim([0.5,9.5]);ylim([-0.52,0.86])
