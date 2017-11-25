function main_prefiltering_analysis(data,labels)
% Applies prefiltering analysis to data (dataSamples x dataVariables), labels: dataVariables labels

% suffix = 'prefilteringAnalysis';
suffix = 'prefilteringAnalysis_testing';
if ~exist(suffix,'dir')
  mkdir(suffix)
end

% Figure settings

set(0,'defaultAxesFontName','Arial');
set(0,'defaultTextFontName','Arial');
w = 560;
h = 420;
monitors = get(0,'MonitorPositions');
x = 0.5*(monitors(1,3)-w);
y = 0.5*(monitors(1,4)-h);
set(0,'DefaultFigurePosition',[x y w h]);
resolutionScaling = 96/get(0,'ScreenPixelsPerInch');
fontSize = 16*resolutionScaling;
lineWidth = 2*resolutionScaling;

% Correlation settings
% corrType = 'Pearson';
corrType = 'Spearman';
corrRows = 'pairwise';

% removing measures that do not show any variation
datamin = min(data);
datamax = max(data);
data(:,datamin == datamax) = [];
labels(datamin == datamax) = [];

% normalizalizing within [0,1]
ndata = standardization(data); % norm data

% matrix of correlations
if strcmp(reportUI,'matlab')
  [ndata_corr,ndatacorr_pval] = corr(ndata,'type',corrType,'rows',corrRows);
else % octave
  pkg load nan;
  corrAlpha = 0.05;
  [ndata_corr,ndatacorr_pval] = corrcoef(ndata,'mode',corrType,'alpha',corrAlpha,'rows',corrRows);
end

% clustering measures by correlation
if strcmp(reportUI,'octave')
  pkg load statistics;
end
Z = linkage(ndata','single','spearman');
cutoff = 0.1; % cutoff: correlation of 0.9

figure % ('color','none') % ,'visible','off')
% set(gca,'layer','top','color','none')
[h,t,outperm] = dendrogram(Z,'colorthreshold',cutoff,'labels',labels);
% one way to decide on which measure to keep is the one that is closer to the center in a significant cluster
if strcmp(reportUI,'matlab')
  set(h,'color','k','lineWidth',lineWidth);
end
xl = get(gca,'XTick');
hold on
plot([min(xl)-0.25 max(xl)+0.25],[cutoff cutoff],'r--','lineWidth',lineWidth)
axis([min(xl)-0.25 max(xl)+0.25 0 1])
set(gca,'fontSize',fontSize,'lineWidth',lineWidth,'TickDir','out','Box','off','YTick',0:.2:1)
ylabel('Dissimilarity (Spearman correlation)','fontSize',fontSize);
plot2svg([suffix,'/dataDendrogram.svg']);
% close all;

% variables not highly correlated (all should be considered)
indZcutin = find(Z(:,3) > cutoff);
tmp = Z(indZcutin,1:2);
ind_varcutin = tmp(tmp <= size(data,2));
% variables highly correlated (only one of each correlated set should be considered)
indZcutoff = find(Z(:,3) <= cutoff);
tmp = Z(indZcutoff,1:2);
ind_varcutoff = tmp(tmp(:,1) <= size(data,2) & tmp(:,2) <= size(data,2),2); % ,1); % Lv and Cv_2 are equally fine since correlations are pairwise (selecting Lv appears in the paper)
data_ind_varset = sort([ind_varcutin;ind_varcutoff]);
rawData = data(:,data_ind_varset);
normData = ndata(:,data_ind_varset);
labelsData = labels(data_ind_varset);

% filtering data by variance contribution
var_normData = var(normData);
[sortedVar_normData,indSortedVar_normData] = sort(var_normData,2,'descend');
explainedVar_normData = 100*sortedVar_normData/sum(sortedVar_normData);
cumExplainedVar_normData = cumsum(explainedVar_normData);

% variance cutoff 90%, keeping only measurements that together account for >=90% of the variance
var_cutoff = 90;
indcutoff = find(cumExplainedVar_normData >= var_cutoff,1);
rawData4Cluster = rawData(:,indSortedVar_normData);
rawData4Cluster = rawData4Cluster(:,1:indcutoff);
normData4Cluster = normData(:,indSortedVar_normData);
normData4Cluster = normData4Cluster(:,1:indcutoff);
normData4Cluster_labels = labels(indSortedVar_normData);
normData4Cluster_labels = normData4Cluster_labels(1:indcutoff);
explainedVar_normData4Cluster = explainedVar_normData(1:indcutoff);
cumExplainedVar_normData4Cluster = cumsum(explainedVar_normData4Cluster);

width = 1/3;
figure % ('color','none') % ,'visible','off')
hold on
% set(gca,'layer','top','color','none')
xbar1 = [1:numel(explainedVar_normData)]-0.5*width;
h1=bar(xbar1,explainedVar_normData,'lineWidth',lineWidth);
set(h1,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'BarWidth',width);
xbar2 = [1:numel(cumExplainedVar_normData)]+0.5*width;
h2=bar(xbar2,cumExplainedVar_normData,'lineWidth',lineWidth);
set(h2,'EdgeColor','k','FaceColor','w','BarWidth',width);
xlimits = [0.5 numel(explainedVar_normData)+0.5];
plot(xlimits,[var_cutoff var_cutoff],'r--','lineWidth',lineWidth)
legend([h1 h2],{'Specific','Cumulative'},'fontSize',fontSize,'Location',[.21,.65,.1,.1]); % left, bottom, width, height
legend boxoff;
axis([xlimits 0 105])
set(gca,'yTick',0:20:100)
set(gca,'xTick',1:8)
set(gca,'fontSize',fontSize,'lineWidth',lineWidth,'TickDir','out','Box','off')
xlabels = labelsData(indSortedVar_normData);
set(gca,'XTickLabel',xlabels)
ylabel('Variance explained (%)','fontSize',fontSize)
filename=[suffix,'/Datavariance.svg'];
plot2svg(filename);
% close all;

% matrix of correlations
if strcmp(reportUI,'matlab')
  [normDataClusterCorr,normDataClusterCorr_pval] = corr(normData4Cluster,'type',corrType,'rows',corrRows);
else % octave
  [normDataClusterCorr,normDataClusterCorr_pval] = corrcoef(normData4Cluster,'mode',corrType,'alpha',corrAlpha,'rows',corrRows);
end
