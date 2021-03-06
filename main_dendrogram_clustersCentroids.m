
clear all;
close all;
clc;

type = 'Kmeans';
resultsdir = 'clusteringResults';
matfile = [resultsdir,'/',type,'ClusteringResults'];
load(matfile);

% calculating centers of k-means clusters
numClusters = 7;
offset = 4; % it starts with numClusters = 5
numOfDim=size(datanorm4Cluster,2);

nC=length(clustFilt(1,:));
k=0;
tmp=1:size(datanorm4Cluster,1);
for i=1:nC
    if(~isempty([clustFilt{numClusters-offset,i}]))
        k = k+1;
        centers(k,:) = mean(datanorm4Cluster([clustFilt{numClusters-offset,i}],:),1);
        centersSD(k,:) = std(datanorm4Cluster([clustFilt{numClusters-offset,i}],:),0,1);
        numElem(k) = length([clustFilt{numClusters-offset,i}]);
        tmp = setdiff(tmp,[clustFilt{numClusters-offset,i}]);
    end
end

% remaining elements if any
if ~isempty(tmp)
  k = k+1;
  clustFilt{numClusters-offset,k} = tmp;
  numElem(k) = size(datanorm4Cluster,1)-sum(numElem);
  centers(k,:) = mean(datanorm4Cluster([clustFilt{numClusters-offset,k}],:),1);
  centersSD(k,:) = std(datanorm4Cluster([clustFilt{numClusters-offset,k}],:),0,1);
end

Ap = [];
for i = 1:k
  % I add a term only to keep the proper structure of the dendrogram according to the neurons in the clusters
  Ap = [Ap;(ones(numElem(i),1)+1e-7*((1:numElem(i))'/numElem(i)-0.5))*centers(i,:)];
end
Y = pdist(Ap);
Z = linkage(Y,'average');
c = cluster(Z,'maxclust',k);

figure('color','none','visible','off');
set(gca,'layer','top','color','none')

% colors
hexcolors={'1D0091','007BFF','03FC5E','FFFF00','FFAE00','FF4000','000000','A8100D'};
% alternatively one can use linspecer to make this automatically distinguishable for up to 12 clusters:
% colorOrder = linspecer(k,'qualitative');
colorOrder=zeros(length(hexcolors),3);
for i=1:length(hexcolors)
  colorOrder(i,:) = rgbconv(hexcolors{i});
end

% threshold color
t = sort(Z(:,3));
Cth = t(size(Z,1)+2-k);

[h,~,outperm] = dendrogram(Z, 0,'ColorThreshold',Cth,'orientation','left');
hcolors = cell2mat(get(h,'Color'));
hy = cell2mat(get(h,'YData'));

% sorted index to colored branches
[~,hy_order] = sort(max([hy(1:end-k+1,1) hy(1:end-k+1,end)],[],2));

% replacing colors
iColor = 1;
colorGroup = hcolors(1,:);
for kk = 1:numel(hy_order)
  % checking if we are in a new color branch
  if ~ismember(colorGroup,hcolors(hy_order(kk),:),'rows')
    colorGroup = [colorGroup;hcolors(hy_order(kk),:)];
    iColor = iColor + 1;
  end
  % replacing each branch color by the desired color sequence from colorOrder
  set(h(hy_order(kk)), 'Color', colorOrder(iColor,:));
end
set(h,'LineWidth',1)
axis off
hold on
ylabel('Neurons','fontsize',16);
xlabel('Distance','fontsize',16);
set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off','XTick',0:.2:1,'YTick',[])
fig2svg([resultsdir,'/dendrogram_',type,'_centroids.svg'])

%% matching the dendrogram with the results of the kmeans to determine the cluster ordering
%% doing so by matching the number of elements (warning!!! it is not reliable if two or more clusters have the same number of elements)
% clustersNumOfElements = find(diff(c(outperm))~=0);
% clustersNumOfElements = [clustersNumOfElements(1);diff([clustersNumOfElements;length(c)])];
% for i = 1:length(clustersNumOfElements)
%   clustOrder(i) = find(clustersNumOfElements(i)==numElem);
% end

%% matching the dendrogram with the results of the kmeans to determine the cluster ordering
%% doing so by checking the distance between centers (note: a reliable way to do it and useful for debugging)
% for i = 1:size(centers,1)
%   X = centers;
%   Y = mean(Ap(c==i,:),1);
%   idx(i) = find(pdist2(X,Y,'euclidean') <= 1e-6);
% end

%% matching the dendrogram with the results of the kmeans to determine the cluster ordering
%% the proper way to do it remapping indexes
[~,idx] = sort(unique(c','stable'));

clustOrder = idx(unique(c(outperm)','stable'));

save(matfile,'clustFilt','clustOrder','colorOrder','-append');
