
clear all;
close all;
clc;

%%%% This is work in progress to ease the use of the Clustering Analysis Toolbox and make it compatible with Matlab and GNU Octave %%%%

% Import data (dataSamples x dataVariables):

[data, labels, celltype, datasetname, spikechannel] = import_dataset; % to start your own analysis, edit dataName and load your data in the import_dataset function

% Prefiltering analysis:

main_prefiltering_analysis(data,labels);

%%%% Next functions to improve %%%%

%% main_estimate_numberOfClusters, main_pairingOfClusterElements, main_kmeans_aic_bic, main_clustering_heatMaps, main_dendrogram_clusterCentroids, main_validation and main_validation_monkey %%

%% keeping this for the case %%
% save([suffix,'/clusterDataset.mat'],'datasetname','spikechannel','rawData4Cluster','normData4Cluster','normData4Cluster_labels','explainedVar_normData4Cluster','cumExplainedVar_normData4Cluster','normDataClusterCorr','normDataClusterCorr_pval','celltype','-v7');
