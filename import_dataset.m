function [data, labels, celltype, datasetname, spikechannel] = import_dataset
% Import data (dataSamples x dataVariables)
% to start your own analysis, edit dataName and load your data in this function

% dataName = 'Waveform';
% dataName = 'Activity';
dataName = 'WF+Act';

if ~any( strcmp(dataName,'Waveform') || strcmp(dataName,'Activity') || strcmp(dataName,'WF+Act') )

  % Add your data here:
  % data (samples x variables), variable labels, celltype tags (if you know them), and references of the recording (if you have them): datasetname and spikechannel
  error('Add your data (samples x variables) here!');

elseif strcmp(dataName,'Waveform')

  %%%% Waveform Data (Ardid et al. JNS 2015) %%%%

  load('wFpreprocessed.mat');

  % peak-to-trough duration, time of 25% repolarization and 1st PCA component of the two
  var1 = wFpreprocessed.par(:,1);
  var1 = var1(~isnan(wFpreprocessed.parPCA1stcomp));
  var2 = wFpreprocessed.par(:,3);
  var2 = var2(~isnan(wFpreprocessed.parPCA1stcomp));
  var3 = wFpreprocessed.parPCA1stcomp(~isnan(wFpreprocessed.parPCA1stcomp));
  dataWF = [var1,var2,var3];

  % cell identification
  celltypeWF = wFpreprocessed.celltype(~isnan(wFpreprocessed.parPCA1stcomp));
  datasetnameWF = wFpreprocessed.datasetname;
  datasetnameWF = datasetnameWF(~isnan(wFpreprocessed.parPCA1stcomp));
  spikechannelWF = wFpreprocessed.spikechannel;
  spikechannelWF = spikechannelWF(~isnan(wFpreprocessed.parPCA1stcomp));

  dataWF_labels{1} = 'P2T'; % 'pk2trgh';
  dataWF_labels{2} = 'T4R'; % 't4rplzn';
  dataWF_labels{3} = '1CmpPCA';

  data = dataWF;
  labels = dataWF_labels;
  celltype = celltypeWF;
  datasetname = datasetnameWF;
  spikechannel = spikechannelWF;

elseif strcmp(dataName,'Activity')

  %%%% Activity Data (Ardid et al. JNS 2015) %%%%

  load('Activity.mat');

  % measures
  var1 = RES.FR;
  var2 = RES.FF;
  var3 = RES.Cv;
  var4 = RES.Cv2;
  var5 = RES.Lv;
  var6 = RES.LvR;

  dataACT = [var1(~isnan(var6)),var2(~isnan(var6)),var3(~isnan(var6)),var4(~isnan(var6)),var5(~isnan(var6)),var6(~isnan(var6))];

  % cell identification
  datasetnameACT = RES.datasetname;
  datasetnameACT = datasetnameACT(~isnan(var6));
  spikechannelACT = RES.spikechannel;
  spikechannelACT = spikechannelACT(~isnan(var6));

  for i = 1:length(datasetnameACT)
    tmp = datasetnameACT{i};
    tmp = tmp(1:findstr(tmp,'-')-1);
    datasetnameACT{i} = tmp;
  end

  clear RES;

  dataACT_labels{1} = 'FR';
  dataACT_labels{2} = 'FF';
  dataACT_labels{3} = 'Cv';
  dataACT_labels{4} = 'Cv_2';
  dataACT_labels{5} = 'Lv';
  dataACT_labels{6} = 'LvR';

  data = dataACT;
  labels = dataACT_labels;
  celltype = [];
  datasetname = datasetnameACT;
  spikechannel = spikechannelACT;

elseif strcmp(dataName,'WF+Act')

  %%%% Waveform and Activity Data together (Ardid et al. JNS 2015) %%%%

  load('wFpreprocessed.mat');

  % peak-to-trough duration and time of 25% repolarization (1st PCA component of the two not used for clustering)
  var1 = wFpreprocessed.par(:,1);
  var2 = wFpreprocessed.par(:,3);
  var1 = var1(~isnan(wFpreprocessed.par(:,3))); % some 25% repolarization measures couldn't be done because the recording was cut before
  var2 = var2(~isnan(wFpreprocessed.par(:,3)));
  dataWF = [var1,var2];

  % cell identification
  celltypeWF = wFpreprocessed.celltype(~isnan(wFpreprocessed.par(:,3)));
  datasetnameWF = wFpreprocessed.datasetname;
  datasetnameWF = datasetnameWF(~isnan(wFpreprocessed.par(:,3)));
  spikechannelWF = wFpreprocessed.spikechannel;
  spikechannelWF = spikechannelWF(~isnan(wFpreprocessed.par(:,3)));

  dataWF_labels{1} = 'P2T'; % 'pk2trgh';
  dataWF_labels{2} = 'T4R'; % 't4rplzn';

  %%%% Activity Data (Ardid et al. JNS 2015) %%%%

  load('Activity.mat');

  % measures
  var1 = RES.FR;
  var2 = RES.FF;
  var3 = RES.Cv;
  var4 = RES.Cv2;
  var5 = RES.Lv;
  var6 = RES.LvR;

  dataACT = [var1(~isnan(var6)),var2(~isnan(var6)),var3(~isnan(var6)),var4(~isnan(var6)),var5(~isnan(var6)),var6(~isnan(var6))];

  % cell identification
  datasetnameACT = RES.datasetname;
  datasetnameACT = datasetnameACT(~isnan(var6));
  spikechannelACT = RES.spikechannel;
  spikechannelACT = spikechannelACT(~isnan(var6));

  for i = 1:length(datasetnameACT)
    tmp = datasetnameACT{i};
    tmp = tmp(1:findstr(tmp,'-')-1);
    datasetnameACT{i} = tmp;
  end

  clear RES;

  dataACT_labels{1} = 'FR';
  dataACT_labels{2} = 'FF';
  dataACT_labels{3} = 'Cv';
  dataACT_labels{4} = 'Cv_2';
  dataACT_labels{5} = 'Lv';
  dataACT_labels{6} = 'LvR';

  % putting the two datasets together
  cont = 0;
  maxSize = max([size(dataWF,1),size(dataACT,1)]);
  dataAll = nan(maxSize,size(dataWF,2)+size(dataACT,2));
  for i = 1:length(datasetnameWF)
    for j = 1:length(datasetnameACT)
      if strcmp(datasetnameWF{i},datasetnameACT{j}) && strcmp(spikechannelWF{i},spikechannelACT{j})
        cont = cont+1;
        dataAll(cont,:) = [dataWF(i,:),dataACT(j,:)];
        celltypeAll{cont} = celltypeWF{i};
        datasetnameAll{cont} = datasetnameWF{i};
        spikechannelAll{cont} = spikechannelWF{i};
      end
    end
  end
  dataAll = dataAll(1:cont,:);
  dataAll_labels = [dataWF_labels,dataACT_labels];

  data = dataAll;
  labels = dataAll_labels;
  celltype = celltypeAll;
  datasetname = datasetnameAll;
  spikechannel = spikechannelAll;

end
