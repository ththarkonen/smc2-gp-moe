
clear
close all

% Include dependencies to path
includeFolders = genpath('include');
addpath( includeFolders );

[ fileName, folderPath] = uigetfile("*.mat");
filePath = fullfile( folderPath, fileName);

predictionSettings = struct();
predictionSettings.nPredictionLocations = 1024;
predictionSettings.nDensityLocations = 1024;
predictionSettings.yPadding = 0.25;

smc2Object = load( filePath );
data = smc2Object.data;

x = data.x;
y = data.y;

xy = [ x, y];
xy = unique( xy, 'rows');

x2 = xy(:,2);
x3 = xy(:,3);

x = xy(:,1:3);
y = xy(:,4);

x2Unique = unique(x2);
x3Unique = unique(x3);
    
x2_ii_jj = x2Unique(28);
x3_ii_jj = x3Unique(1);

inds = x2 == x2_ii_jj;
inds = inds & ( x3 == x3_ii_jj );

xSubset = x( inds, :);
ySubset = y( inds );

xData = xSubset;
yData = ySubset;

smc2Object.data.xSubset = xSubset;
smc2Object.data.ySubset = ySubset;

smc2Object = computePredictionsNasa( smc2Object, predictionSettings);
save( filePath, '-struct', 'smc2Object', '-v7.3');
