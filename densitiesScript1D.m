
clear
close all

% Include dependencies to path
includeFolders = genpath('include');
addpath( includeFolders );

[ fileName, folderPath] = uigetfile("*.mat");
filePath = fullfile( folderPath, fileName);

nX = 1024;
nY = 1024;

predictionSettings = struct();
predictionSettings.nPredictionLocations = 1024;
predictionSettings.nDensityLocations = 1024;
predictionSettings.yPadding = 0.25;

smc2Object = load( filePath );
smc2Object = computePredictions1D( smc2Object, predictionSettings);
save( filePath, '-struct', 'smc2Object', '-v7.3');