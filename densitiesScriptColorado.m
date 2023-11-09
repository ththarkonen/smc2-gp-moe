
clear
close all

% Include dependencies to path
includeFolders = genpath('include');
addpath( includeFolders );

[ fileName, folderPath] = uigetfile("*.mat");
filePath = fullfile( folderPath, fileName);

nX = 1024;
nY = 2048;

smc2Object = load( filePath );
data = smc2Object.data;
    
heightInterpolation = scatteredInterpolant( data.x(:,1), data.x(:,2), data.x(:,3));

lonStar = linspace( 0, 1, nX)';
lanStar = 0.5 + 0 * lonStar;
heigthStar = heightInterpolation( lonStar, lanStar);
timeStar = 1 / 2 + 0 * lonStar;

xStar = [ lonStar, lanStar, heigthStar, timeStar];

predictionSettings = struct();
predictionSettings.nPredictionLocations = nX;
predictionSettings.nDensityLocations = nY;
predictionSettings.yPadding = 1.50;
predictionSettings.xStar = xStar;

smc2Object = computePredictionsColorado( smc2Object, predictionSettings);
save( filePath, '-struct', 'smc2Object', '-v7.3');