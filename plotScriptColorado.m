
clear
close all

% Include dependencies to path
includeFolders = genpath('include');
addpath( includeFolders );

[ fileName, folderPath] = uigetfile("*.mat");
filePath = fullfile( folderPath, fileName);
psmPath = fullfile( folderPath, "psm.txt");

result = load( filePath );
PSM = readmatrix( psmPath );

X = result.densities.grids.X;
Y = result.densities.grids.Y;
C = result.densities.data;

data = result.data;
inds = ( data.x(:,2) <= 1 ) & ( data.x(:,2) >= 0.0 ); 

[ nY, nX] = size( X );

xStar = X(1,:);
xStar = xStar(:);
yStar = Y(:,1);

highDensityMatrix = computeHighDensityArea( result );

figure();
hold on

hData = plot( result.data.x( inds, 1), result.data.y( inds ), '.');
hData.MarkerSize = 45;

imagesc( xStar, yStar, highDensityMatrix)

cmap = flipud( gray );
colormap( cmap );

h = plot( xStar, result.quantiles.data.median);
h.Color = [ 51, 51, 51] / 255;
h.LineWidth = 5;

hData = plot( result.data.x( inds, 1), result.data.y( inds ), '.', "Color", hData.Color);
hData.MarkerSize = 45;

xlim([0,1])
ylim([min(yStar), max(yStar)])

h = xlabel("$x$");
h.Interpreter = "latex";

h = ylabel("$y$");
h.Interpreter = "latex";

ax = gca();
ax.Layer = "top";
ax.YDir = "normal";
ax.FontSize = 40;

clim([0, 4])
pbaspect([ 1920, 1080, 1])

figure();

image( 255 * flipud( PSM ) );
colormap( cmap );

ax = gca();
ax.XTick = [];
ax.YTick = [];
axis tight
axis square
ax.Color = "none";
ax.Color = 'none';
ax.Visible = "off";

figure();

h = bar( result.posteriorNumberOfExperts );

h = xlabel("Non-empty experts");
h.Interpreter = "latex";

axis square

ax = gca();
ax.Box = "off";
ax.FontSize = 50;
