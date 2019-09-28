% Recursive Least Squares Dictionary Learning - Image Inpainting Analysis

%% General Parameters and Initialization

clear();
close('all');

set(0, 'DefaultFigureWindowStyle', 'docked');
defaultLooseInset = get(0, 'DefaultAxesLooseInset');
set(0, 'DefaultAxesLooseInset', [0.05, 0.05, 0.05, 0.05]);

titleFontSize   = 14;
axisFotnSize    = 12;
stringFontSize  = 12;

thinLineWidth   = 2;
normalLineWidth = 3;
thickLineWidth  = 4;

smallSizeData   = 36;
mediumSizeData  = 48;
bigSizeData     = 60;

randomNumberStream = RandStream('mlfg6331_64', 'NormalTransform', 'Ziggurat');
subStreamNumber = 57162;
set(randomNumberStream, 'Substream', subStreamNumber);
RandStream.setGlobalStream(randomNumberStream);

addpath(genpath('RawData'));
addpath(genpath('AuxiliaryFunctions'));


%% Setting Constants

FALSE   = 0;
TRUE    = 1;

OFF = 0;
ON  = 1;


%% Setting Analysis Parameters


%% Loading Data


%% Image Inpainting - ILS-DLA


%% Image Inpainting - K-SVD-DLA


%% Image Inpainting - RLS-DLA


%% Display Results

hFigure         = figure();
hAxes           = axes();
hBarSeries     = bar(vLsRmsError);
set(get(hAxes, 'Title'), 'String', ['Category - RMS Error'], ...
    'FontSize', titleFontSize);
set(get(hAxes, 'XLabel'), 'String', 'Category', ...
    'FontSize', axisFotnSize);
set(hAxes, 'XTick', [1:(numCategories + 1)]);
set(hAxes, 'XTickLabel', xTickLabelStr);
set(hAxes, 'XLim', [0.5, (numCategories + 1.5)]);
set(get(hAxes, 'YLabel'), 'String', 'Error [$]', ...
    'FontSize', axisFotnSize);


%% Restore Defaults
set(0, 'DefaultFigureWindowStyle', 'normal');
set(0, 'DefaultAxesLooseInset', defaultLooseInset);

