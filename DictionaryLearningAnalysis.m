% Recursive Least Squares Dictionary Learning - Dictionary Learning Analysis

%% General Parameters and Initialization

clear();
close('all');

% set(0, 'DefaultFigureWindowStyle', 'docked');
defaultLooseInset = get(0, 'DefaultAxesLooseInset');
set(0, 'DefaultAxesLooseInset', [0.05, 0.05, 0.05, 0.05]);

titleFontSize   = 14;
axisFotnSize    = 12;
stringFontSize  = 12;

thinLineWidth   = 2;
normalLineWidth = 3;
thickLineWidth  = 4;

smallMarkerSize     = 6;
mediumMarkerSize    = 8;
largeMarkerSize     = 10;

smallSizeData   = 36;
mediumSizeData  = 48;
bigSizeData     = 60;

randomNumberStream = RandStream('mlfg6331_64', 'NormalTransform', 'Ziggurat');
subStreamNumber = 57162;
subStreamNumber = 59132;
% subStreamNumber = round(sum(clock()));
set(randomNumberStream, 'Substream', subStreamNumber);
RandStream.setGlobalStream(randomNumberStream);

addpath(genpath('RawData'));
addpath(genpath('AuxiliaryFunctions'));
addpath(genpath('../AuxiliaryFunctions'));


%% Setting Constants

FALSE   = 0;
TRUE    = 1;

OFF = 0;
ON  = 1;

NUM_TO_PCT_FACTOR = 100;

OMP_MODE_CARDINALITY    = 1;
OMP_MODE_NOISE_THR      = 2;

LAMBDA_FACTOR_MODE_POLYNOMIAL   = 1;
LAMBDA_FACTOR_MODE_HYPERBOLA    = 2;
LAMBDA_FACTOR_MODE_EXPONENTIAL  = 3;

RLS_CONFIG_CONSTANT         = 1;
RLS_CONFIG_LINEAR_100       = 2;
RLS_CONFIG_LINEAR_085       = 3;
RLS_CONFIG_LINEAR_075       = 4;
RLS_CONFIG_LINEAR_050       = 5;
RLS_CONFIG_QUADRATIC_100    = 6;
RLS_CONFIG_QUADRATIC_085    = 7;
RLS_CONFIG_QUADRATIC_075    = 8;
RLS_CONFIG_CUBIC_100        = 9;
RLS_CONFIG_CUBIC_085        = 10;
RLS_CONFIG_CUBIC_075        = 11;
RLS_CONFIG_HYPERBOLA        = 12;
RLS_CONFIG_EXP              = 13;


%% Setting Analysis Parameters

simString = 'Dictionary Analysis';

% Parameters for Data Generation
numAtoms        = 128;
atomLength      = 64;
dataCardinality = 5;
numSamples      = 8192;
noiseStd        = 0.0;

% Parameters for Analysis
dataCardinalityDla  = 5;
noiseStdThr         = 1.01 * sqrt((noiseStd ^ 2) * atomLength);
% ompMode             = OMP_MODE_NOISE_THR;
ompMode             = OMP_MODE_CARDINALITY;
dicErrThr           = 0;
simThr              = 0.99;
rlsDlaConfig        = RLS_CONFIG_LINEAR_050;

switch(ompMode)
    case(OMP_MODE_CARDINALITY)
        ompThr       = dataCardinalityDla;
    case(OMP_MODE_NOISE_THR)
        ompThr       = noiseStdThr;
end

% Lambda:
[lambdaFctrInitialValue, lambdaFctrEffIterFactor, lambdaFctrPwrFctr, lambdaFctrMode] = GetLambdaParams(rlsDlaConfig);

rldDlaRepFctr = 4;
numIterations = rldDlaRepFctr * numSamples;
numEffIter    = round(lambdaFctrEffIterFactor * numIterations);
vLambdaFctr   = CreateLambdaFactor(lambdaFctrInitialValue, numIterations, numEffIter, lambdaFctrPwrFctr, lambdaFctrMode);


%% Generating Data

% Reference Dictionary
mD      = CreateRandomDictionary(atomLength, numAtoms);
% Dictionary for Initialization
mDInit  = CreateRandomDictionary(atomLength, numAtoms);

% Create Data - mX = mD * mA
[mX, mA] = CreateDataSamples(mD, dataCardinality, numSamples, noiseStd);


cMethodString = {'ILS DLA', 'K-SVD DLA', 'RLS DLA'};
cMethodDla    = {@ApplyIlsDla, @ApplyKSvdDla, @ApplyRlsDla};
% cMethodString = {'RLS DLA'};
% cMethodDla    = {@ApplyRlsDla};


%% Dictionary Analysis

parfor ii = 1 : length(cMethodDla)
    methodString = cMethodString{ii};
    disp([' ']);
    disp([simString, ' - ', methodString]);
    disp([' ']);

    hTicTimer               = tic();
    RunDla                  = cMethodDla{ii};
    switch (methodString)
        case({'ILS DLA', 'K-SVD DLA'})
            cParams = {mX, mDInit, ompThr, ompMode, dicErrThr};
        case('RLS DLA')
            cParams = {repmat(mX, [1, rldDlaRepFctr]), mDInit, vLambdaFctr, ompThr, ompMode, dicErrThr};
    end
    [mTrainD, mTrainA, vDlaErr{ii}]  = RunDla(cParams{:});
    
    runTime = toc(hTicTimer);
    
    disp([' ']); disp([simString, ' - ', methodString, ' - Run Time - ', num2str(runTime), ' [Sec]']);
    disp([' ']);

    %-- Evaluating Results
    mMatchAtoms{ii}  = AnalyzeDictionarySimilarity(mD, mTrainD, simThr);
end


%% Display Results

% Dictionary Atoms Correlation
mAtomCorr = abs(mD.' * mD);

hFigure         = figure();
set(hFigure, 'Units', 'pixels');
set(hFigure, 'Position', [100, 100, 720, 700]);
hAxes           = axes();
hImageObject = imagesc(mAtomCorr);
hColorbarObject = colorbar(hAxes);
set(hAxes, 'Units', 'pixels');
set(hAxes, 'Position', [50, 50, 600, 600]);
set(get(hAxes, 'Title'), 'String', ['Dictionary Atoms Correlations'], ...
    'FontSize', titleFontSize);
set(get(hAxes, 'XLabel'), 'String', 'Atom Number', ...
    'FontSize', axisFotnSize);
set(get(hAxes, 'YLabel'), 'String', 'Atom Number', ...
    'FontSize', axisFotnSize);

% The Distribution of Atom Usage
vAtomUsage = sum((mA ~= 0), 2);

hFigure         = figure();
hAxes           = axes();
set(hFigure, 'Units', 'pixels');
set(hFigure, 'Position', [100, 100, 1200, 600]);
hBarSeries      = bar(vAtomUsage);
set(hBarSeries, 'BarWidth', 0.7);
set(get(hAxes, 'Title'), 'String', ['Atom Usage'], ...
    'FontSize', titleFontSize);
set(get(hAxes, 'XLabel'), 'String', 'Atom Number', ...
    'FontSize', axisFotnSize);
set(hAxes, 'XLim', [0.5, (numAtoms + 0.5)]);
set(get(hAxes, 'YLabel'), 'String', 'Frequency of Usage', ...
    'FontSize', axisFotnSize);


% The Reproduction of the Dictionary
% mBarColor = [0.0000, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];

mMaxCorrelation = [mMatchAtoms{1}(:, 3), mMatchAtoms{2}(:, 3), mMatchAtoms{3}(:, 3)];
titleString = {['Atom Reproduction']; ['Success Rate - ILS - ', num2str(round(NUM_TO_PCT_FACTOR * mean(mMatchAtoms{1}(:, 4)))), ...
    '%, K-SVD - ', num2str(round(NUM_TO_PCT_FACTOR * mean(mMatchAtoms{2}(:, 4)))), ...
    '%, RLS - ', num2str(round(NUM_TO_PCT_FACTOR * mean(mMatchAtoms{3}(:, 4)))), '%']};

vSuccessThr = simThr * ones([numAtoms, 1]);

hFigure         = figure();
hAxes           = axes();
set(hFigure, 'Units', 'pixels');
set(hFigure, 'Position', [100, 100, 900, 600]);
hLineSeries = plot([vSuccessThr, mMaxCorrelation]);
set(hLineSeries(1), 'LineStyle', ':', 'LineWidth', thinLineWidth);
set(hLineSeries(2:4), 'LineStyle', 'none');
set(hLineSeries(2), 'Marker', 'o');
set(hLineSeries(3), 'Marker', '+');
set(hLineSeries(4), 'Marker', '*');
set(hLineSeries(2), 'LineWidth', thinLineWidth);
set(hLineSeries(3), 'LineWidth', thinLineWidth);
set(hLineSeries(4), 'LineWidth', thinLineWidth);
set(hLineSeries, 'MarkerSize', mediumMarkerSize);
% hBarSeries      = bar(mMaxCorrelation);
% set(hBarSeries, 'BarWidth', 0.5);
% set(hBarSeries, 'EdgeColor', [0, 0, 0]);
% set(hBarSeries, 'LineWidth', 1);
% set(hBarSeries(1), 'FaceColor', mBarColor(1, :));
% set(hBarSeries(2), 'FaceColor', mBarColor(2, :));
% set(hBarSeries(3), 'FaceColor', mBarColor(3, :));
set(get(hAxes, 'Title'), 'String', titleString, ...
    'FontSize', titleFontSize);
set(get(hAxes, 'XLabel'), 'String', 'Atom Number', ...
    'FontSize', axisFotnSize);
set(hAxes, 'XLim', [0.5, (numAtoms + 0.5)]);
set(get(hAxes, 'YLabel'), 'String', 'Max Correlation', ...
    'FontSize', axisFotnSize);
hLegend = ClickableLegend({['Success Thr'], ['ILS'], ['K-SVD'], ['RLS']});
% hLegend = legend({['Success Thr'], ['ILS'], ['K-SVD'], ['RLS']});

% Dictionary Learning Error

vDlaErrRlsIdx                       = [1:size( vDlaErr{3}, 1)];
vDlaErrRlsIdx(isnan(vDlaErr{3}))    = [];
vDlaErr{3}(isnan(vDlaErr{3}))       = [];

hFigure         = figure();
hAxes           = axes();
set(hFigure, 'Units', 'pixels');
set(hFigure, 'Position', [100, 100, 800, 600]);
hLineSeries     = plot([1:size(vDlaErr{1}, 1)], log10(vDlaErr{1}), ...
    [1:size(vDlaErr{2}, 1)], log10(vDlaErr{2}), ...
    vDlaErrRlsIdx, log10(vDlaErr{3}));
set(hLineSeries, 'LineWidth', normalLineWidth);
set(get(hAxes, 'Title'), 'String', {['Dictionary Learning Error']; ['Normalized to Similar Run Time']}, ...
    'FontSize', titleFontSize);
set(get(hAxes, 'XLabel'), 'String', 'Iteration Number', ...
    'FontSize', axisFotnSize);
set(get(hAxes, 'YLabel'), 'String', 'Dictionery Learning Error log_10', ...
    'FontSize', axisFotnSize);
hLegend = ClickableLegend({['ILS'], ['K-SVD'], ['RLS']});
% hLegend = legend({['ILS'], ['K-SVD'], ['RLS']});


%% Restore Defaults
set(0, 'DefaultFigureWindowStyle', 'normal');
set(0, 'DefaultAxesLooseInset', defaultLooseInset);

