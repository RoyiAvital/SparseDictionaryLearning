% Recursive Least Squares Dictionary Learning - Parameters Analysis
% TODO: Add a case fo Adaptive Weight (Slowly adapting at the begining,
% then staying constant).

%% General Parameters and Initialization

clear();
% close('all');

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

simString = 'RLS DLA - Parameter Analysis';

% Parameters for Data Generation
numAtoms        = 128;
atomLength      = 64;
dataCardinality = 4;
numSamples      = 8192;
noiseStd        = 0.0;

% Parameters for Analysis
dataCardinalityDla  = 4;
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

vRlsDlaConfig   = [RLS_CONFIG_CONSTANT, RLS_CONFIG_LINEAR_100, RLS_CONFIG_QUADRATIC_100, RLS_CONFIG_CUBIC_100];
rlsConfigString = {['Constant'], ['Linear'], ['Quadratic'], ['Cubic']};
numConfig       = length(vRlsDlaConfig);

rldDlaRepFctr   = 4;
numIterations   = rldDlaRepFctr * numSamples;
mLambdaFctr     = zeros([numIterations, numConfig]);

for ii = 1:numConfig
    % Lambda:
    rlsDlaConfig = vRlsDlaConfig(ii);
    [lambdaFctrInitialValue, lambdaFctrEffIterFactor, lambdaFctrPwrFctr, lambdaFctrMode] = GetLambdaParams(rlsDlaConfig);
    numEffIter    = round(lambdaFctrEffIterFactor * numIterations);
    mLambdaFctr(:, ii)   = CreateLambdaFactor(0.95, numIterations, numEffIter, lambdaFctrPwrFctr, lambdaFctrMode);
end

mLambdaFctr(:, 1)   = CreateLambdaFactor(1, numIterations, numEffIter, 0, 1);


%% Generating Data

% Reference Dictionary
mD      = CreateRandomDictionary(atomLength, numAtoms);
% Dictionary for Initialization
mDInit  = CreateRandomDictionary(atomLength, numAtoms);

% Create Data - mX = mD * mA
[mX, mA] = CreateDataSamples(mD, dataCardinality, numSamples, noiseStd);

methodString = 'RLS DLA';


%% Dictionary Analysis

mDlaErr         = zeros([(numIterations + 1), numConfig]);
mMaxCorr        = zeros([numAtoms, numConfig]);
mMatchedAtoms   = zeros([numAtoms, numConfig]);

for ii = [1:size(mLambdaFctr, 2)]
    
    disp([' ']);
    disp([simString, ' - ', methodString]);
    disp([' ']);
    
    vLambdaFctr = mLambdaFctr(:, ii);

    hTicTimer               = tic();
    [mTrainD, mTrainA, vDlaErr]  = ApplyRlsDla(repmat(mX, [1, rldDlaRepFctr]), mDInit, vLambdaFctr, ompThr, ompMode, dicErrThr);
    runTime = toc(hTicTimer);
    
    disp([' ']); disp([simString, ' - ', methodString, ' - Run Time - ', num2str(runTime), ' [Sec]']);
    disp([' ']);
    
    mDlaErr(:, ii) = vDlaErr;

    %-- Evaluating Results
    mMatchAtoms             = AnalyzeDictionarySimilarity(mD, mTrainD, simThr);
    mMaxCorr(:, ii)         = mMatchAtoms(:, 3);
    mMatchedAtoms(:, ii)    = mMatchAtoms(:, 4);
end


%% Display Results

% The Reproduction of the Dictionary
% mBarColor = [0.0000, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];

vSuccessThr = simThr * ones([numAtoms, 1]);

hFigure         = figure();
set(hFigure, 'Units', 'pixels');
set(hFigure, 'Position', [100, 100, 800, 600]);

hAxes           = subplot(3, 1, 1);
set(hFigure, 'Units', 'pixels');
set(hFigure, 'Position', [100, 100, 900, 900]);
hLineSeries = plot(mLambdaFctr);
set(hLineSeries,'LineWidth', normalLineWidth);
set(get(hAxes, 'Title'), 'String', 'Lambda Factor Per Iteration', ...
    'FontSize', titleFontSize);
set(get(hAxes, 'XLabel'), 'String', 'Iteration Number', ...
    'FontSize', axisFotnSize);
set(get(hAxes, 'YLabel'), 'String', '\lambda', ...
    'FontSize', axisFotnSize);
hLegend = ClickableLegend(rlsConfigString);


hAxes           = subplot(3, 1, 2);

titleString = {['Atom Reproduction']; ['Success Rate - Constant - ', num2str(round(NUM_TO_PCT_FACTOR * mean(mMatchedAtoms(:, 1)))), ...
    '%, Linear - ', num2str(round(NUM_TO_PCT_FACTOR * mean(mMatchedAtoms(:, 2)))), ...
    '%, Quadratic - ', num2str(round(NUM_TO_PCT_FACTOR * mean(mMatchedAtoms(:, 3)))), ...
    '%, Cubic - ', num2str(round(NUM_TO_PCT_FACTOR * mean(mMatchedAtoms(:, 4)))), '%']};

hLineSeries = plot([vSuccessThr, mMaxCorr]);
set(hLineSeries(1), 'LineStyle', ':', 'LineWidth', thinLineWidth);
set(hLineSeries(2:end), 'LineStyle', 'none');
set(hLineSeries(2), 'Marker', 'o');
set(hLineSeries(3), 'Marker', '+');
set(hLineSeries(4), 'Marker', '*');
set(hLineSeries(5), 'Marker', 'x');
set(hLineSeries(2:end), 'LineWidth', thinLineWidth);
set(hLineSeries(2:end), 'MarkerSize', mediumMarkerSize);
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
hLegend = ClickableLegend({['Success Thr'], rlsConfigString{:}});

% Dictionary Learning Error
hAxes           = subplot(3, 1, 3);

vDlaErrRlsIdx  = [1:size(mDlaErr, 1)];
vRowNanIdx = false([size(mDlaErr, 1), 1]);

for ii = 1:length(vDlaErrRlsIdx)
    if(any(isnan(mDlaErr(ii, :))))
        vRowNanIdx(ii) = true();
    end
end

vDlaErrRlsIdx(vRowNanIdx)   = [];
mDlaErr(vRowNanIdx, :)      = [];

hLineSeries     = plot(vDlaErrRlsIdx, log10(mDlaErr));
set(hLineSeries, 'LineWidth', normalLineWidth);
set(get(hAxes, 'Title'), 'String', {['Dictionary Learning Error']}, ...
    'FontSize', titleFontSize);
set(get(hAxes, 'XLabel'), 'String', 'Iteration Number', ...
    'FontSize', axisFotnSize);
set(get(hAxes, 'YLabel'), 'String', 'Dictionery Learning Error log_10', ...
    'FontSize', axisFotnSize);
hLegend = ClickableLegend(rlsConfigString);


%% Restore Defaults
% set(0, 'DefaultFigureWindowStyle', 'normal');
set(0, 'DefaultAxesLooseInset', defaultLooseInset);

