% Recursive Least Squares Dictionary Learning - LS and RLS Analysis
% TODO: Add a case fo Adaptive Weight (Slowly adapting at the begining,
% then staying constant).

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

smallSizeData   = 36;
mediumSizeData  = 48;
bigSizeData     = 60;

randomNumberStream = RandStream('mlfg6331_64', 'NormalTransform', 'Ziggurat');
subStreamNumber = 57162;
subStreamNumber = 2638;
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


%% Setting Analysis Parameters

modelOrder      = 3;
numMeasurements = 500;
vModel          = [0:(numMeasurements - 1)].';
mModel          = bsxfun(@power, vModel, [0:(modelOrder - 1)]);
vParams         = [0.7; 0.4; 0.1];

noiseStd    = 0.1;
noiseVar    = noiseStd * noiseStd;
vNoiseVar   = noiseVar * ones([numMeasurements, 1]);


%% Creating Data

vMeasurements = (mModel * vParams) + (noiseStd * randn([numMeasurements, 1]));


%% LS Estimation

[vEstParamsLs, mEstCovLs] = ApplyLs(vMeasurements, mModel, vNoiseVar);


%% RLS Estimation

vNoiseVar(:)    = 1;
dampingFactor   = 0.998;
vDampingFactor  = dampingFactor .^ [0:(numMeasurements - 1)];

% Initializtion of RLS using LS -> LS results
% Pay attention that for IID noise any conastant Noise Vector will have the
% same Estimated Parameters yet different Covariance.
[vEstParams, mEstCov] = ApplyLs(vMeasurements(1:3), mModel(1:3, :), vNoiseVar(1:3));
[vEstParamsRls, mEstCovRls] = ApplyRls(vMeasurements(4:end), mModel(4:end, :), vNoiseVar(4:end), vEstParams, mEstCov);

% Damping and Weighted are the same when the Weights are chosen right.
vEstParams  = zeros([modelOrder, 1]);
mEstCov     = noiseVar * eye(modelOrder);
% Exponentialy Damped - Wikipedia Style
[vEstParamsRls, mEstCovRls] = ApplyRlsDamped(vMeasurements, mModel, dampingFactor, vEstParams, mEstCov);
% Weighted LS Style
[vEstParamsRls, mEstCovRls] = ApplyRls(vMeasurements, mModel, vDampingFactor, vEstParams, mEstCov);


%% Comparison

numIterations = numMeasurements - modelOrder;

[vEstParamsInit01, mEstCovInit01] = ApplyLs(vMeasurements(1:modelOrder), mModel(1:modelOrder, :), vNoiseVar(1:modelOrder));
vNoiseVar01 = vNoiseVar;

vEstParamsInit02    = zeros([modelOrder, 1]);
mEstCovInit02       = noiseVar * eye(modelOrder);
vNoiseVar02         = vNoiseVar;

dampingFactor       = 0.98;
vEstParamsInit03    = zeros([modelOrder, 1]);
mEstCovInit03       = noiseVar * eye(modelOrder);
vNoiseVar03         = dampingFactor .^ [0:(numMeasurements - 1)];

vModelParamsRef     = zeros([numIterations, modelOrder]);
vModelParamsLs      = zeros([numIterations, modelOrder]);
vModelParamsRls01   = zeros([numIterations, modelOrder]);
vModelParamsRls02   = zeros([numIterations, modelOrder]);
vModelParamsRls03   = zeros([numIterations, modelOrder]);

rlsFirstMsmntIdx = modelOrder + 1;

for iIter = 1:numIterations
    measurementIdx = iIter + 3;
    
    vEstParamsLs    = ApplyLs(vMeasurements(1:measurementIdx), mModel(1:measurementIdx, :), vNoiseVar(1:measurementIdx));
    vEstParamsRls01 = ApplyRls(vMeasurements(rlsFirstMsmntIdx:measurementIdx), mModel(rlsFirstMsmntIdx:measurementIdx, :), vNoiseVar01(rlsFirstMsmntIdx:measurementIdx), vEstParamsInit01, mEstCovInit01);
    vEstParamsRls02 = ApplyRls(vMeasurements(rlsFirstMsmntIdx:measurementIdx), mModel(rlsFirstMsmntIdx:measurementIdx, :), vNoiseVar02(rlsFirstMsmntIdx:measurementIdx), vEstParamsInit02, mEstCovInit02);
    vEstParamsRls03 = ApplyRls(vMeasurements(rlsFirstMsmntIdx:measurementIdx), mModel(rlsFirstMsmntIdx:measurementIdx, :), vNoiseVar03(rlsFirstMsmntIdx:measurementIdx), vEstParamsInit03, mEstCovInit03);
    
    vModelParamsRef(iIter, :)   = vParams.';
    vModelParamsLs(iIter, :)    = vEstParamsLs;
    vModelParamsRls01(iIter, :) = vEstParamsRls01;
    vModelParamsRls02(iIter, :) = vEstParamsRls02;
    vModelParamsRls03(iIter, :) = vEstParamsRls03;
    
    disp(['Finished Iteration #', num2str(iIter, '%04d'), ' Out of ', num2str(numIterations), ' Iterations']);
end


%% Display Results

hFigure         = figure('Units', 'pixels', 'Position', [100, 100, 1200, 800]);

axesIdx         = 1;
hAxes(axesIdx)  = subplot(3, 1, axesIdx);
hLineSeries     = plot(1:numIterations, [vModelParamsRef(:, axesIdx), vModelParamsLs(:, axesIdx), vModelParamsRls01(:, axesIdx), vModelParamsRls02(:, axesIdx), vModelParamsRls03(:, axesIdx)]);
set(get(hAxes(axesIdx), 'Title'), 'String', ['Polynomial Parameters Estimation - Parameter #', num2str(axesIdx, '%03d')], ...
    'FontSize', axisFotnSize);
set(hLineSeries(1), 'LineWidth', thinLineWidth, 'LineStyle', ':');
set(hLineSeries(2:5), 'LineWidth', normalLineWidth);
set(get(hAxes(axesIdx), 'XLabel'), 'String', 'Iteration Number', ...
    'FontSize', axisFotnSize);
set(get(hAxes(axesIdx), 'YLabel'), 'String', 'Parameter Value', ...
    'FontSize', axisFotnSize);
hLegend = legend({['Ground Truth'], ['Batch LS'], ['RLS - LS Init'], ['RLS - Zero Init'], ['RLS - Zero Init, Damped']});

axesIdx         = 2;
hAxes(axesIdx)  = subplot(3, 1, axesIdx);
hLineSeries     = plot(1:numIterations, [vModelParamsRef(:, axesIdx), vModelParamsLs(:, axesIdx), vModelParamsRls01(:, axesIdx), vModelParamsRls02(:, axesIdx), vModelParamsRls03(:, axesIdx)]);
set(get(hAxes(axesIdx), 'Title'), 'String', ['Polynomial Parameters Estimation - Parameter #', num2str(axesIdx, '%03d')], ...
    'FontSize', axisFotnSize);
set(hLineSeries(1), 'LineWidth', thinLineWidth, 'LineStyle', ':');
set(hLineSeries(2:5), 'LineWidth', normalLineWidth);
set(get(hAxes(axesIdx), 'XLabel'), 'String', 'Iteration Number', ...
    'FontSize', axisFotnSize);
set(get(hAxes(axesIdx), 'YLabel'), 'String', 'Parameter Value', ...
    'FontSize', axisFotnSize);
hLegend = legend({['Ground Truth'], ['Batch LS'], ['RLS - LS Init'], ['RLS - Zero Init'], ['RLS - Zero Init, Damped']});

axesIdx         = 3;
hAxes(axesIdx)  = subplot(3, 1, axesIdx);
hLineSeries     = plot(1:numIterations, [vModelParamsRef(:, axesIdx), vModelParamsLs(:, axesIdx), vModelParamsRls01(:, axesIdx), vModelParamsRls02(:, axesIdx), vModelParamsRls03(:, axesIdx)]);
set(get(hAxes(axesIdx), 'Title'), 'String', ['Polynomial Parameters Estimation - Parameter #', num2str(axesIdx, '%03d')], ...
    'FontSize', axisFotnSize);
set(hLineSeries(1), 'LineWidth', thinLineWidth, 'LineStyle', ':');
set(hLineSeries(2:5), 'LineWidth', normalLineWidth);
set(get(hAxes(axesIdx), 'XLabel'), 'String', 'Iteration Number', ...
    'FontSize', axisFotnSize);
set(get(hAxes(axesIdx), 'YLabel'), 'String', 'Parameter Value', ...
    'FontSize', axisFotnSize);
hLegend = legend({['Ground Truth'], ['Batch LS'], ['RLS - LS Init'], ['RLS - Zero Init'], ['RLS - Zero Init, Damped']});

hSupTitle = suptitle(['Comparison Between Batch LS and RLS']);
set(hSupTitle, 'FontSize', titleFontSize);

% Print Figure to EPS format
% print(hFigure, 'Figures\LSRLSAnalysis', '-depsc');


%% Restore Defaults
% set(0, 'DefaultFigureWindowStyle', 'normal');
set(0, 'DefaultAxesLooseInset', defaultLooseInset);

