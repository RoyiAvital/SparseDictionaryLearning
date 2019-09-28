% Recursive Least Squares Dictionary Learning - Image Denoising Analysis

%% General Parameters and Initialization

clear();
close('all');

% set(0, 'DefaultFigureWindowStyle', 'docked');
% defaultLooseInset = get(0, 'DefaultAxesLooseInset');
% set(0, 'DefaultAxesLooseInset', [0.05, 0.05, 0.05, 0.05]);

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
addpath(genpath('../AuxiliaryFunctions'));


%% Setting Constants

FALSE   = 0;
TRUE    = 1;

OFF = 0;
ON  = 1;

OPERATION_MODE_POLYNOMIAL   = 1;
OPERATION_MODE_HYPERBOLA    = 2;
OPERATION_MODE_EXPONENTIAL  = 3;

OMP_MODE_CARDINALITY    = 1;
OMP_MODE_NOISE_THR      = 2;

LAMBDA_FACTOR_MODE_POLYNOMIAL   = 1;
LAMBDA_FACTOR_MODE_HYPERBOLA    = 2;
LAMBDA_FACTOR_MODE_EXPONENTIAL  = 3;

RLS_CONFIG_CONSTANT         = 1;
RLS_CONFIG_LINEAR_100       = 2;
RLS_CONFIG_LINEAR_085       = 3;
RLS_CONFIG_LINEAR_075       = 4;
RLS_CONFIG_QUADRATIC_100    = 5;
RLS_CONFIG_QUADRATIC_085    = 6;
RLS_CONFIG_QUADRATIC_075    = 7;
RLS_CONFIG_CUBIC_100        = 8;
RLS_CONFIG_CUBIC_085        = 9;
RLS_CONFIG_CUBIC_075        = 10;
RLS_CONFIG_HYPERBOLA        = 11;
RLS_CONFIG_EXP              = 12;

%% Setting Analysis Parameters
patch_width        = 8;
patch_height       = 8;
vPatch_size        = [patch_height, patch_width];
atomLength         = patch_width * patch_height;
numAtoms           = 300;
dataCardinalityDla = 5;
training_set_size  = 5000;

%%
noiseStdThr = 0;
ompMode     = OMP_MODE_CARDINALITY;
dicErrThr   = 0;

switch(ompMode)
    case(OMP_MODE_CARDINALITY)
        ompThr       = dataCardinalityDla;
    case(OMP_MODE_NOISE_THR)
        ompThr       = noiseStdThr;
end

%% Loading Data
mInputImage = imread('barbara.png');

mInputImage = mean(double(mInputImage), 3);
mInputImage = mInputImage(301:500, 301:500);
figure; imshow(mInputImage, []);

mSuper_set = im2col(mInputImage, vPatch_size);

%-- Remove mean:
vSuper_set_mean = mean(mSuper_set, 1);
mSuper_set      = bsxfun(@minus, mSuper_set, vSuper_set_mean);

vTrain_set_idx  = randperm(length(mSuper_set), training_set_size);
mTrain          = mSuper_set(:, vTrain_set_idx);

mDInit          = CreateRandomDictionary(atomLength, numAtoms);

%% Lambda:
rlsDlaConfig = RLS_CONFIG_LINEAR_100;
[lambdaFctrInitialValue, lambdaFctrEffIterFactor, lambdaFctrPwrFctr, lambdaFctrMode] = ...
                                                                        GetLambdaParams(rlsDlaConfig);
numIterations = training_set_size;
numEffIter    = round(lambdaFctrEffIterFactor * numIterations);
vLambda       = CreateLambdaFactor(lambdaFctrInitialValue, numIterations, numEffIter, lambdaFctrPwrFctr, lambdaFctrMode);

if 1
cMethodString = {'ILS DLA', 'K-SVD DLA', 'RLS DLA'};
cMethodDla    = {@ApplyIlsDla, @ApplyKSvdDla, @ApplyRlsDla};

figure;
for ii = 1 : length(cMethodDla)
    methodString = cMethodString{ii};
    disp([' ']); disp(methodString);

    hTicTimer               = tic();
    Dla                     = cMethodDla{ii};
    switch (methodString)
        case {'ILS DLA', 'K-SVD DLA'}, cParams = {mTrain, mDInit, ompThr, ompMode, dicErrThr};
        case 'RLS DLA',                cParams = {mTrain, mDInit, vLambda, ompThr, ompMode, dicErrThr};
    end
    [mD, mTrainA, vDlaErr]  = Dla(cParams{:});
    runTime                 = toc(hTicTimer);
    
    disp([' ']); disp([methodString, ' - Run Time - ', num2str(runTime), ' [Sec]']);

    %-- Evaluating Results
    vW = sqrt( sum(mD.^2, 1) );
    mD = bsxfun(@rdivide, mD, vW);
    mG = mD' * mD;
    mA = omp(mD' * mSuper_set, mG, dataCardinalityDla);
    
    vActivity  = sum(abs(mA), 2);
    mR         = mD * mA - mSuper_set;
    RMSE_batch = sqrt( mean( mR(:).^2 ) );
    
    mP = bsxfun(@plus, mD * mA, vSuper_set_mean);
    mRecImage = Col_To_Im(mP, size(mInputImage), vPatch_size);
    
    
%     subplot(3,3,3*ii-2); imshow([mInputImage, mRecImage],[]); title(['RMSE = ', num2str(RMSE_batch)]);
%     subplot(3,3,3*ii-1); stem(vActivity);                     title('Activity');
%     subplot(3,3,3*ii-0); Display_D(mD);                       title([cMethodString, ' Dic.']);
    
    subplot(3,2,2*ii-1); imshow([mInputImage, mRecImage],[]); title(['RMSE = ', num2str(RMSE_batch)]);
    subplot(3,2,2*ii-0); Display_D(mD);                       title([cMethodString, ' Dic.']);
end

%% Image Denoising - ILS-DLA
[mDIls, mAIls]  = ApplyIlsDla(mTrain, mDInit, ompThr, ompMode, dicErrThr);

%% Image Denoising - K-SVD-DLA
[mDKsvd, mAKsvd] = ApplyKSvdDla(mTrain, mDInit, ompThr, ompMode, dicErrThr);

%% Image Denoising - RLS-DLA
vLambda = CreateLambdaFactor(initLambda, numIterations, numEffIter, pwrFctr, opMode);
[mDRls, mARls]   = ApplyRlsDla(mTrain, mDInit, vLambda, ompThr, ompMode, dicErrThr);

%% 
figure;
cD = {mDIls, mDKsvd, mDRls};
% cA = {mAIls, mAKsvd, mARls};
for ii = 1 : length(cD)
    mD = cD{ii};
    
    vW = sqrt( sum(mD.^2, 1) );
    mD = bsxfun(@rdivide, mD, vW);
    mG = mD' * mD;
    mA = omp(mD' * mSuper_set, mG, dataCardinalityDla);
    
    vActivity  = sum(abs(mA), 2);
    mR         = mD * mA - mSuper_set;
    RMSE_batch = sqrt( mean( mR(:).^2 ) );
    
    mP = bsxfun(@plus, mD * mA, vSuper_set_mean);
    mRecImage = Col_To_Im(mP, size(mInputImage), vPatch_size);
    
    subplot(3,3,3*ii-2); imshow([mInputImage, mRecImage],[]);
    subplot(3,3,3*ii-1); stem(vActivity); title('Activity');
    subplot(3,3,3*ii-0); Display_D(mD);   title(['RMSE = ', num2str(RMSE_batch)]);
end
end

%%
noise_std = 10;
I_noisy   = mInputImage + noise_std * randn(size(mInputImage));

%-- Create Super Set from the Image:
mSuper_set = im2col(I_noisy, vPatch_size);
%-- Remove mean:
vSuper_set_mean = mean(mSuper_set, 1);
mSuper_set      = bsxfun(@minus, mSuper_set, vSuper_set_mean);

vTrain_set_idx    = randperm(length(mSuper_set), training_set_size);
mTrain            = mSuper_set(:, vTrain_set_idx);


ompThr    = sqrt(1.15 * atomLength * noise_std^2);
ompMode   = OMP_MODE_NOISE_THR;
dicErrThr = 0;


if 0
%% Image Denoising - ILS-DLA
[mDIls, mAIls] = ApplyIlsDla(mTrain, mDInit, ompThr, ompMode, dicErrThr);

%% Image Denoising - K-SVD-DLA
[mDKsvd, mAKsvd] = ApplyKSvdDla(mTrain, mDInit, ompThr, ompMode, dicErrThr);

%% Image Denoising - RLS-DLA
initLambda    = 0.9;
numIterations = training_set_size;
numEffIter    = training_set_size / 2;
pwrFctr       = 2;
opMode        = OPERATION_MODE_POLYNOMIAL;

vLambda = CreateLambdaFactor(initLambda, numIterations, numEffIter, pwrFctr, opMode);
vLambda = ones(1, numIterations);
[mDRls, mARls] = ApplyRlsDla(mTrain, mDInit, vLambda, ompThr, ompMode, dicErrThr);

%% 

cD = {mDIls, mDKsvd, mDRls};
% cA = {mAIls, mAKsvd, mARls};
for ii = 1 : length(cD)
    figure;
    mD = cD{ii};
    
    vW = sqrt( sum(mD.^2, 1) );
    mD = bsxfun(@rdivide, mD, vW);
    mG = mD' * mD;
    mA = omp2((mD.' * mSuper_set), sum(mSuper_set .* mSuper_set), (mD.' * mD), ompThr);
%     mA = omp(mD' * mSuper_set, mG, dataCardinalityDla);
    vActivity  = sum(abs(mA), 2);
    
    mR         = mD * mA - mSuper_set;
    RMSE_batch = sqrt( mean( mR(:).^2 ) );
    
    mP = bsxfun(@plus, mD * mA, vSuper_set_mean);
    mRecImage = Col_To_Im(mP, size(mInputImage), vPatch_size);
    
    subplot(1,3,1); imshow([I_noisy, mRecImage],[]);
    subplot(1,3,2); stem(vActivity); title('Activity');
    subplot(1,3,3); Display_D(mD);   title(['RMSE = ', num2str(RMSE_batch)]);
end
end

%% Display Results

% hFigure    = figure();
% hAxes      = axes();
% hBarSeries = bar(vLsRmsError);
% set(get(hAxes, 'Title'), 'String', ['Category - RMS Error'], ...
%     'FontSize', titleFontSize);
% set(get(hAxes, 'XLabel'), 'String', 'Category', ...
%     'FontSize', axisFotnSize);
% set(hAxes, 'XTick', [1:(numCategories + 1)]);
% set(hAxes, 'XTickLabel', xTickLabelStr);
% set(hAxes, 'XLim', [0.5, (numCategories + 1.5)]);
% set(get(hAxes, 'YLabel'), 'String', 'Error [$]', ...
%     'FontSize', axisFotnSize);


%% Restore Defaults
% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLooseInset);

