function [ mD, mA, vDicErr ] = ApplyIlsDla( mX, mDInit, ompThr, ompMode, dicErrThr )
% ----------------------------------------------------------------------------------------------- %
%[ mDictionary, mWeights ] = ApplyRlsDla( mInputData )
% Applies Dictionary Learning using Recursive LEast Squares Method.
% Input:
%   - mInputImage       -   Input Image.
%                           Structure: Image Matrix (1 / 3 Channels).
%                           Type: 'Single' / 'Double'.
%                           Range: [0, 1].
%   - spatialRadius     -   Spatial Radius.
%                           The Spatial Radius of the Guided Filter.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range: {1, 2, 3, ...}.
%   - rangeRadius       -   Range Radius.
%                           Sets the sensitivity of the Guided Filter to
%                           Tonal (Range) Differences.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range {1, 2, 3, ...}.
%   - luminosityMode    -   Lumninosity Mode.
%                           Lumninosity Mode Binary Flag.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range {0, 1}.
% Output:
%   - mOutputImage      -   Output Image.
%                           Structure: Image Matrix (1 / 3 Channels).
%                           Type: 'Single' / 'Double'.
%                           Range: [0, 1].

% References
%   1.  Recursive Least Squares Dictionary Learning Algorithm.
% Remarks:
%   1.  Prefixes:
%       -   'm' - Matrix.
%       -   'v' - Vector.
%   2.  Colofd
% TODO:
%   1.  Add "Forgetting Factor".
%   2.  Use ORMP instead of OMP.
%   3.  Pre Calculate the Gram Matrix for the OMP (See 2.4.1 on the
%       article).
% Release Notes:
%   -   1.0.000    15/02/2016
%       *   First realease version.
% ----------------------------------------------------------------------------------------------- %

%% Initialize Constants

FALSE   = 0;
TRUE    = 1;

OFF = 0;
ON  = 1;

OMP_MODE_CARDINALITY  = 1;
OMP_MODE_NOISE_THR    = 2;

NUM_ITERATIONS_FACTOR = 1 / 10;

DEBUG_MODE_GRANULARITY_FACTOR = 10;


%% Initialize Parameters

debugMode = ON;

atomLength    = size(mDInit, 1);
numSamples    = size(mX, 2);
maxIterations = ceil(numSamples * NUM_ITERATIONS_FACTOR);

switch(ompMode)
    case(OMP_MODE_CARDINALITY)
        hBatchOmp = @(mD) omp((mD.' * mX), (mD.' * mD), ompThr);
    case(OMP_MODE_NOISE_THR)
        vNormX = sum(mX .* mX);
        hBatchOmp = @(mD) omp2((mD.' * mX), vNormX, (mD.' * mD), ompThr);
end

vDicErr = nan([(maxIterations + 1), 1]);


%% Running Algorithm

% Initialization
% Updtae Dictionary
mD = mDInit;
% Update Weights
mA = hBatchOmp(mD);

mE          = (mD * mA) - mX;
vDicErr(1)  = mean(mE(:) .^ 2);

if(debugMode == ON)
    disp([' ']);
    disp(['Initial Dictionary Error - ', num2str(vDicErr(1))]);
end

for ii = 1:maxIterations
    % Update Dictionary
    mD = mX * mA.' / (mA * mA.');
    mD = bsxfun(@rdivide, mD, sqrt(sum((mD .* mD), 1)));
    
    % Update Weights
    mA = hBatchOmp(mD);
    
    % Error
    mE          = (mD * mA) - mX;
    repError    = (norm(mE, 'fro') ^ 2);
    vDicErr(ii + 1) = repError / (atomLength * numSamples); %<! mean(mE(:) .^ 2);
    if(repError <= (numSamples * dicErrThr))
        break;
    end
    
    if(debugMode == ON)
        if(mod(ii, DEBUG_MODE_GRANULARITY_FACTOR) == 0)
            disp(['Dictionary Error at Iteration -  ', num2str(ii), ' Out of - ', num2str(maxIterations), ' is - ', num2str(vDicErr(ii + 1))]);
        end
    end
end

if(debugMode == ON)    
    disp(['Final Dictionary Error - ', num2str(vDicErr(ii + 1))]);
    disp([' ']);
end


end

