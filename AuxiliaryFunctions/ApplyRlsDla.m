function [ mD, mA, vDicErr ] = ApplyRlsDla( mX, mDInit, vLambdaFctr, ompThr, ompMode, dicErrThr )
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
%   3.  Pre Calculate the Gram Matrix for the OMP (See 2.4.1 on the6
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

OMP_MODE_CARDINALITY    = 1;
OMP_MODE_NOISE_THR      = 2;

DEBUG_MODE_GRANULARITY_FACTOR = 100;


%% Initialize Parameters

debugMode = ON;

numAtoms   = size(mDInit, 2);
numSamples = size(mX, 2);

switch(ompMode)
    case(OMP_MODE_CARDINALITY)
        hBatchOmp   = @(mD)     omp( (mD.' * mX), (mD.' * mD), ompThr);
        hSingleOmp  = @(mD, vX) omp(mD, vX, [], ompThr);
    case(OMP_MODE_NOISE_THR)
        vNormX      = sum(mX .* mX);
        hBatchOmp   = @(mD)     omp2((mD.' * mX), vNormX, (mD.' * mD), ompThr);
        hSingleOmp  = @(mD, vX) omp2(mD, vX, [], ompThr);
end

vDicErr = nan([(numSamples + 1), 1]);
    

%% Running Algorithm

% Updtae Dictionary
mD = mDInit;
% Update Eeights
mA = hBatchOmp(mD);
mC = eye(numAtoms);

mE          = (mD * mA) - mX;
vDicErr(1)  = mean(mE(:) .^ 2);

if(debugMode == ON)
    disp([' ']);
    disp(['Initial Dictionary Error - ', num2str(vDicErr(1))]);
end

for ii = 1 : numSamples
    vX = mX(:, ii);
    vA = hSingleOmp(mD, vX);
    % vA = OrthogonalMatchingPursuit(vX, mD, 1000, ompThr);
    
    % Representation Error
    vR = vX - (mD * vA);
    
    mC_star = mC / vLambdaFctr(ii);
    vU      = mC_star * vA;
    % vV = mDictionary.' * vRepErr;
    
    alphaFactor = 1 / (1 + (vA.' * vU));
    
    
    % Update the Dictionary
    % Though it is faster to first multiply the vector with the factor, it
    % is done to keep symmetry.
    mC = mC_star - (alphaFactor * (vU * vU.'));
    mD = mD + (alphaFactor * (vR * vU.'));
    
    mD = bsxfun(@rdivide, mD, sqrt(sum((mD .* mD), 1)));
    
    if(mod(ii, DEBUG_MODE_GRANULARITY_FACTOR) == 0)
        mA                  = hBatchOmp(mD);
        mE                  = (mD * mA) - mX;
        vDicErr(ii + 1)     = mean(mE(:) .^ 2);
    end
    
    % Error
    if(debugMode == ON)
        if(mod(ii, DEBUG_MODE_GRANULARITY_FACTOR) == 0)
            disp(['Dictionary Error at Iteration -  ', num2str(ii), ' Out of - ', num2str(numSamples), ' is - ', num2str(vDicErr(ii + 1))]);
        end
    end
end

if(debugMode == ON)    
    disp(['Final Dictionary Error - ', num2str(vDicErr(ii + 1))]);
    disp([' ']);
end


end

