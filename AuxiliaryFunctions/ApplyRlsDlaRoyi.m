function [ mD, mA ] = ApplyRlsDla( mX, mDInit, ompThr, ompMode, dicErrThr )
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

FALSE   = 0;
TRUE    = 1;

OFF = 0;
ON  = 1;

OMP_MODE_CARDINALITY    = 1;
OMP_MODE_NOISE_THR      = 2;

debugMode                   = ON;
debugModeGranularityFactor  = 20;

numAtoms   = size(mDInit, 2);
numSamples = size(mX, 2);

% Initialization
% Updtae Dictionary
mD = mDInit;
% Update Eeights
mA = omp((mD.' * mX), (mD.' * mD), ompThr);
mC = eye(numAtoms);

if(debugMode == ON)
    mE      = (mD * mA) - mX;
    dicErr  = mean(mE(:) .^ 2);
    
    disp(['Initial Dictionary Error - ', num2str(dicErr)]);
end

for ii = 1:numSamples
    vX = mX(:, ii);
    vA = omp(mD, vX, [], ompThr);
    
    % Representation Error
    vR = vX - (mD * vA);
    
    vU = mC * vA;
    % vV = mDictionary.' * vRepErr;
    
    alphaFactor = 1 / (1 + (vA.' * vU));
    
    
    % Update the Dictionary
    % Though it is faster to first multiply the vector with the factor, it
    % is done to keep symmetry.
    mC = mC - (alphaFactor * (vU * vU.'));
    mD = mD + (alphaFactor * (vR * vU.'));
    
    mD = bsxfun(@rdivide, mD, sqrt(sum((mD .* mD), 1)));
    
    % Error
    if(debugMode == ON)
        if(mod(ii, debugModeGranularityFactor) == 0)
            mA      = omp((mD.' * mX), (mD.' * mD), ompThr);
            mE      = (mD * mA) - mX;
            dicErr  = mean(mE(:) .^ 2);
            
            disp(['Dictionary Error at Iteration -  ', num2str(ii), ' Out of - ', num2str(numSamples), ' is - ', num2str(dicErr)]);
        end
    end
end

mA = omp((mD.' * mX), (mD.' * mD), ompThr);

if(debugMode == ON)
    mE      = (mD * mA) - mX;
    dicErr  = mean(mE(:) .^ 2);
    
    disp(['Final Dictionary Error - ', num2str(dicErr)]);
end


end

