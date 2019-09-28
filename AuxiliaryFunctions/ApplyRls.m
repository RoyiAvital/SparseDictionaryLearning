function [ vEstParams, mEstCov ] = ApplyRls( vMeasurements, mModel, vNoiseVar, vEstParams, mEstCov )
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

numIterations   = size(vMeasurements, 1);
modelOrder      = size(mModel, 2);
mI              = eye(modelOrder);

for iIter = 1:numIterations
    vModelRow   = mModel(iIter, :);
    vGainFactor = (mEstCov * vModelRow.') / (vNoiseVar(iIter) + (vModelRow * mEstCov * vModelRow.'));
    vEstParams  = vEstParams + (vGainFactor * (vMeasurements(iIter) - (vModelRow * vEstParams)));
    mEstCov     = (mI - (vGainFactor * vModelRow)) * mEstCov;
end


end

