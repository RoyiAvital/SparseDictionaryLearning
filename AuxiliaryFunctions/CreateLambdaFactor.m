function [ vLambdaFactor ] = CreateLambdaFactor( initLambda, numIterations, numEffIter, pwrFctr, opMode )
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

OPERATION_MODE_POLYNOMIAL   = 1;
OPERATION_MODE_HYPERBOLA    = 2;
OPERATION_MODE_EXPONENTIAL  = 3;

switch(opMode)
    case(OPERATION_MODE_POLYNOMIAL)
        vLambdaFactor               = ones([numIterations, 1]);
        vLambdaFactor(1:numEffIter) = 1 - ((1 - initLambda) * ((1 - ([1:numEffIter].' / numEffIter)) .^ pwrFctr));
    case(OPERATION_MODE_HYPERBOLA)
        vLambdaFactor = 1 - ((1 - initLambda) ./ (1 + ([1:numIterations] / numEffIter)));
    case(OPERATION_MODE_EXPONENTIAL)
        vLambdaFactor = 1 - ((1 - initLambda) .* (0.5 .^ ([1:numIterations] / numEffIter)));
end




end