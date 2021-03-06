function [ mMatchAtoms ] = AnalyzeDictionarySimilarity( mD1, mD2, simThr )
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

NUM_FEATURES = 4;

ATOM_IDX_DIC_1_IDX  = 1;
ATOM_IDX_DIC_2_IDX  = 2;
ATOM_CORR_IDX       = 3;
ATOM_SIM_FLAG       = 4;

numAtoms = size(mD1, 2);

% Coherence
mMutualCoherence = mD2.' * mD1;

[vMaxCoherence, vMaxCoherenceIdx] = max(abs(mMutualCoherence));

mMatchAtoms = zeros(numAtoms, NUM_FEATURES);

for ii = 1:numAtoms
    mMatchAtoms(ii, ATOM_IDX_DIC_1_IDX) = ii;
    mMatchAtoms(ii, ATOM_IDX_DIC_2_IDX) = vMaxCoherenceIdx(ii);
    mMatchAtoms(ii, ATOM_CORR_IDX)      = vMaxCoherence(ii);
    mMatchAtoms(ii, ATOM_SIM_FLAG)      = (vMaxCoherence(ii) >= simThr);
end


end