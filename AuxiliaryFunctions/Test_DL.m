
close all;
clear;

% addpath('../ompbox10/');
addpath( genpath('../') );

numberOfAtom = 300;
atomLength   = 10;
cardinality  = 3;

%% Generate D:
mD = randn(atomLength, numberOfAtom);
mD = bsxfun(@rdivide, mD, sqrt(sum(mD.^2, 1)));

% mX = CreateSamples(mD, cardinality, 5000, 0);
mX = CreateDataSamples(mD, cardinality, 1000, 0);

%% DL:
mD_init = randn(atomLength, numberOfAtom);
mD_init = bsxfun(@rdivide, mD_init, sqrt(sum(mD_init.^2, 1)));

% %-- MOD:
% mD_est  = DL_MOD(mD_init, mX, cardinality);
% mA      = omp(mD_est' * mX, mD_est' * mD_est, cardinality);
% mX_est  = mD_est * mA;
% A       = mean(abs(mX(:) - mX_est(:)))
% 
% %-- KSVD:
% mD_est = DL_KSVD(mD_init, mX, cardinality);
% mA     = omp(mD_est' * mX, mD_est' * mD_est, cardinality);
% mX_est = mD_est * mA;
% B      = mean(abs(mX(:) - mX_est(:)))

%-- KSVD:
mD_est = DL_RLS(mD_init, mX, cardinality, ones(1, size(mX,2)));
mA     = omp(mD_est' * mX, mD_est' * mD_est, cardinality);
mX_est = mD_est * mA;
C      = mean(abs(mX(:) - mX_est(:)))

