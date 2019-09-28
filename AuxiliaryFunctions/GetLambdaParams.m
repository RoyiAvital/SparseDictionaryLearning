function [lambdaFctrInitialValue, lambdaFctrEffIterFactor, lambdaFctrPwrFctr, lambdaFctrMode] = ...
                                GetLambdaParams(rlsDlaConfig)
                            
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

LAMBDA_FACTOR_MODE_POLYNOMIAL   = 1;
LAMBDA_FACTOR_MODE_HYPERBOLA    = 2;
LAMBDA_FACTOR_MODE_EXPONENTIAL  = 3;

switch(rlsDlaConfig)
    case(RLS_CONFIG_CONSTANT)
        lambdaFctrInitialValue  = 1;
        lambdaFctrEffIterFactor = 1;
        lambdaFctrPwrFctr       = 0;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_POLYNOMIAL;
    case(RLS_CONFIG_LINEAR_100)
        lambdaFctrInitialValue  = 0.95;
        lambdaFctrEffIterFactor = 1.00;
        lambdaFctrPwrFctr       = 1;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_POLYNOMIAL;
    case(RLS_CONFIG_LINEAR_085)
        lambdaFctrInitialValue  = 0.85;
        lambdaFctrEffIterFactor = 0.85;
        lambdaFctrPwrFctr       = 1;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_POLYNOMIAL;
    case(RLS_CONFIG_LINEAR_075)
        lambdaFctrInitialValue  = 0.98;
        lambdaFctrEffIterFactor = 0.75;
        lambdaFctrPwrFctr       = 1;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_POLYNOMIAL;
    case(RLS_CONFIG_LINEAR_050)
        lambdaFctrInitialValue  = 0.99;
        lambdaFctrEffIterFactor = 0.50;
        lambdaFctrPwrFctr       = 1;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_POLYNOMIAL;
    case(RLS_CONFIG_QUADRATIC_100)
        lambdaFctrInitialValue  = 0.95;
        lambdaFctrEffIterFactor = 1.00;
        lambdaFctrPwrFctr       = 2;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_POLYNOMIAL;
    case(RLS_CONFIG_QUADRATIC_085)
        lambdaFctrInitialValue  = 0.95;
        lambdaFctrEffIterFactor = 0.85;
        lambdaFctrPwrFctr       = 2;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_POLYNOMIAL;
    case(RLS_CONFIG_QUADRATIC_075)
        lambdaFctrInitialValue  = 0.95;
        lambdaFctrEffIterFactor = 0.75;
        lambdaFctrPwrFctr       = 2;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_POLYNOMIAL;
    case(RLS_CONFIG_CUBIC_100)
        lambdaFctrInitialValue  = 0.95;
        lambdaFctrEffIterFactor = 1.00;
        lambdaFctrPwrFctr       = 3;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_POLYNOMIAL;
    case(RLS_CONFIG_CUBIC_085)
        lambdaFctrInitialValue  = 0.95;
        lambdaFctrEffIterFactor = 0.85;
        lambdaFctrPwrFctr       = 3;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_POLYNOMIAL;
    case(RLS_CONFIG_CUBIC_075)
        lambdaFctrInitialValue  = 0.95;
        lambdaFctrEffIterFactor = 0.75;
        lambdaFctrPwrFctr       = 3;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_POLYNOMIAL;
    case(RLS_CONFIG_HYPERBOLA)
        lambdaFctrInitialValue  = 0.95;
        lambdaFctrEffIterFactor = 0.45;
        lambdaFctrPwrFctr       = 1;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_HYPERBOLA;
    case(RLS_CONFIG_EXP)
        lambdaFctrInitialValue  = 0.95;
        lambdaFctrEffIterFactor = 0.45;
        lambdaFctrPwrFctr       = 1;
        lambdaFctrMode          = LAMBDA_FACTOR_MODE_EXPONENTIAL;
end

end