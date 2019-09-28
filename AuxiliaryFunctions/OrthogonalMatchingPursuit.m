function [ vX_hat ] = OrthogonalMatchingPursuit( vB, mA, cardinalityK, noiseThr )

vX_hat         = zeros(size(mA, 2), 1);
vSupport       = zeros(size(mA, 2), 1);
vResidual      = vB;
support_length = 0;

if(nargin <= 3)
    noiseThr = 0;
end

while((support_length < cardinalityK) && (norm(vResidual) > noiseThr))
    %-- Sweep:
    vProjections = mA' * vResidual;
    [~, j0]      = max( abs(vProjections) );

    %-- Update Support:
    support_length           = support_length + 1;
    vSupport(support_length) = j0;

    %-- Update Provisional Solution:
    vA_j0 = mA(:, j0);
    if(support_length == 1)
        mM_inv = 1; % 1 / (vA_j0' * vA_j0);
        mA_s   = mA(:, j0);
        vX_s   = vProjections(j0);
    else
        vD     = mA_s' * vA_j0;
        vMd    = mM_inv * vD;
        p      = 1 / (1 - vD' * vMd);
        
        mM_inv = [mM_inv + (p * vMd) * vMd', -p * vMd;
                  -p * vMd',                 p        ];
        mA_s   = mA(:, vSupport(1:support_length));
        vX_s   = mM_inv * (mA_s'  * vB);
                         
%         vX_s = mA(:, vSupport(1:support_length)) \ vB
    end
    %-- Update Residual:
    vResidual = vB - mA_s * vX_s;
end

vSupport(vSupport == 0) = [];
vX_hat(vSupport) = vX_s;


end

