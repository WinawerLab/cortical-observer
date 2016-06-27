function sigmaP = stdObs_convertS2P(sigmaS, n)
% CONVERT SIGMA S-TO-P
% sigmaS = stdObs_convertS2P(sigmaP, n)
%
% Convert the SOC model's pRF size parameter to the equivalent pRF model
% pRF size parameter 

    sigmaP = 0.161*sigmaS*sqrt(n).^(-1)  + 0.249*sigmaS.^(-0.5) - 0.059;
end
