function sigmaS = stdObs_convertP2S(sigmaP, n)
% CONVERT SIGMA P-TO-S
% sigmaS = stdObs_convertP2S(sigmaP, n)
%
% Convert the pRF model's pRF size parameter to the equivalent SOC model's
% pRF size parameter
%
% Can handle a list of sigmaP, but currently only one n value at a time
    assert(length(n)==1)
    assert(any(size(sigmaP)==1))
    
    % Arithmetic:
    % sigmaP = 0.161*sigmaS*sqrt(n).^(-1)  + 0.249*sigmaS.^(-0.5) - 0.059;
    % 0 = 0.161*sigmaS*sqrt(n).^(-1)  + 0.249*sigmaS.^(-0.5) - 0.059 - sigmaP;
    
    sigmaS = zeros(1, length(sigmaP));
    for ii = 1:length(sigmaP)
        a = 0.161*sqrt(n).^(-1);
        b = 0.249;
        c = -0.059 - sigmaP(ii);

        % Arithmetic continued:
        % 0 = a*sigmaS + b*sigmaS.^(-0.5) + c
        % (-a*sigmaS - c)/b = sigmaS.^(-0.5)
        % (-a/b) * (sigmaS^1.5) + (-c/b) * (sigmaS^0.5) = 1
        % (-a/b) * (sigmaS^0.5)^3 + (-c/b) * (sigmaS^0.5) - 1 = 0

        sqrtSigmaS = roots([-a/b, 0, -c/b, -1]);
        sqrtSigmaS = sqrtSigmaS(sqrtSigmaS > 1); % only one out of three is right
        sigmaS(ii) = sqrtSigmaS.^2;
    end
end

% To check for round tripping:
% sigmaSVals = [2, 4, 6, 8, 10, 12, 14];
% p = stdObs_convertS2P(sigmaSVals, 0.3)
% roundTrip = stdObs_convertP2S(p, 0.3)
% assert(all(abs(roundTrip - sigmaSVals) < 0.00001))