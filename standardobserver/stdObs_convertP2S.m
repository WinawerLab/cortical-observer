function sigmaS = stdObs_convertP2S(sigmaP, n)
% CONVERT SIGMA P-TO-S
% sigmaS = stdObs_convertP2S(sigmaP, n)
%
% Convert the pRF model's pRF size parameter to the equivalent SOC model's
% pRF size parameter

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
        sigmaS(ii) = sqrtSigmaS(1).^2; % only take the first, positive one
    end
end

