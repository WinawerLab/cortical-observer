%% Conversion from sigma_s_pix to sigma_p_pix
% This script measures the relationship between "sigma_s" (the
% parameter in the SOC equation, also called D some places)
% and "sigma_p", the empirical measurement of the width of the pRF,
% as measured by measuring the response to point stimuli at different
% locations in the image.

% The procedure is roughly as follows:
%   1) Measure responses to point stimuli
%   2) Fit the responses with a Gaussian and measure its std. dev.
%   3) Repeat for varying values of the nonlinearity n
%   4) Describe the relationship between sigma_s and _p as a function of n
%
% The punchline:
% sigma_p = 0.161*sigma_s*sqrt(n).^(-1)  + 0.249*sigma_s.^(-0.5) - 0.059
%
% This relationship is enshrined in stdObs_convertS2P and stdObs_convertP2S

%% Create or load the Gabor spots stimuli (point stimuli)
gaborSpotsFname = fullfile(cortical_obs_rootpath, 'stimulusgen', 'gaborSpots.mat');

if exist(gaborSpotsFname, 'file')
    load(gaborSpotsFname, 'gaborSpots');
else
    gaborSpots = createSpots(gaborSpotsFname);
end

gaborFlat = stackToFlat(gaborSpots.gaborStack);

%% Choose parameters to use when predicting responses to point stimuli
R = 1; S = .5; % haven't yet tested with different R,S values
X = 45; Y = 45;
G = 1;
C = 0.93; % haven't yet tested with different C values

sigmaSVals = [2, 4, 6, 8, 10, 12, 14]; % in pixels
nVals = linspace(0.1,2,20);

%% Compute or load sigma P fit to the sigma S

sigmaPFitName = fullfile(cortical_obs_rootpath, 'stimulusgen', 'sigmaFitResults.mat');
if false %exist(sigmaPFitName, 'file')
    load(sigmaPFitName, 'sigmaPFit')
else
    sigmaPFit = struct();
    sigmaPFit.sigmaSVals = sigmaSVals;
    sigmaPFit.nVals = nVals;
    sigmaPFit.sigmaPFit = zeros(length(sigmaSVals), length(nVals));
    
    for ss = 1:length(sigmaSVals)
        disp(ss);
        for nn = 1:length(nVals)
            disp(nn);

            N = nVals(nn);  
            D = sigmaSVals(ss);
            params = [R, S, X, Y, D, G, N, C];

            predictions = socmodel_nogaborstep(params, gaborFlat);

            %
            predIm = reshape(predictions, length(gaborSpots.dotPosPx), length(gaborSpots.dotPosPx));
            %figure; imshow(predIm, []);

            %
            [xPts,yPts] = meshgrid(gaborSpots.dotPosDeg);

            my2dGauss = fittype(@(b,s,X,Y)(b*exp(-(X.^2+Y.^2)/(2*s^2))), ...
                            'independent', {'X', 'Y'},...
                            'coefficients', {'b', 's'});
            opt = fitoptions(my2dGauss);
            opt.startpoint = [1, 1];
            fitobj = fit([xPts(:),yPts(:)], predIm(:), my2dGauss, opt);

            %figure, plot(fitobj), hold on,
            %plot3(X(:), Y(:), predIm(:), '.')

            sigmaPFit.sigmaPFit(ss, nn) = abs(fitobj.s);
        end
    end
    save(sigmaPFitName, 'sigmaPFit')
end

%% Fit sigma_s/sqrt(n) model for each sigma_s separately at a wide spread of n's
aFit = zeros(size(sigmaSVals));
bFit = zeros(size(sigmaSVals));

for ss = 1:length(sigmaSVals)
    nsFit = fit(sigmaSVals(ss)*sqrt(nVals).^(-1)', sigmaPFit.sigmaPFit(ss, :)', 'poly1');
    aFit(ss) = nsFit.p1;
    bFit(ss) = nsFit.p2;
    
    %figure; hold on;
    %plot(sigmaSVals(sigma_s)*sqrt(nVals).^(-1)',sigmaPFit.sigmaPFit(sigma_s, :)', 'o-')
    %plot(nsFit);
end

%% The slope is constant, but the intercept varies smoothly with n, so zoom in on it
figure; hold on;
plot(sigmaSVals, bFit, 'bo-');
xlabel('sigmaSvals'), ylabel('bFit'); 

%pFit = fit(sigmaSVals'.^(-0.5), bFit', 'poly1');
plot(sigmaSVals, 0.249*sigmaSVals.^(-0.5)-0.059, 'r');

%% Plot slices
figure; hold all;

for ss = 1:length(sigmaSVals)
    plot(sigmaPFit.sigmaPFit(ss, :), 'o-')
end

for ss = 1:length(sigmaSVals)
    %plot(aFit(dd)*dVals(dd)*sqrt(nVals).^(-1) + bFit(dd), '-')
    %plot(0.16*dVals(dd)*sqrt(nVals).^(-1)  - 0.05*dVals(dd) + 0.23, '-')
       % ^^^^ old guess, but at larger values it's clearly not linear, try
       % again...
       
    % THE PUNCHLINE:
    % plot(0.161*sigmaSVals(ss)*sqrt(nVals).^(-1)  + 0.249*sigmaSVals(ss).^(-0.5) - 0.059, '-')
    plot(arrayfun(@(n)(stdObs_convertS2P(sigmaSVals(ss), n)), nVals))
    
end
legend([arrayfun(@(x)(['sigma\_s = ', num2str(x)]), sigmaSVals, 'UniformOutput', false), arrayfun(@(x)(['sigma\_s = ', num2str(x)]), sigmaSVals, 'UniformOutput', false)]);

xlabel('n parameter')
ylabel('measured pRF size (sigma\_p)')
title('Relationship between n and sigma\_p for varying sigma\_s');


%% Plot the other slices
figure; hold all;
for nn = 1:length(nVals)
    plot(sigmaSVals, sigmaPFit.sigmaPFit(:, nn), 'o-')
    plot(sigmaSVals, stdObs_convertS2P(sigmaSVals, nVals(nn)), 'o-')
end
xlabel('parameter for pRF size (sigma\_s)');
ylabel('measured pRF size (sigma\_p)');
title('Relationship between sigma\_p and sigma\_s for varying n');
legend(arrayfun(@(x)(['n = ', num2str(x)]), nVals, 'UniformOutput', false))
