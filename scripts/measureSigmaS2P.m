%% Conversion from sigma_n to sigma_p

%% Create or load the Gabor spots stimuli
gaborSpotsFname = fullfile(cortical_obs_rootpath, 'stimulusgen', 'gaborSpots.mat');

if exist(gaborSpotsFname, 'file')
    load(gaborSpotsFname, 'gaborSpots');
else
    gaborSpots = createSpots(gaborSpotsFname);
end

gaborFlat = stackToFlat(gaborSpots.gaborStack);

%% Choose parameters for sResults
R = 1; S = .5;
X = 45; Y = 45;
G = 1;
C = 0.93; % for example

dVals = [1, 1.5, 2, 2.5, 3];
nVals = linspace(0.1,2,20);

%% Compute sResults

sResultsName = fullfile(cortical_obs_rootpath, 'stimulusgen', 'sResults.mat');
if exist(sResultsName, 'file')
    load(sResultsName, 'sResults')
else
    sResults = zeros(length(dVals), length(nVals));
    for dd = 1:length(dVals)
        disp(dd);
        for nn = 1:length(nVals)
            disp(nn);

            N = nVals(nn);  
            D = dVals(dd);
            params = [R, S, X, Y, D, G, N, C];

            predictions = socmodel_nogaborstep(params, gaborFlat);

            %
            predIm = reshape(predictions, length(dotPosPx), length(dotPosPx));
            %figure; imshow(predIm, []);

            %
            [xPts,yPts] = meshgrid(dotPosDeg);

            my2dGauss = fittype(@(b,s,X,Y)(b*exp(-(X.^2+Y.^2)/(2*s^2))), ...
                            'independent', {'X', 'Y'},...
                            'coefficients', {'b', 's'});
            opt = fitoptions(my2dGauss);
            opt.startpoint = [1, 1];
            fitobj = fit([xPts(:),yPts(:)], predIm(:), my2dGauss, opt);

            %figure, plot(fitobj), hold on,
            %plot3(X(:), Y(:), predIm(:), '.')

            sResults(dd, nn) = fitobj.s;
        end
    end
end

%% Fit D/sqrt(n) model for each D
aFit = zeros(size(dVals));
bFit = zeros(size(dVals));

for dd = 1:length(dVals)
    nsFit = fit(dVals(dd)*sqrt(nVals).^(-1)', sResults(dd, :)', 'poly1');
    aFit(dd) = nsFit.p1;
    bFit(dd) = nsFit.p2;
end

%% Plot slices
figure; hold all;

for dd = 1:length(dVals)
    plot(sResults(dd, :), 'o-')
end

for dd = 1:length(dVals)
    %plot(aFit(dd)*dVals(dd)*sqrt(nVals).^(-1) + bFit(dd), '-')
    plot(0.16*dVals(dd)*sqrt(nVals).^(-1)  - 0.05*dVals(dd) + 0.23, '-')
end
legend([arrayfun(@(x)(['sigma\_s = ', num2str(x)]), dVals, 'UniformOutput', false), arrayfun(@(x)(['sigma\_s = ', num2str(x)]), dVals, 'UniformOutput', false)]);

xlabel('n parameter')
ylabel('measured pRF size (sigma\_p)')
title('Relationship between n and sigma\_p for varying sigma\_s');


%% Plot the other slices
figure; hold all;
for nn = 1:length(nVals)
    plot(sResults(:, nn), 'o-')
end
xlabel('parameter for pRF size (sigma\_s)')
ylabel('measured pRF size (sigma\_p)')
title('Relationship between sigma\_p and sigma\_s for varying n');
legend(arrayfun(@(x)(['n = ', num2str(x)]), nVals, 'UniformOutput', false))
