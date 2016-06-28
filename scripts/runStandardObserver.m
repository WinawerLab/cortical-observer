%%
%%%%%%%%%% PART 1 - parameters %%%%%%%%%%
%
% epaRaw: matrix containing eccentricity (deg visual angle), polar angle
%           (deg from upper vertical), and visual area (1-3)
%    epa: same as epaRaw but trimmed for border voxels, etc
%
%% Load retinotopy template as deposited by Noah's script
templateOutput = load('/Volumes/server/Projects/SOC/data/retinotopy/wl_subj001.mat');

epaRaw = cat(2, templateOutput.eccentricity', templateOutput.polar_angle', templateOutput.visual_area');
nVoxRaw = size(epaRaw,1);
    % Out: voxs x 3 (ecc, polar angle, area)

% Remove visual areas 1.5 and 2.5, which are voxels on the borders
discardAreas = or(epaRaw(:,3) == 1.5, epaRaw(:,3) == 2.5);

% Remove eccentricities that are large (we might skip this or do it later
% for visualization)
discardEcc = epaRaw(:,1) > 20;
discardVoxels = or(discardAreas, discardEcc);
epa = epaRaw(~discardVoxels, :);
    
% %% Sanity check hemifield (example/dummy)
% % Convert eccentricity and polar angle (in degrees) to x,y (in pixels)
% imSzPxDemo = 400; % images assumed to be square
% pxPerDegDemo = 10; % these are just dummy values
% xyDemo = stdObs_epToXy(epa(:,1:2), imSzPxDemo, pxPerDegDemo);
%     % In: voxs x 2 (ecc, polar angle), 1x2 (size), 1
%     % Out: voxs x 2 (x, y)
%     
% figure(1); clf; hold on;
% scatter(xyDemo(epa(:,2) >=0, 1), xyDemo(epa(:,2) >=0, 2), 'ro');
% scatter(xyDemo(epa(:,2) < 0, 1), xyDemo(epa(:,2) < 0, 2), 'bo');
% 
% axis ij; xlim([0, imSzPxDemo]); ylim([0, imSzPxDemo]);

%%
%%%%%%%%%% PART 2 - predictions %%%%%%%%%%
%

%% Prepare data
data = load_subj001_2015_10_22();
imStack = data.stimuli.imStack(:,:,:,1);
imStack = (double(imStack)- 127)/255;
imStack = imresize(imStack, [150, 150]);
imFlat = stackToFlat(imStack);
    % TODO: just one frame/exemplar per class for now

cpd = data.stimuli.cpd;
totalFovDeg = data.stimuli.totalfov;
cpIm = cpd * totalFovDeg;

imSzPx = sqrt(size(imFlat,1));
pxPerDeg = imSzPx / totalFovDeg;

%% Prepare parameters
params = stdObs_epaToParams(epa, imSzPx, pxPerDeg);
    % In: voxs x 3 (ecc, polar angle, area)
    % Out: voxs x 5 (x, y, sigma_s, n, c)
  
% (Sanity check: do pRF sizes get larger with area?)
% mean(params(epa(:,3) == 1, 3)) - 24.38 (pixels)
% mean(params(epa(:,3) == 2, 3)) - 26.26 
% mean(params(epa(:,3) == 3, 3)) - 31.29

%% Make predictions
if exist('gaborOutput', 'var')
    preds = stdObs_predict(params, imFlat, cpIm, gaborOutput);
else
    [preds, gaborOutput] = stdObs_predict(params, imFlat, cpIm);
end
    % In: voxs x 5 (x, y, sigma_s, n, c), px x nIms, 1
    % Out: voxs x nIms 
    
%% Plot predictions versus data
rois = {'RV1', 'RV2', 'RV3'};
for area = 1:3
    roiIdx = strInCellArray(rois{area}, data.roiNames);
    roiData = data.roiBetamn{roiIdx};
    
    predStdObs = mean(preds(epa(:,3) == area, :), 1);
    predStdObs_scaled = predStdObs * mean(roiData(:)) / mean(predStdObs);

    setupBetaFig;
    plotWithColors(roiData, data.plotOrder, data.plotNames, data.catColors)
    plot(predStdObs_scaled(data.plotOrder), 'go');
end
 
