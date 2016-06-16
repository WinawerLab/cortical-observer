%%
%%%%%%%%%% PART 1 - parameters %%%%%%%%%%
%
%% Load fMRI data as deposited by Noah's script
retData = load('/Volumes/server/Projects/SOC/data/retinotopy/wl_subj001.mat');

epaRaw = cat(2, retData.eccentricity', retData.polar_angle', retData.visual_area');
nVoxRaw = size(epaRaw,1);
    % Out: voxs x 3 (ecc, polar angle, area)

% Remove visual areas 1.5 and 2.5, which are voxels on the borders
discardAreas = or(epaRaw(:,3) == 1.5, epaRaw(:,3) == 2.5);
discardEcc = epaRaw(:,1) > 20;
discardVoxels = or(discardAreas, discardEcc);
epa = epaRaw(~discardVoxels, :);

%% Compute parameters

% Obtain SOC parameters for voxels
params = stdObs_epaToParams(epa);
    % In: voxs x 3 (ecc, polar angle, area)
    % Out: voxs x 3 (sigma_s, n, c)
  
% (Sanity check: do pRF sizes get larger with area?)
% mean(params(epa(:,3) == 1, 1)) - 5.57
% mean(params(epa(:,3) == 2, 1)) - 5.31
% mean(params(epa(:,3) == 3, 1)) - 7.71
    
%% Sanity check hemifield (demo)
% Convert eccentricity and polar angle (in degrees) to x,y (in pixels)
imSzPxDemo = 300; % images assumed to be square
pxPerDegDemo = 10; % TODO get real values
xyDemo = stdObs_epToXy(epa(:,1:2), imSzPxDemo, pxPerDegDemo);
    % In: voxs x 2 (ecc, polar angle), 1x2 (size), 1
    % Out: voxs x 2 (x, y)
    
figure(1); hold on;
scatter(xyDemo(epa(:,2) >=0, 1), xyDemo(epa(:,2) >=0, 2), 'ro');
scatter(xyDemo(epa(:,2) < 0, 1), xyDemo(epa(:,2) < 0, 2), 'bo');
axis ij; xlim([0, imSzPxDemo]); ylim([0, imSzPxDemo]);

%%
%%%%%%%%%% PART 2 - predictions %%%%%%%%%%
%

%% Prepare data and parameters
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

xy = stdObs_epToXy(epa(:,1:2), imSzPx, pxPerDeg);
xyParams = [xy, params];

%% Make predictions
if exist('gaborOutput', 'var')
    preds = stdObs_predict(xyParams, imFlat, cpIm, gaborOutput);
else
    [preds, gaborOutput] = stdObs_predict(xyParams, imFlat, cpIm);
end
    % In: voxs x 5 (x, y, sigma_s, n, c), px x nIms, 1
    % Out: voxs x nIms 
    
%% Plot *all* predictions versus data
rois = {'RV1', 'RV2', 'RV3'};
for area = 1:3 % just V1 for now
    roiIdx = strInCellArray(rois{area}, data.roiNames);
    roiData = data.roiBetamn{roiIdx};
    
    predStdObs = mean(preds(epa(:,3) == area, :), 1);
    predStdObs_scaled = predStdObs * mean(roiData(:)) / mean(predStdObs);

    setupBetaFig;
    plotWithColors(roiData, data.plotOrder, data.plotNames, data.catColors)
    plot(predStdObs_scaled(data.plotOrder), 'go');
end
 
