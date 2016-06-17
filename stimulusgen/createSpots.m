function gaborSpots = createSpots(outputdir)
% Optional outputdir parameter

    %% Set up image size parameters
    res = 400;                     % native resolution that we construct at
    totalfov = 12;                 % total number of degrees for image
    cpd = 3;
    pxPerDeg = res/totalfov;
    degPerPx = totalfov/res;

    maskFrac = 11.5/12; 

    radius = totalfov/2;           % radius of image in degrees
    cpim = totalfov*cpd;           % cycles per image that we are aiming for
    spacing = res./cpim;            % pixels to move from one cycle to the next

    span = [-0.5, 0.5];             % dynamic range, to put into 'imshow' etc.

    %% Choose which bandpass filter to use
    bandwidth = 1;
    fltsz = 31;
    flt = mkBandpassCosine(res, cpim, bandwidth, fltsz, 0);

    %% Parameters for spots spaced apart
    dotSizeDeg = 0.25; 
    dotSizePx = dotSizeDeg * pxPerDeg;

    mid = ceil(([res,res]+1)/2);

    dotPosDeg = -totalfov/2+dotSizeDeg/2:dotSizeDeg:totalfov/2-dotSizeDeg/2;
    dotPosPx = dotPosDeg * pxPerDeg + res/2;

    %% Create spots spaced apart
    imStack = zeros(res, res, length(dotPosPx).^2);
    whichIm = 1;
    for x = 1:length(dotPosPx)
        disp(x)
        for y = 1:length(dotPosPx)
            imStack(:,:,whichIm) = mkDisc([res,res], dotSizePx/2, [dotPosPx(x), dotPosPx(y)], 0, [0.5, 0]);
            whichIm = whichIm + 1;
        end
    end

    %%
    outputSz = 150; padSz = 30;
    imStack = imresize(imStack, [outputSz, outputSz]);
    imStack = padarray(imStack, [padSz/2, padSz/2, 0], 0, 'both');

    imFlat = stackToFlat(imStack);
    %imFilt = imfilter(im, flt, 'circular');

    resizedPxPerDeg = outputSz/totalfov;
    resizedDegPerPx = 1/resizedPxPerDeg;

    %% Put through gabor pipeline
    numor = 8; numph = 2;
    gaborFlat = gaborenergy(imFlat, numor, numph, cpim);
    gaborStack = flatToStack(gaborFlat, 1);

    % Save out
    if ~exist(outputdir, 'var')
        outputdir = fullfile(cortical_obs_rootpath, 'stimulusgen');
    end
        % not currently bothering with a dedicated data directory,
        % or with subdirectories called datestr(now,'yyyy-mm-dd'),
        % but you can do that if you want to avoid clobbering gaborSpots.mat
        
    if ~exist(outputdir, 'dir')
        mkdir(outputdir);
    end
    gaborSpots = {};
    gaborSpots.bandwidth = bandwidth;
    gaborSpots.numor = numor;
    gaborSpots.numph = numph;
    gaborSpots.function = 'gaborenergy';
    gaborSpots.inputImStack = imStack;
    gaborSpots.gaborStack = gaborStack;
    gaborSpots.generatingFile = 'createSpots.m';
    gaborSpots.dateSaved = datestr(now);

    outputFile = 'gaborSpots.mat';
    save(fullfile(outputdir, outputFile), 'gaborSpots');
end


