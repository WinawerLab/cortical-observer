function params = stdObs_epaToParams(epa, imSzPx, pxPerDeg)
% STANDARD CORTICAL OBSERVER: (eccentricity, polar angle, area), imageSize, pixels PerDegree -> SOC parameters
%
% Convert the pRF model's pRF size parameter to the equivalent SOC model's
% pRF size parameter
%
% In: epa = voxs x 3 (ecc, polar angle, area)
%      imSz = 1 number (assumed square)
%      pxPerDeg = 1 number, pixels per degree
% Out: params = voxs x 5 (x, y, sigma_s, n, c)

    % Assert that all "area" numbers are 1, 2, or 3
    assert(all(any(cat(2, epa(:,3)==1, epa(:,3)==2, epa(:,3)==3), 2)), ...
        'Area must be 1, 2, or 3');
    
    nVox = size(epa,1);
    ecc = epa(:,1);
    area = epa(:,3);
    
    % N - data from Kay et al., PLoS CB 2013, Fig 9
    n = zeros(nVox,1);
    n(area == 1) = 0.18;
    n(area == 2) = 0.13;
    n(area == 3) = 0.12;

    % c - data from Kay et al., PloS CB 2013, Fig 9
    c = zeros(nVox,1);
    c(area == 1) = 0.93;
    c(area == 2) = 0.99;
    c(area == 3) = 0.99;
    
    % sigma_s - data from Harvey & Dumoulin, J.Neuro 2011, Fig 4
    sigma_p_deg = zeros(nVox, 1);
    sigma_p_deg(area == 1) = 0.5 + 0.10*ecc(area == 1);
    sigma_p_deg(area == 2) = 0.5 + 0.15*ecc(area == 2);
    sigma_p_deg(area == 3) = 0.5 + 0.27*ecc(area == 3);
    
    sigma_p_px = sigma_p_deg * pxPerDeg;
    sigma_s_px = stdObs_convertP2S(sigma_p_px, n);
    
    % x and y
    xy = stdObs_epToXy(epa(:, 1:2), imSzPx, pxPerDeg);
    
    % Put it all together
    params = cat(2, xy, sigma_s_px, n, c);
end