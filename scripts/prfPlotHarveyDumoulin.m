% Recreate Fig 4 in Harvey & Dumoulin J.Neuro 2011
% http://www.ncbi.nlm.nih.gov/pubmed/21940451

% This is used as our estimate of sigma_p, in degrees,
% by visual area

ecc = linspace(0.5, 5.5, 21);
prf = struct();
prf.('V1') = 0.5 + 0.10*ecc;
prf.('V2') = 0.5 + 0.15*ecc;
prf.('V3') = 0.5 + 0.27*ecc;
prf.('hV4') = 0.5 + 0.4*ecc;

figure(1); clf; hold on;
plot(ecc, prf.('V1'), 'ko-')
plot(ecc, prf.('V2'), 'ro-')
plot(ecc, prf.('V3'), 'go-')
plot(ecc, prf.('hV4'), 'bo-')
xlim([0,6]); ylim([0,5]);
xlabel('Eccentricity'); ylabel('sigma p (degrees of visual angle)')
title('Modeled sigma p relationship')
axis('square');

n.('V1') = 0.18;
n.('V2') = 0.13;
n.('V3') = 0.12;
n.('hV4') = 0.12;

figure(2); clf; hold on;
plot(ecc, stdObs_convertS2P(prf.('V1'), n.('V1')), 'ko-')
plot(ecc, stdObs_convertS2P(prf.('V2'), n.('V2')), 'ro-')
plot(ecc, stdObs_convertS2P(prf.('V3'), n.('V3')), 'go-')
plot(ecc, stdObs_convertS2P(prf.('hV4'), n.('hV4')), 'bo-')
xlim([0,6]); ylim([0,5]);
xlabel('Eccentricity'); ylabel('sigma s')
title('Modeled sigma s relationship')
axis('square');