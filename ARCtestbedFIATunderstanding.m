%%

pupilRadiusMm = 20;
trueFocusDistD = 3.5;
distPupil2camCorrectMm = 100;
distPupil2camActualMm = 180;

radDotCircleMm = (1000/trueFocusDistD-distPupil2camActualMm)*pupilRadiusMm/(1000/trueFocusDistD);

estimatedFocusDistMm = pupilRadiusMm*distPupil2camCorrectMm/(pupilRadiusMm-radDotCircleMm);

estimatedFocusDistD = 1000/estimatedFocusDistMm
