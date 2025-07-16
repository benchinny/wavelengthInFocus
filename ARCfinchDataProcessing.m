%%

Finch_primaries_red = readtable('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/Finch_et_al_primaries_Red.csv');

%%

wave = round(Finch_primaries_red.Wavelength,1);

power = round(Finch_primaries_red.Intensity,3);

waveUnq = unique(wave);
powerUnq = [];

for i = 1:length(waveUnq)
    ind = abs(wave-waveUnq(i))<0.0001;
    powerUnq(i) = mean(power(ind));
end

waveInterp = 616:1:700;
powerInterp = interp1(waveUnq,powerUnq,waveInterp,'pchip');
powerInterp(powerInterp<0) = 0;

%%

waveBluePunq = unique(waveBlueP);
waveRedPunq = unique(waveRedP(1:end-1));
powerBluePunq = [];
powerRedPunq = [];

for i = 1:length(waveBluePunq)
    ind = waveBlueP==waveBluePunq(i);
    powerBluePunq(i) = mean(powerBlueP(ind));
end

for i = 1:length(waveRedPunq)
    ind = waveRedP==waveRedPunq(i);
    powerRedPunq(i) = mean(powerRedP(ind));
end

waveBluePunq = waveBluePunq';
waveRedPunq = waveRedPunq';

waveBluePunq = [380 waveBluePunq 780];
waveRedPunq = [380 waveRedPunq 780];
powerBluePunq = [0 powerBluePunq 0];
powerRedPunq = [0 powerRedPunq 0];
wave = 380:4:780;
d = struct;
d.wave = wave;
d.spd(:,1) = interp1(waveRedPunq,powerRedPunq,wave,'linear');
d.spd(:,3) = interp1(waveBluePunq,powerBluePunq,wave,'linear');
d.spd(:,2) = zeros(size(wave));

%%

colorCell = {'RedVsGreen' 'RedVsBlue' 'RedVsViolet' 'OrangeVsBlue' 'OrangeVsViolet' 'GreenVsViolet'};
dataFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/';
lumRatiosTrue = 0:1.25:10;
accommodationFinch = [];

for i = 1:length(colorCell)
    FinchData = readtable([dataFolder 'Finch_data_' colorCell{i} '.csv']);
    for j = 1:length(lumRatiosTrue)
        ind = abs(FinchData.LuminanceRatio-lumRatiosTrue(j))<0.1;
        accommodationFinch(j,i) = mean(FinchData.Accommodation(ind));
    end
end

FinchData = struct;
FinchData.LumRatio = fliplr(lumRatiosTrue./10)';
FinchData.Accommodation = accommodationFinch;

save([dataFolder 'Finch_data_all.mat'],'FinchData');

%%  

clear all;
