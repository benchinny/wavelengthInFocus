%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
subjNumAll = [11 13 15 20 26 27 28 30];
saveFolder = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\data\packagedAccommodationData\';

for i = 1:length(subjNumAll)
    [zCoeffCell,colorConditionsRGB,optDistance] = ARCnlz_packageData(subjNumAll(i),dataPath);
    fileName = ['packagedDataS' num2str(i)];
    save([saveFolder fileName],'colorConditionsRGB','optDistance','zCoeffCell');
end