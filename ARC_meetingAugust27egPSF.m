%%

foldername = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/ARC/';

filenames = {'S1013V10_AFC_RightACL0_2408221213.mat' ...
              'S1013V10_AFC_RightACL0_2408221218.mat' ...
              'S1013V10_AFC_RightACL0_2408221218.mat'};

for i = 1:length(filenames)
    load([foldername filenames{i}]);
    AFCp.rspAcu = AFCp.rspAcu';
    if i==1
        AFCpAll = AFCp;
    else
        AFCpAll = structmerge(AFCpAll,AFCp);
    end
end