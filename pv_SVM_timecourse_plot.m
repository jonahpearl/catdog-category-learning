% plot SVM timecourse results

%% Paramd and load SVM record
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
svmRecordPath = 'XMA2/Monkey_structs/SVM_Records.mat';
pv_path = 'XMA2/Monkey_structs';

% Load record
fullRecordPath = fullfile(EXT_HD, svmRecordPath);
load(fullRecordPath, 'Record');

%% Find desired SVM data and load
ID = 439280;
fNames = fields(Record);
row = find([Record.ID] == ID);
for f = 1:length(fNames)
    if strcmp(fNames{f}, 'MatchInputParams')
        fprintf('%s: %s \n', fNames{f}, mat2str(Record(row).(fNames{f})))
    else
        fprintf('%s: %s \n', fNames{f}, string(Record(row).(fNames{f})))
    end
    
end
