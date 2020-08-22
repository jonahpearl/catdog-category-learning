% pre-process data for category training, once it's been converted from REX
% to MATLAB.

%% Define data paths
clearvars
close all

EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';

%% (Run once) Combine monkeys' data

% path to the category training ("skip" XMA) data
skip_path = 'SXMA/Monkey_structs/';

% load Marta's data first
Monkeys = load(fullfile(EXT_HD, skip_path, 'Marta_skip_fix_cat_xma_raw_behav.mat'), 'Monkeys');
Monkeys = Monkeys.Monkeys;

% append Max's data
Data = load(fullfile(EXT_HD, skip_path, 'Max_skip_fix_cat_xma_raw_behav.mat'), 'Monkeys');
Data = Data.Monkeys;
Monkeys(2).Name = Data.Name;
Monkeys(2).Sessions = Data.Sessions;
clear Data

% Save combined data
save(fullfile(EXT_HD, skip_path, 'MaxMarta_skip_fix_cat_xma_raw_behav.mat'), 'Monkeys')

%% Load combined data
% path to the category training ("skip" XMA) data
skip_path = 'SXMA/Monkey_structs/';

% load combined data
load(fullfile(EXT_HD, skip_path, 'MaxMarta_skip_fix_cat_xma_raw_behav.mat'), 'Monkeys')

%% Add trial-unique data to end of each Sessions struct

% path to trial-unique day of data
uni_path = 'USXMA/Monkey_structs';

% load data
uni = load(fullfile(EXT_HD, uni_path, 'MaxMarta_skip_fix_cat_UNI_xma_raw_behav.mat'), 'Monkeys');
uni = uni.Monkeys;

for m = 1:length(Monkeys)
    
    % check for ordering issues
    mk_name = regexp(Monkeys(m).Name, '(\w*)?_skip\w*', 'tokens', 'once');
    mk_name = mk_name{1};
    if isempty(regexp(uni(m).Name, mk_name, 'once'))
        error('Monkeys in wrong order in Monkeys vs uni structs')
    end
    
    % add UNI data at end of category training data
    end_plus1 = length(Monkeys(m).Sessions)+1;
    Monkeys(m).Sessions(end_plus1) = uni(m).Sessions;
end

%% Get date strs and short names of sessions
marta_xls = fullfile(EXT_HD, 'RecordingMarta.xlsx');
max_xls = fullfile(EXT_HD, 'RecordingMax.xlsx');

for m = 1:length(Monkeys)
    % get short session names
    date_strs = {Monkeys(m).Sessions.DateStr};
    if regexp(Monkeys(m).Name, 'Marta\w*')
        short_names = get_short_names(date_strs, marta_xls, '\w*cats\w*');
    else
        short_names = get_short_names(date_strs, max_xls, '\w*cats\w*');
    end
    for i = 1:length(Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).ShortName = short_names{i};
    end
end


%% Add "Categories" field to sessions (20/20 vs 240/240)
for m = 1:length(Monkeys)
    for i = 1:length(Monkeys(m).Sessions)
        if regexp(Monkeys(m).Sessions(i).ShortName, 'UNI')
            Monkeys(m).Sessions(i).Categories = {1:240, 241:480};
        else
            Monkeys(m).Sessions(i).Categories = {1:20, 21:40};
        end
    end
end


%% Remove extra UNI trials
for m = 1:length(Monkeys)
    
    % error check
    completed = [Monkeys(m).Sessions(end).TrialInfo.Completed];
    completed_imgs = [Monkeys(m).Sessions(end).TrialInfo(completed == 1).ImageID];
    if ~all(1:480 == sort(completed_imgs(1:480)))
        error('%s UNI has repeat of an image before all are shown once', Monkeys(m).Name)
    end
    
    % find last trial of first img presentations
    for j = 480:length(Monkeys(m).Sessions(end).TrialInfo)
        if sum([Monkeys(m).Sessions(end).TrialInfo(1:j).Completed]) == 480
            final_UNI_ind = j;
            break
        else
            continue
        end
    end
    
    % delete trials after that
    Monkeys(m).Sessions(end).TrialInfo(final_UNI_ind+1:end) = [];
    Monkeys(m).Sessions(end).TimesBT(final_UNI_ind+1:end) = [];
    Monkeys(m).Sessions(end).CodesBT(final_UNI_ind+1:end) = [];
    
    % only look at first presentations in CueInfo
    for k = 1:length(Monkeys(m).Sessions(end).CueInfo)
        times_on = Monkeys(m).Sessions(end).CueInfo(k).Times_on;
        trial_nums = Monkeys(m).Sessions(end).CueInfo(k).TrialNum;
        Monkeys(m).Sessions(end).CueInfo(k).Times_on = times_on(1);
        Monkeys(m).Sessions(end).CueInfo(k).TrialNum = trial_nums(1);
        Monkeys(m).Sessions(end).CueInfo(k).NumApp = 1;
    end
    
    % error check again
    if max([Monkeys(m).Sessions(end).CueInfo.TrialNum]) ~= final_UNI_ind
        error('Final trial and max trial in CueInfo do not align')
    end
end

%% Get reaction time info

for m = 1:length(Monkeys) 
    for i = 1:length(Monkeys(m).Sessions) 
        for j = 1:length(Monkeys(m).Sessions(i).CodesBT)
            
            % relevant trial info
            codes = Monkeys(m).Sessions(i).CodesBT{j};
            times = Monkeys(m).Sessions(i).TimesBT{j} * 1000; % to match with ms of rasters
            
            if Monkeys(m).Sessions(i).TrialInfo(j).ROG == 1
                Monkeys(m).Sessions(i).TrialInfo(j).RT = ...
                    times(codes == 1141) - times(codes == 1121);
                
            elseif Monkeys(m).Sessions(i).TrialInfo(j).ROR == 1
                Monkeys(m).Sessions(i).TrialInfo(j).RT = ...
                    times(codes == 1036) - times(codes == 1100);
                
            else
               Monkeys(m).Sessions(i).TrialInfo(j).RT = NaN; 
            end
        end
    end
end

%% Save data
save(fullfile(EXT_HD, skip_path, 'MaxMarta_skip_preProcessed.mat'), 'Monkeys')

