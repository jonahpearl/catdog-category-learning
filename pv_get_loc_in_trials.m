% get loc in trial and trial num data (based on script
% pv_get_interval_spike_counts)

%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';

% Load the data, broken up into sessions for storage
Data = load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_ni.mat'));

% Stitch the sessions back into one big structure
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
glmRecordPath = 'XMA2/Monkey_structs/GLM_Records.mat';
pv_path = 'XMA2/Monkey_structs';

% Load behavioral data
load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_behav_and_metaNI.mat')) % behavior and neural summaries, but w/o spike times



%% Get date strs and short names of sessions
marta_xls = fullfile(EXT_HD, 'RecordingMarta.xlsx');
max_xls = fullfile(EXT_HD, 'RecordingMax.xlsx');

for m = 1:length(Monkeys)
    date_strs = {Monkeys(m).Sessions.DateStr};
    if regexp(Monkeys(m).Name, 'Marta\w*')
        shortNames = get_short_names(date_strs, marta_xls, '\w*cats\w*');
    else
        shortNames = get_short_names(date_strs, max_xls, '\w*cats\w*');
    end
    for i = 1:length(Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).ShortName = shortNames{i};
    end
end

%% Collect loc in trial data

% Parameter for organizing the code. It's important that we decide on how
% the spike counts will be organized, so that we can align the correct
% spike counts with the correct trial information later on!
catgs = {1:260, 261:520};
spikeCountPath = 'XMA2/Spike_count_mats';

for m = 1:length(Monkeys)
    rSessions = 1:length(Monkeys(m).Sessions); % "relevant sessions"
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        
        
        % === Positions images were shown, sorted by category ===
        fname_base = sprintf('%%s_imgLocInTrials.mat');
        catg1 = Monkeys(m).Sessions(sessn).CueInfo(catgs{1});
        catg2 = Monkeys(m).Sessions(sessn).CueInfo(catgs{2});
        catg1_posns = vertcat(catg1.Loc_in_trial);
        catg2_posns = vertcat(catg2.Loc_in_trial); 
        X = zeros(length(catg1_posns)+length(catg2_posns),1);
        X(:) = [catg1_posns; catg2_posns];
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf(fname_base, MonkID);
        save(fullfile(EXT_HD, spikeCountPath, fileName), 'X');
        
        
        % === Repeat above but for trial number ===
        fname_base = sprintf('%%s_trialNum.mat');
        catg1 = Monkeys(m).Sessions(sessn).CueInfo(catgs{1});
        catg2 = Monkeys(m).Sessions(sessn).CueInfo(catgs{2});
        catg1_tn = vertcat(catg1.TrialNum);
        catg2_tn = vertcat(catg2.TrialNum);
        X = zeros(length(catg1_tn)+length(catg2_tn),1);
        X(:) = [catg1_tn; catg2_tn];
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf(fname_base, MonkID);
        save(fullfile(EXT_HD, spikeCountPath, fileName), 'X');
    end
end






