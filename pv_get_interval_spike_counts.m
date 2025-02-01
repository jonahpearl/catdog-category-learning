% get spike counts for all intervals of interest and store them

%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';

% Use this dir for the 260/260 PF and main cat/dog training
% pv_path = 'XMA2/Monkey_structs';
% spikeCountPath = 'XMA2/Spike_count_mats';  % main data
% Data = load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_ni.mat'));
% catgs = {1:260, 261:520};

% Use this dir for the initial 20/20 PF (control analysis)
% pv_path = 'XMA/Monkey_structs';
% spikeCountPath = 'XMA/Spike_count_mats';  % initial 20/20 PF
% Data = load(fullfile(EXT_HD, pv_path, 'Marta_fix_cat_xma_ni.mat'));
% Data = load(fullfile(EXT_HD, pv_path, 'Max_fix_cat_xma_ni.mat'));
% catgs = {1:20, 21:40};

% Use this dir for car/truck data
pv_path = 'XMA2/Monkey_structs';
spikeCountPath = 'XMA2/Spike_count_mats_CARTRUCK';  % car truck skip
Data = load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_CARTRUCK_ni.mat'));
catgs = {1:260, 261:520};

% Stitch the sessions back into one big structure
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data
% Monkeys = Data.Monkeys;

%% Get date strs and short names of sessions
marta_xls = fullfile(EXT_HD, 'RecordingMarta.xlsx');
max_xls = fullfile(EXT_HD, 'RecordingMax.xlsx');

for m = 1:length(Monkeys)
    date_strs = {Monkeys(m).Sessions.DateStr};
    if regexp(Monkeys(m).Name, 'Marta\w*')
%         shortNames = get_short_names(date_strs, marta_xls, '\w*cats\w*');
        shortNames = get_short_names(date_strs, marta_xls, '\w*cars\w*');
    else
%         shortNames = get_short_names(date_strs, max_xls, '\w*cats\w*');
        shortNames = get_short_names(date_strs, max_xls, '\w*cars\w*');
    end
    for i = 1:length(Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).ShortName = shortNames{i};
    end
end

%% Parameters

% *** For time course ***
first_time = -100;
last_time = 450;
width = 100;
step = 5;
starts = first_time:step:last_time;
rIntervals = cell(1, length(starts)); % "relevant Intervals"
fname_base = sprintf('%%s_allNeurons_step%d_wd%d.mat', step, width); % double %% escapes the first %s
for i = 1:length(starts)
    rIntervals{i} = [starts(i) (starts(i) + width)];
end

% For simpler analyses
% rIntervals = {[75, 175], [175 225], [175 275], [175 350]};
% fname_base = sprintf('%%s_allNeurons_variableBin_1.mat');


%% Collect spike count data for all neurons

% Parameter for organizing the code. It's important that we decide on how
% the spike counts will be organized, so that we can align the correct
% spike counts with the correct trial information later on!

% for m = 1:length(Monkeys)
for m = 2
%     rSessions = 1:length(Monkeys(m).Sessions); % "relevant sessions"
    rSessions = [4 5];
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        rUnits = 1:length(Monkeys(m).Sessions(sessn).UnitInfo);
        
        % Times images were shown, sorted by category.
        catg1 = Monkeys(m).Sessions(sessn).CueInfo(catgs{1});
        catg2 = Monkeys(m).Sessions(sessn).CueInfo(catgs{2});
        catg1_timeson = vertcat(catg1.Times_on);
        catg2_timeson = vertcat(catg2.Times_on);
            
        % Pre allocate X. X is (num imgs shown) x (num units) x (num
        % intervals). So X(A, B, C) is the number of spikes fired on trial A 
        % (img presentation A) by unit B, in interval C, relative to
        % cue-on.
        X = zeros(length(catg1_timeson)+length(catg2_timeson), length(rUnits), length(rIntervals));
            
        for j = 1:length(rUnits)
            unit = rUnits(j);
            
            % Collect all spike times
            allspike_times = sort(Monkeys(m).Sessions(sessn).UnitInfo(unit).Spike_times);
            
            for p = 1:length(rIntervals)
                interval = rIntervals{p};
                
                % Spikes must be sorted!
                % getSpikesInInt(list of spike times, times catg1 cues on,
                % times catg2 cues on, current interval relative to cue on). 
                X(:,j,p) = getSpikesInInt(allspike_times, catg1_timeson, catg2_timeson, interval);
            end
            fprintf('Done with %s, sessn %d, unit %d \n', Monkeys(m).Name, i, j)
        end
        
        % Generate file name.
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf(fname_base, MonkID);
        
        % Check if file will be too big. If so, split it.
        t = whos('X');
        if t.bytes > 2e9
            % Split into two parts. Split along the interval dimension so
            % we can load one at a time, if necessary later.
            len = round(size(X,3)/2);
            X1 = X(:, :, 1:len);
            X2 = X(:, :, len+1:end);
            save(fullfile(EXT_HD, spikeCountPath, fileName), 'X1', 'X2', 'rIntervals');
        else
            % Save normally.
            save(fullfile(EXT_HD, spikeCountPath, fileName), 'X', 'rIntervals');
        end
    end
end



%% Functions
function vals = getSpikesInInt(sortedSpikeTimes, timesOn1, timesOn2, int)
    
    % collect spike counts for catg1
    catgSC1 = zeros(length(timesOn1), 1);
    for k = 1:length(timesOn1)
        [lower_idx, upper_idx] = binarySearch_window(sortedSpikeTimes, timesOn1(k) + int(1), timesOn1(k) + int(2)); % allspike_times must be sorted!
        if lower_idx == -1 || upper_idx == -1 || upper_idx - lower_idx == -1 % ie no spikes in window
            catgSC1(k) = 0;
            continue
        else
            catgSC1(k) = upper_idx - lower_idx + 1;
        end
    end

     
    % collect spike counts for catg2
    catgSC2 = zeros(length(timesOn2), 1);
    for k = 1:length(timesOn2)
        [lower_idx, upper_idx] = binarySearch_window(sortedSpikeTimes, timesOn2(k) + int(1), timesOn2(k) + int(2)); % allspike_times must be sorted!
        if lower_idx == -1 || upper_idx == -1 || upper_idx - lower_idx == -1 % ie no spikes in window
            catgSC2(k) = 0;
            continue
        else
            catgSC2(k) = upper_idx - lower_idx + 1;
        end
    end
    
    vals = [catgSC1; catgSC2];
end
