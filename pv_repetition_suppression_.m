% repetition suppression analysis

%% Load behavioral data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';

% Use this for main analyses
pv_path = 'XMA2/Monkey_structs';
spikeCountPath = 'XMA2/Spike_count_mats';  % main data
Data = load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_ni.mat'));
catgs = {1:260, 261:520};

fullSVMPath = fullfile(EXT_HD, pv_path, 'SVM_results_%g.mat');

rArrayLocs = {'te'};
% rArrayLocs = {'anterior', 'middle', 'posterior'};
rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)


ignoreVal = 20; % if neuron has less than this num spikes, do not use it.
runShuffle = false; % run the shuffled condition?
    nShuffles = 5;

%% Load the data
% Stitch the sessions back into one big structure
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Repetition suppression analysis

%{
-- sort trials into three types: 1) fixn err, 2) directly after fixn err,
3) all others.
-- using all others, generate mean/std per img/unit combo.
-- use that to z-score responses for fixn err / directly after fixn err
trials.
-- compare distributions.
%}


%% 1. sort trials into types

% ok actually it's a bit more complicated: there are fixation errors that
% can happen before cue onset at all, and there are some trials where there
% is no "fixation error" registered, but there's also nothing else that
% happened, so presumably the monkey never fixated in the first place and
% the trial start timed out or something.
% Let's just start with the trials that are clearly directly after a
% fixation error.

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        
        % Find fixn err trials, on which at least one cue was shown
        fixn_err_idx = find([Monkeys(m).Sessions(sessn).TrialInfo.Fixn_err] == 1 & ...
            cellfun(@(v) sum(~isnan(v))>0, {Monkeys(m).Sessions(sessn).TrialInfo.Cues_seen}));
        
        % TODO: refactor this so it just steps forward from each fixn err
        % trial, and finds the next completed trial, if any. Then we can add info
        % like "imgs that were seen on fe tr" to the completed trial, and
        % also distinguish between trials that directly followed FE's vs
        % those that "indirectly" followed, ie likely had some time in
        % between
        
        % Find completed trials directly following those trials
        following_comp_idx = fixn_err_idx + 1;
        following_comp_idx = following_comp_idx(following_comp_idx <= length(Monkeys(m).Sessions(sessn).TrialInfo));
        following_comp_idx = following_comp_idx([Monkeys(m).Sessions(sessn).TrialInfo(following_comp_idx).Completed]==1);
        [Monkeys(m).Sessions(sessn).TrialInfo.Directly_follows_FE] = deal(0);
        [Monkeys(m).Sessions(sessn).TrialInfo(following_comp_idx).Directly_follows_FE] = deal(1);
        
    end
end

%% Get mean/std per image/unit combo

% takes ~2 minutes
intervals_to_test = {[75 175], [175 275]};

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        
        % Find trials of interest
        fe_bool = [Monkeys(m).Sessions(sessn).TrialInfo.Fixn_err] == 1;
        post_fe_bool = [Monkeys(m).Sessions(sessn).TrialInfo.Directly_follows_FE] == 1;
        other_bool = ([Monkeys(m).Sessions(sessn).TrialInfo.Directly_follows_FE] == 0) & ...
            ([Monkeys(m).Sessions(sessn).TrialInfo.Fixn_err] == 0) & ... 
            ([Monkeys(m).Sessions(sessn).TrialInfo.Completed] == 1);
        
        % Get long lists of paired imgs/times for use later
        % TODO: need to filter this by images that were actually seen in
        % the corresponding fixn err trial (b/c eg sometimes only the 1st
        % of 5 imgs are actually seen, so we would only expect suppression
        % in that one img).
        all_post_fe_imgs = [Monkeys(m).Sessions(sessn).TrialInfo(post_fe_bool).Cues_seen];
        all_post_fe_img_times = [Monkeys(m).Sessions(sessn).TrialInfo(post_fe_bool).Cue_times];
        
        all_other_imgs = [Monkeys(m).Sessions(sessn).TrialInfo(other_bool).Cues_seen];
        all_other_img_times = [Monkeys(m).Sessions(sessn).TrialInfo(other_bool).Cue_times];
        
        % Only bother looking at images that were seen on post-fe trs
        imgs_seen_on_post_fe_trs = unique([Monkeys(m).Sessions(sessn).TrialInfo(post_fe_bool).Cues_seen]);
        imgs_seen_on_post_fe_trs = imgs_seen_on_post_fe_trs(imgs_seen_on_post_fe_trs ~= 0);
        
        for iInterval = 1:length(intervals_to_test)
            interval = intervals_to_test{iInterval};
            
            % pre-allocate matrix for the data
            response_mat = zeros(length(imgs_seen_on_post_fe_trs), length(Monkeys(m).Sessions(sessn).UnitInfo)); % imgs x units x z-score
        
            for iImg = 1:length(imgs_seen_on_post_fe_trs)
                img = imgs_seen_on_post_fe_trs(iImg);

                % Find times when this image was seen either normally, or
                % specifically after a fixation error.
                post_fe_on_times = all_post_fe_img_times(all_post_fe_imgs == img);
                other_on_times = all_other_img_times(all_other_imgs == img);

                for iUnit = 1:length(Monkeys(m).Sessions(sessn).UnitInfo)
                    
                    % Get all spike times for this unit
                    allspike_times = sort(Monkeys(m).Sessions(sessn).UnitInfo(iUnit).Spike_times);
                    
                    % Get this unit's spike counts for these times
                    other_scs = arrayfun(@(t) spikeCountInInterval(allspike_times, t, interval), other_on_times);
                    post_fe_scs = arrayfun(@(t) spikeCountInInterval(allspike_times, t, interval), post_fe_on_times);
                    
                    % Z-score wrt the normal presentations
                    mu = mean(other_scs);
                    sig = std(other_scs);
                    z = mean((post_fe_scs - mu)/sig);
                    
                    % Store the data
                    response_mat(iImg, iUnit) = z;
                end
            end
            
            id = get_good_interval_name2(interval, "", "Post_FE_response_supp_mat");
            Monkeys(m).Sessions(sessn).(id) = response_mat;
            fprintf("Done with %s, session %d, interval %d to %d \n", Monkeys(m).Name, sessn, interval(1), interval(2));
        end
    end
end
        

%% Functions

function sc = spikeCountInInterval(allSortedSpikeTimes, time, interval)
% Given all of a unit's spike times, plus a time of interest and an
% window (interval) around it, calculate the unit's spike count in that window.
% Interval is represented as (a,b) --> (t+a, t+b).
    [lower_idx, upper_idx] = binarySearch_window(allSortedSpikeTimes, time + interval(1), time + interval(2)); % allspike_times must be sorted!
    if lower_idx == -1 || upper_idx == -1 || upper_idx - lower_idx == -1 % ie no spikes in window
        sc = 0;
    else
        sc = upper_idx - lower_idx + 1;
    end
end

