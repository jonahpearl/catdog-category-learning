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
        
        % TODO: try w ignoring "fe interludes"
        
        % Find completed trials directly following those trials
        [Monkeys(m).Sessions(sessn).TrialInfo.Directly_follows_FE] = deal(0);
        [Monkeys(m).Sessions(sessn).TrialInfo.Imgs_seen_on_prev_FE] = deal(0);
        for iFe = 1:length(fixn_err_idx)
            idx = fixn_err_idx(iFe);
            following_completed_idx = idx + 1;
            if (following_completed_idx > length(Monkeys(m).Sessions(sessn).TrialInfo) || ...
                    Monkeys(m).Sessions(sessn).TrialInfo(following_completed_idx).Completed ~= 1)
                continue
            else
                Monkeys(m).Sessions(sessn).TrialInfo(following_completed_idx).Directly_follows_FE = 1;
                Monkeys(m).Sessions(sessn).TrialInfo(following_completed_idx).Imgs_seen_on_prev_FE = Monkeys(m).Sessions(sessn).TrialInfo(idx).Cues_seen;
            end
        end
    end
end

%% Get mean/std per image/unit combo

% takes ~2 minutes
% intervals_to_test = {[75 175], [175 275]};
intervals_to_test = {[175 275]};

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        
        % Find trials of interest
        fe_bool = [Monkeys(m).Sessions(sessn).TrialInfo.Fixn_err] == 1;
        post_fe_bool = [Monkeys(m).Sessions(sessn).TrialInfo.Directly_follows_FE] == 1;
        baseline_bool = ([Monkeys(m).Sessions(sessn).TrialInfo.Directly_follows_FE] == 0) & ...
            ([Monkeys(m).Sessions(sessn).TrialInfo.Fixn_err] == 0) & ... 
            ([Monkeys(m).Sessions(sessn).TrialInfo.Completed] == 1);
        
        % Get long lists of paired imgs/times for use later
        all_post_fe_imgs = [];
        all_post_fe_img_times = [];
        post_fe_idx = find(post_fe_bool==1);
        for iTr = 1:length(post_fe_idx)
            cues_seen = Monkeys(m).Sessions(sessn).TrialInfo(post_fe_idx(iTr)).Cues_seen;
            times = Monkeys(m).Sessions(sessn).TrialInfo(post_fe_idx(iTr)).Cue_times;
            n_cues_on_fe = length(Monkeys(m).Sessions(sessn).TrialInfo(post_fe_idx(iTr)).Imgs_seen_on_prev_FE);
            all_post_fe_imgs = [all_post_fe_imgs cues_seen(1:n_cues_on_fe)];
            all_post_fe_img_times = [all_post_fe_img_times times(1:n_cues_on_fe)];
        end
        
        assert(length(all_post_fe_imgs) == length(all_post_fe_img_times))
        
        all_other_imgs = [Monkeys(m).Sessions(sessn).TrialInfo(baseline_bool).Cues_seen];
        all_other_img_times = [Monkeys(m).Sessions(sessn).TrialInfo(baseline_bool).Cue_times];
        
        % Only bother looking at images that were seen on post-fe trs
        imgs_seen_on_post_fe_trs = unique(all_post_fe_imgs);
        
        for iInterval = 1:length(intervals_to_test)
            interval = intervals_to_test{iInterval};
            
            % pre-allocate matrix for the data
            baseline_mean_responses = zeros(length(imgs_seen_on_post_fe_trs), length(Monkeys(m).Sessions(sessn).UnitInfo)); % imgs x units x z-score
            postfe_mean_responses = zeros(length(imgs_seen_on_post_fe_trs), length(Monkeys(m).Sessions(sessn).UnitInfo)); % imgs x units x z-score
        
            baseline_all_responses = [];
            postfe_all_responses = [];
            
            for iImg = 1:length(imgs_seen_on_post_fe_trs)
                img = imgs_seen_on_post_fe_trs(iImg);

                % Find times when this image was seen either normally, or
                % specifically after a fixation error.
                post_fe_on_times = all_post_fe_img_times(all_post_fe_imgs == img);
                baseline_on_times = all_other_img_times(all_other_imgs == img);

                for iUnit = 1:length(Monkeys(m).Sessions(sessn).UnitInfo)
                    
                    if strcmp(Monkeys(m).Sessions(sessn).UnitInfo(iUnit).Location, "teo")
                        baseline_mean_responses(iImg, iUnit) = NaN;
                        postfe_mean_responses(iImg, iUnit) = NaN;
                        continue
                    end
                    
                    % Get all spike times for this unit
                    allspike_times = sort(Monkeys(m).Sessions(sessn).UnitInfo(iUnit).Spike_times);
                    
                    % Get this unit's spike counts for these times
                    baseline_scs = arrayfun(@(t) spikeCountInInterval(allspike_times, t, interval), baseline_on_times);
                    post_fe_scs = arrayfun(@(t) spikeCountInInterval(allspike_times, t, interval), post_fe_on_times);
                    
                    % Z-score wrt the normal presentations
                    mu = mean(baseline_scs);
                    sig = std(baseline_scs);
                    z_post_fe = mean((post_fe_scs - mu)/sig);
                    z_baseline = mean((baseline_scs - mu)/sig);
                    
                    % Store the data
                    postfe_mean_responses(iImg, iUnit) = z_post_fe;
                    baseline_mean_responses(iImg, iUnit) = z_baseline;
                    
                    baseline_all_responses = [baseline_all_responses (baseline_scs - mu)/sig];
                    postfe_all_responses = [postfe_all_responses (post_fe_scs - mu)/sig];
                end
            end
            
            id = get_good_interval_name2(interval, "", "Post_FE_response_supp_mat");
            Monkeys(m).Sessions(sessn).(id) = postfe_mean_responses;
            
            id = get_good_interval_name2(interval, "", "Post_FE_response_supp_all");
            Monkeys(m).Sessions(sessn).(id) = postfe_all_responses;
            
            id = get_good_interval_name2(interval, "", "Baseline_response_supp_mat");
            Monkeys(m).Sessions(sessn).(id) = baseline_mean_responses;
            
            id = get_good_interval_name2(interval, "", "Baseline_response_supp_all");
            Monkeys(m).Sessions(sessn).(id) = baseline_all_responses;
            
            fprintf("Done with %s, session %d, interval %d to %d \n", Monkeys(m).Name, sessn, interval(1), interval(2));
        end
    end
end

%% Plot results pre vs post

intervals_to_test = {[175 275]};

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    figure
    hold on
    for iInterval = 1:length(intervals_to_test)
        interval = intervals_to_test{iInterval};

        for i = 1:length(rSessions)
            sessn = rSessions(i);

            % pre-allocate matrix for the data
            id = get_good_interval_name2(interval, "", "Post_FE_response_supp_all");
            postfe_mean_responses = Monkeys(m).Sessions(sessn).(id);
            postfe_mean_responses = postfe_mean_responses(~isnan(postfe_mean_responses) & ~isinf(postfe_mean_responses));
            histogram(postfe_mean_responses, 'BinEdges', -4:0.1:6, 'Normalization', 'probability')
            xlabel("Responses on post-FE trials (Z per image/unit pair)")
            ylabel("Probability")
            title(Monkeys(m).Name, "Interpreter", "none")
        end
        legend(["Pre", "Post"])
        
        pre = Monkeys(m).Sessions(rSessions(1)).(id)(:);
        post = Monkeys(m).Sessions(rSessions(2)).(id)(:);
        pre = pre(~isnan(pre) & ~isinf(pre));
        post = post(~isnan(post) & ~isinf(post));
%         [p, h] = ranksum(pre, post)  % ranksum is super sensitive to the long tails

        disp([mean(pre); mean(post)])
        [h, p] = ttest2(pre, post)
    end
end

%% Plot results baseline vs post fe, within sessions

intervals_to_test = {[175 275]};

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    figure
    hold on
    
    for iInterval = 1:length(intervals_to_test)
        interval = intervals_to_test{iInterval};

        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            subplot(2,1,i)
            hold on
            
            % get data
            id = get_good_interval_name2(interval, "", "Post_FE_response_supp_all");
            postfe_mean_responses = Monkeys(m).Sessions(sessn).(id);
            postfe_mean_responses = postfe_mean_responses(~isnan(postfe_mean_responses) & ~isinf(postfe_mean_responses));
            
            id = get_good_interval_name2(interval, "", "Baseline_response_supp_all");
            baseline_mean_responses = Monkeys(m).Sessions(sessn).(id);
            baseline_mean_responses = baseline_mean_responses(~isnan(baseline_mean_responses) & ~isinf(baseline_mean_responses));
            
            % plot
            histogram(postfe_mean_responses, 'BinEdges', -4:0.1:6, 'Normalization', 'probability')
            histogram(baseline_mean_responses, 'BinEdges', -4:0.1:6, 'Normalization', 'probability')
            
            xlabel("Responses on post-FE trials (Z per image/unit pair)")
            ylabel("Probability")
            title(sprintf("Monk %s, %s", Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName), "Interpreter", "none")
            legend(["Post-FE", "Typical"])
            set(gca, "YScale", "log")
            
            fprintf("Monkey %s, %s", Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
%             [p, h] = ranksum(postfe_mean_responses, baseline_mean_responses)
            [h, p] = ttest2(postfe_mean_responses, baseline_mean_responses)
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

