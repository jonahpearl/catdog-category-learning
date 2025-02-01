%% Load behavioral data
clearvars
close all

% Define paths to data
% EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';

% CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
% 
% % Use this for main analyses
% mdRecordPath = 'XMA2/Monkey_structs/MD_Records.mat';
% pv_path = 'XMA2/Monkey_structs';
% behav_file = 'MaxMarta_xma2_behav_and_metaNI.mat';
% spikeCountPath = 'XMA2/Spike_count_mats';
% 
% % Load behavioral data
% load(fullfile(EXT_HD, pv_path, behav_file)) % behavior and neural summaries, but w/o spike times

% interval_data_path = fullfile(EXT_HD, spikeCountPath, fileName);

% Working from locally copied data, run this + drag/drop data .mat files
interval_data_path  = "/Users/jonahpearl/Documents/BJR group/Catdog_paper/20240908_revisions_v2/data/Spike_count_mats/";


%% Set analysis parameters

rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)
% manualIntervals = {[75 175], [175 275]};
manualIntervals = {[175 350]};
rArrayLocs = {'te'};

% path to spike counts. See pv_get_interval_spike_counts if you want to add
% more intervals.
step = 5; % Interval parameters (for loading the correct spike counts file)
width = 100;
% X_fname_base = sprintf('%%s_allNeurons_step%d_wd%d.mat', step, width); % contains full time-course of spike counts (Fig 2b)
X_fname_base = sprintf('%%s_allNeurons_variableBin_1.mat'); % contains 75-175, 175-225, 175-275, and 175-350.
ignoreVal = 20; % if neuron has less than this num spikes, do not use it.

% Other Parameters
TE_LOCS = {'anterior', 'middle', 'posterior'};
all_imgs = 1:520;
catg2_ind1 = 261;
random_seed = 10; % for reproducibility 


paramStruct = struct('RandomSeed', random_seed, ...
    'IgnoreVal', ignoreVal,...
    'SessionsUsed', {rSessionsByMonk},...
    'IntervalsUsed', {manualIntervals});

%% Collect Session Y (image id and catg id)
% very fast

catgs = {1:260, 261:520};
for m = 1:length(Monkeys)
    sessions_to_use = 1:length(Monkeys(m).Sessions);
        
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        catg1 = Monkeys(m).Sessions(sessn).CueInfo(catgs{1});
        catg2 = Monkeys(m).Sessions(sessn).CueInfo(catgs{2});
        catg1_timeson = vertcat(catg1.Times_on);
        catg2_timeson = vertcat(catg2.Times_on);

        % make Y for catg and image ID
        Y = vertcat(repelem(1,length(catg1_timeson))' , repelem(2,length(catg2_timeson))');
        GEO(m).Sessions(sessn).Session_Y_catg = Y;
        
        Y = [];
        for j = 1:length(Monkeys(m).Sessions(sessn).CueInfo)
            Y = vertcat(Y, repelem(Monkeys(m).Sessions(sessn).CueInfo(j).CueID, Monkeys(m).Sessions(sessn).CueInfo(j).NumApp)');
        end
        if length(Y) ~= length(GEO(m).Sessions(sessn).Session_Y_catg)
            error('Error compiling session Y''s ')
        end
        GEO(m).Sessions(sessn).Session_Y_imageID = Y;
    end
end

%% Run the analysis with the real data

% want to use the same num neurons pre vs. post, so we will grab some fraction of the minimum
% lower vals here also seem to make pdist fail less often.
fraction_min_n_units = 0.5;  

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    
    % Determine how many neurons to use based on min across sessions
    % TODO: this will have to change if we want to do array-specific MD
    min_units = inf;
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        n_units = sum(([Monkeys(m).Sessions(sessn).UnitInfo.SpikeNum] > ignoreVal) & ...
            (ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, TE_LOCS))); 
        if n_units < min_units
            min_units = n_units;
        end
    end
    n_units_to_use = round(min_units * fraction_min_n_units);
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Get data from appropriate storage place
        % X_full is (trials) x (units) x (intervals)
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf(X_fname_base, MonkID);
        [X_full, rIntervals] = load_interval_data(fullfile(interval_data_path, fileName)); % X is spike counts, rIntervals is list of each interval
        Y = GEO(m).Sessions(sessn).Session_Y_catg; % list of categories for each trial (image) in X (either 1 or 2)
        imgs_shown = GEO(m).Sessions(sessn).Session_Y_imageID;  % list of imgs for each trial
        
        % For finding intervals in X_full
        rIntervals_original = rIntervals;
        
        for iLoc = 1:length(rArrayLocs)
            array = rArrayLocs{iLoc};
            
            % Find all potential units to be used
            if strcmp(array, 'te')
                viable_units_bool = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, TE_LOCS) & ...
                    [Monkeys(m).Sessions(sessn).UnitInfo.SpikeNum] > ignoreVal;
            else
                viable_units_bool = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, array) & ...
                    [Monkeys(m).Sessions(sessn).UnitInfo.SpikeNum] > ignoreVal;
            end
        
            for iInt = 1:length(manualIntervals)
                interval = manualIntervals{iInt};
                fprintf("Monkey %s, session %s, int %d \n", Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName, iInt)

                % Find index in X's 3rd dim for requested interval
                int_idx = find(cellfun(@(a) all(a == interval), rIntervals_original));

                % Subset to given interval
                X_subset = X_full(:, :, int_idx);  % resulting mat is trials x units
                Y_subset = Y;
                
                % Get results of this bootstrap
                % Perform bootstrapping
                results = calculateDistances(X_subset, Y_subset, imgs_shown, all_imgs, viable_units_bool, n_units_to_use);
                
                % Store results in GEO struct
                md_id = get_good_interval_name2(interval, array, "md_real");
                GEO(m).Sessions(sessn).(md_id) = results;

            end
        end
    end
end

%% (if not loading in from saved) Run the bootstrap

% want to use the same num neurons pre vs. post, so we will grab some fraction of the minimum
% lower vals here also seem to make pdist fail less often.
fraction_min_n_units = 0.5;  
n_bootstrap_iters = 100;

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    
    % Determine how many neurons to use based on min across sessions
    % TODO: this will have to change if we want to do array-specific MD
    min_units = inf;
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        n_units = sum(([Monkeys(m).Sessions(sessn).UnitInfo.SpikeNum] > ignoreVal) & ...
            (ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, TE_LOCS))); 
        if n_units < min_units
            min_units = n_units;
        end
    end
    n_units_to_use = round(min_units * fraction_min_n_units);
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Get data from appropriate storage place
        % X_full is (trials) x (units) x (intervals)
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf(X_fname_base, MonkID);
        [X_full, rIntervals] = load_interval_data(fullfile(EXT_HD, spikeCountPath, fileName)); % X is spike counts, rIntervals is list of each interval
        Y = GEO(m).Sessions(sessn).Session_Y_catg; % list of categories for each trial (image) in X (either 1 or 2)
        imgs_shown = GEO(m).Sessions(sessn).Session_Y_imageID;  % list of imgs for each trial
        
        % For finding intervals in X_full
        rIntervals_original = rIntervals;
        
        for iLoc = 1:length(rArrayLocs)
            array = rArrayLocs{iLoc};
            
            % Find all potential units to be used
            if strcmp(array, 'te')
                viable_units_bool = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, TE_LOCS) & ...
                    [Monkeys(m).Sessions(sessn).UnitInfo.SpikeNum] > ignoreVal;
            else
                viable_units_bool = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, array) & ...
                    [Monkeys(m).Sessions(sessn).UnitInfo.SpikeNum] > ignoreVal;
            end
        
            for iInt = 1:length(manualIntervals)
                interval = manualIntervals{iInt};
                fprintf("Monkey %s, session %s, int %d \n", Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName, iInt)

                % Find index in X's 3rd dim for requested interval
                int_idx = find(cellfun(@(a) all(a == interval), rIntervals_original));

                % Subset to given interval
                X_subset = X_full(:, :, int_idx);  % resulting mat is trials x units
                Y_subset = Y;
                
                % Get results of this bootstrap
                % Perform bootstrapping
                bootstrap_results = cell(n_bootstrap_iters, 1);
                for b = 1:n_bootstrap_iters
                    fprintf("Boot iter %d\n", b);
                    bootstrap_results{b} = calculateDistances_resample(X_subset, Y_subset, imgs_shown, all_imgs, viable_units_bool, n_units_to_use);
                end
                
                % Aggregate bootstrap results
                aggregated_results = struct();
                fields = fieldnames(bootstrap_results{1});
                for f = 1:length(fields)
                    field = fields{f};
                    for b = 1:n_bootstrap_iters
                        aggregated_results(b).(field) = bootstrap_results{b}.(field);
                    end
                end
                
                % Store results in GEO struct
                md_id = get_good_interval_name2(interval, array, "md_bootstrap");
                GEO(m).Sessions(sessn).(md_id) = aggregated_results;

            end
        end
    end
end

%% Save the data

% Save the data and add a row to the Record
fullMDPath = fullfile(EXT_HD, pv_path, 'MD_bootstrap_results_%g.mat');
fullRecordPath = fullfile(EXT_HD, mdRecordPath);
save_SVM_data(GEO, paramStruct, fullMDPath, fullRecordPath);
fprintf('File saved \n')

%% Load the bootstrapped data if needed

ID = 210715;

fullMDPath = fullfile(EXT_HD, pv_path, 'MD_bootstrap_results_%g.mat');
load(sprintf(fullMDPath, ID), 'data');
GEO = data;
clear data

%% Plot the real data, MD w/in categories

figure("Position", [500 500 500 800]) % lin plot
% figure("Position", [500 500 280 460]); % log plot
hold on

xlim_by_monk = {[6 20], [7, 22]};
% outlier_MD_bounds = [5 25];
outlier_MD_bounds = [0 50];

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    disp(Monkeys(m).Name)
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};
        
        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
            
            variances = zeros(2, 2);  % sessions x categories
            no_outlier_variances = zeros(2,2);
            for i = 1:length(rSessions)
                sessn = rSessions(i);
                
                md_id = get_good_interval_name2(interval, array, "md_real");
                results = GEO(m).Sessions(sessn).(md_id);
                within_catg_mds = [results.within_cats(:)' results.within_dogs(:)'];
%                 betwn_catg_mds = results.betw_catgs(:)';
                within_catg_mds_masked = within_catg_mds((within_catg_mds > outlier_MD_bounds(1)) & (within_catg_mds < outlier_MD_bounds(2)));
                
                if contains(Monkeys(m).Sessions(sessn).ShortName, "Pre")
                    pre_within_dists_masked = within_catg_mds;
                elseif contains(Monkeys(m).Sessions(sessn).ShortName, "Post")
                    post_within_dists_masked = within_catg_mds;
                end
                
                subplot(2,1,m)
                hold on
                bins = 0:0.1:50;  % linear yscale
%                 bins = 0:0.5:50;  % log yscale
                histogram(within_catg_mds, 'BinEdges', bins, 'Normalization', 'probability', 'FaceColor', mlc(i))
                
                % Draw mean + std
                mu = mean(within_catg_mds);
                sigma = std(within_catg_mds);
                
                yl = ylim;
                scatter(mu, yl(2) * 1.1, 100, mlc(i), 'filled', 'v',  'MarkerEdgeColor', 'k')
                plot([mu-sigma mu+sigma], repelem(yl(2) * 1.1, 2), '-', 'Color', mlc(i), 'LineWidth', 1.5)
                
                xlim(xlim_by_monk{m})  % to show var redn. on linear-y plots
                formatSVMPlot(gca, gcf)
            end
            
%             legend(["Pre", "Post"])
%             title(sprintf("Monkey %s", Monkeys(m).Name), "Interpreter", "none")
            xlabel("Mahal. dist.")
            ylabel("Prob.")
%             set(gca, 'YScale', 'log')
            title("pre vs post")
            
            [h, p] = vartest2(pre_within_dists_masked, post_within_dists_masked);
            fprintf("Monkey %s, int %d, vartest pre/within vs post/within: p=%0.4g \n", Monkeys(m).Name, iInt, p)
            disp([var(pre_within_dists_masked) var(post_within_dists_masked)]);
            
            [h, p] = ttest2(pre_within_dists_masked, post_within_dists_masked);
            fprintf("Monkey %s, int %d, ttest pre/within vs post/within p=%0.4g \n", Monkeys(m).Name, iInt, p)
            disp([mean(pre_within_dists_masked) mean(post_within_dists_masked)])
            fprintf("")
        end
    end
end

%% Plot the real data, MD betw categories

% figure("Position", [500 500 500 800]) % lin plot
figure("Position", [500 500 280 460]); % log plot
hold on

xlim_by_monk = {[6 20], [7, 22]};
% outlier_MD_bounds = [5 25];
outlier_MD_bounds = [0 50];

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    disp(Monkeys(m).Name)
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};
        
        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
            
            variances = zeros(2, 2);  % sessions x categories
            no_outlier_variances = zeros(2,2);
            for i = 1:length(rSessions)
                sessn = rSessions(i);
                
                md_id = get_good_interval_name2(interval, array, "md_real");
                results = GEO(m).Sessions(sessn).(md_id);
                betw_catg_mds = results.betw_catgs(:)';
                betw_catg_mds_masked = betw_catg_mds((betw_catg_mds > outlier_MD_bounds(1)) & (betw_catg_mds < outlier_MD_bounds(2)));
                
                if contains(Monkeys(m).Sessions(sessn).ShortName, "Pre")
                    pre_betw_dists_masked = betw_catg_mds;
                elseif contains(Monkeys(m).Sessions(sessn).ShortName, "Post")
                    post_betw_dists_masked = betw_catg_mds;
                end
                
                subplot(2,1,m)
                hold on
%                 bins = 0:0.1:50;  % linear yscale
                bins = 0:0.5:50;  % log yscale
                histogram(betw_catg_mds, 'BinEdges', bins, 'Normalization', 'probability', 'FaceColor', mlc(i))
                
                % Draw mean + std
                mu = mean(betw_catg_mds);
                sigma = std(betw_catg_mds);
                yl = ylim;
                scatter(mu, yl(2) * 1.1, 100, mlc(i), 'filled', 'v',  'MarkerEdgeColor', 'k')
                plot([mu-sigma mu+sigma], repelem(yl(2) * 1.1, 2), '-', 'Color', mlc(i), 'LineWidth', 1.5)
                
%                 xlim(xlim_by_monk{m})  % to show var redn. on linear-y plots
                formatSVMPlot(gca, gcf)
            end
            
%             legend(["Pre", "Post"])
%             title(sprintf("Monkey %s", Monkeys(m).Name), "Interpreter", "none")
            xlabel("Mahal. dist.")
            ylabel("Prob.")
            set(gca, 'YScale', 'log')
            title("pre vs post")
            
            [h, p] = vartest2(pre_betw_dists_masked, post_betw_dists_masked);
            fprintf("Monkey %s, int %d, vartest pre/betw vs post/betw: p=%0.4g \n", Monkeys(m).Name, iInt, p)
            disp([var(pre_betw_dists_masked) var(post_betw_dists_masked)]);
            
            [h, p] = ttest2(pre_betw_dists_masked, post_betw_dists_masked);
            fprintf("Monkey %s, int %d, ttest pre/betw vs post/betw p=%0.4g \n", Monkeys(m).Name, iInt, p)
            disp([mean(pre_betw_dists_masked) mean(post_betw_dists_masked)])
            fprintf("")
        end
    end
end

%% Compare bootstrap distbns of (dog - cat) means, pre vs post

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    assert(length(rSessions) == 2);
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
                
            md_id = get_good_interval_name2(interval, array, "md_bootstrap");
            pre_vals = [GEO(m).Sessions(rSessions(1)).(md_id).within_dogs_median] - [GEO(m).Sessions(rSessions(1)).(md_id).within_cats_median];
            post_vals = [GEO(m).Sessions(rSessions(2)).(md_id).within_dogs_median] - [GEO(m).Sessions(rSessions(2)).(md_id).within_cats_median];
            [h,p] = ttest2(pre_vals, post_vals);
            
            fprintf("Monkey %s, pre vs. post (dog - cat): %0.4g \n", Monkeys(m).Name, p)
            fprintf("\tPRE: %0.3g +/- %0.3g \n", mean(pre_vals), std(pre_vals))
            fprintf("\tPOST: %0.3g +/- %0.3g \n", mean(post_vals), std(post_vals))
            
            figure; hold on; [~,edges] = histcounts([pre_vals post_vals]); histogram(pre_vals, 'BinEdges', edges); histogram(post_vals, 'BinEdges', edges)
            title("dogs - cats")
        end
    end
end

%% Compare bootstrap distbns of cat means, pre vs post

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    assert(length(rSessions) == 2);
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
                
            md_id = get_good_interval_name2(interval, array, "md_bootstrap");
            pre_vals = [GEO(m).Sessions(rSessions(1)).(md_id).within_cats_median];
            post_vals = [GEO(m).Sessions(rSessions(2)).(md_id).within_cats_median];
            [h,p] = ttest2(pre_vals, post_vals);
            
            fprintf("Monkey %s, pre vs. post, within-cats MD: %0.4g \n", Monkeys(m).Name, p)
            fprintf("\tPRE: %0.4g +/- %0.3g \n", mean(pre_vals), std(pre_vals))
            fprintf("\tPOST: %0.4g +/- %0.3g \n", mean(post_vals), std(post_vals))
            
            figure; hold on; [~,edges] = histcounts([pre_vals post_vals]); histogram(pre_vals, 'BinEdges', edges); histogram(post_vals, 'BinEdges', edges)
            title("cats")
        end
    end
end

%% Compare bootstrap distbns of dog means, pre vs post

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    assert(length(rSessions) == 2);
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
                
            md_id = get_good_interval_name2(interval, array, "md_bootstrap");
            pre_vals = [GEO(m).Sessions(rSessions(1)).(md_id).within_dogs_median];
            post_vals = [GEO(m).Sessions(rSessions(2)).(md_id).within_dogs_median];
            [h,p] = ttest2(pre_vals, post_vals);
            
            fprintf("Monkey %s, pre vs. post, within-dogs MD: %0.4g \n", Monkeys(m).Name, p)
            fprintf("\tPRE: %0.4g +/- %0.3g \n", mean(pre_vals), std(pre_vals))
            fprintf("\tPOST: %0.4g +/- %0.3g \n", mean(post_vals), std(post_vals))
            
            figure; hold on; [~,edges] = histcounts([pre_vals post_vals]); histogram(pre_vals, 'BinEdges', edges); histogram(post_vals, 'BinEdges', edges)
            title("dogs")
        end
    end
end

%% Compare bootstrap distbns of within-category var, pre vs post

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    assert(length(rSessions) == 2);
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
                
            md_id = get_good_interval_name2(interval, array, "md_bootstrap");
            pre_vals = [GEO(m).Sessions(rSessions(1)).(md_id).within_catg_variance];
            post_vals = [GEO(m).Sessions(rSessions(2)).(md_id).within_catg_variance];
            [h,p] = ttest2(pre_vals, post_vals);
            
            fprintf("Monkey %s, pre vs. post within-category variance: %0.4g \n", Monkeys(m).Name, p)
            fprintf("\tPRE: %0.3g +/- %0.3g \n", mean(pre_vals), std(pre_vals))
            fprintf("\tPOST: %0.3g +/- %0.3g \n", mean(post_vals), std(post_vals))
            
            figure; hold on; [~,edges] = histcounts([pre_vals post_vals]); histogram(pre_vals, 'BinEdges', edges); histogram(post_vals, 'BinEdges', edges)
        end
    end
end

%% Compare bootstrap distbns of (betw - within) mean/median, pre vs post

figure("Position", [500 500 200 350])
for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    assert(length(rSessions) == 2);
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
                
            md_id = get_good_interval_name2(interval, array, "md_bootstrap");
            pre_vals = [GEO(m).Sessions(rSessions(1)).(md_id).betw_catgs_median] - [GEO(m).Sessions(rSessions(1)).(md_id).within_catgs_median];
            post_vals = [GEO(m).Sessions(rSessions(2)).(md_id).betw_catgs_median] - [GEO(m).Sessions(rSessions(2)).(md_id).within_catgs_median];
            [h,p] = ttest2(pre_vals, post_vals);
            
            fprintf("Monkey %s, pre vs. post (betw - within) medians: %0.4g \n", Monkeys(m).Name, p)
            fprintf("\tPRE: %0.3g +/- %0.3g \n", mean(pre_vals), std(pre_vals))
            fprintf("\tPOST: %0.3g +/- %0.3g \n", mean(post_vals), std(post_vals))
            
            n_bins = 20;
            subplot(2,1,m); hold on; [~,edges] = histcounts([pre_vals post_vals], "NumBins", n_bins); histogram(pre_vals, 'BinEdges', edges, "Normalization", "probability"); histogram(post_vals, 'BinEdges', edges, "Normalization", "probability")
            xline(0);
            xlabel("med(b)-med(w)")
            ylabel("Density")
            formatSVMPlot(gca, gcf, 20);
        end
    end
end

%% Compare bootstrap distbns of betw-catg var, pre vs post

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    assert(length(rSessions) == 2);
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
                
            md_id = get_good_interval_name2(interval, array, "md_bootstrap");
            pre_vals = [GEO(m).Sessions(rSessions(1)).(md_id).betw_catg_variance];
            post_vals = [GEO(m).Sessions(rSessions(2)).(md_id).betw_catg_variance];
            [h,p] = ttest2(pre_vals, post_vals);
            
            fprintf("Monkey %s, pre vs. post betw-catg variance: %0.4g \n", Monkeys(m).Name, p)
            fprintf("\tPRE: %0.3g +/- %0.3g \n", mean(pre_vals), std(pre_vals))
            fprintf("\tPOST: %0.3g +/- %0.3g \n", mean(post_vals), std(post_vals))
            
            figure; hold on; [~,edges] = histcounts([pre_vals post_vals]); histogram(pre_vals, 'BinEdges', edges); histogram(post_vals, 'BinEdges', edges)
        end
    end
end

%% Compare bootstrap paired (within var - betw var) pre vs post

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    assert(length(rSessions) == 2);
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
                
            md_id = get_good_interval_name2(interval, array, "md_bootstrap");
            pre_vals = ([GEO(m).Sessions(rSessions(1)).(md_id).within_catg_variance] - ...
                [GEO(m).Sessions(rSessions(1)).(md_id).betw_catg_variance]);
            
            post_vals = ([GEO(m).Sessions(rSessions(2)).(md_id).within_catg_variance] - ...
                [GEO(m).Sessions(rSessions(2)).(md_id).betw_catg_variance]);
            
            [h,p] = ttest2(pre_vals, post_vals);
            
            fprintf("Monkey %s, pre vs. post (within-var - betw-var): %0.4g \n", Monkeys(m).Name, p)
            fprintf("\tPRE: %0.3g +/- %0.3g \n", mean(pre_vals), std(pre_vals))
            fprintf("\tPOST: %0.3g +/- %0.3g \n", mean(post_vals), std(post_vals))
            
            figure; hold on; [~,edges] = histcounts([pre_vals post_vals]); histogram(pre_vals, 'BinEdges', edges); histogram(post_vals, 'BinEdges', edges)
        end
    end
end

%% Compare bootstrap distbns of JSD, pre vs post

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    assert(length(rSessions) == 2);
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
                
            md_id = get_good_interval_name2(interval, array, "md_bootstrap");
            pre_vals = [GEO(m).Sessions(rSessions(1)).(md_id).jsd_within_vs_betw];
            post_vals = [GEO(m).Sessions(rSessions(2)).(md_id).jsd_within_vs_betw];
            [h,p] = ttest2(pre_vals, post_vals);
            
            fprintf("Monkey %s, pre vs. post, MD JSD(w/in vs betw): %0.4g \n", Monkeys(m).Name, p)
            fprintf("\tPRE: %0.3g +/- %0.3g \n", mean(pre_vals), std(pre_vals))
            fprintf("\tPOST: %0.3g +/- %0.3g \n", mean(post_vals), std(post_vals))
            
            figure; hold on; [~,edges] = histcounts([pre_vals post_vals]); histogram(pre_vals, 'BinEdges', edges); histogram(post_vals, 'BinEdges', edges)
        end
    end
end

%% Compare bootstrap distbns of dprime, pre vs post

figure("Position", [500 500 200 350])
for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    assert(length(rSessions) == 2);
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
                
            md_id = get_good_interval_name2(interval, array, "md_bootstrap");
            pre_vals = [GEO(m).Sessions(rSessions(1)).(md_id).dprime_within_vs_betw];
            post_vals = [GEO(m).Sessions(rSessions(2)).(md_id).dprime_within_vs_betw];
            [h,p] = ttest2(pre_vals, post_vals);
            
            % computeCohen_d(x1, x2) --> if mean(x1) > mean(x2), then d' >
            % 0, and vice versa.
%             results.dprime_within_vs_betw = computeCohen_d(MD_within_all, dists_betw_catgs);
            % so negative d' --> mean(betw) > mean(within).
            fprintf("Monkey %s, pre vs. post, MD d'(w/in vs betw): %0.4g \n", Monkeys(m).Name, p)
            fprintf("\tPRE: %0.3g +/- %0.3g \n", mean(pre_vals), std(pre_vals))
            fprintf("\tPOST: %0.3g +/- %0.3g \n", mean(post_vals), std(post_vals))
            
            subplot(2,1,m); hold on; [~,edges] = histcounts([pre_vals post_vals], "NumBins", 25); histogram(pre_vals, 'BinEdges', edges, "Normalization", "probability"); histogram(post_vals, 'BinEdges', edges, "Normalization", "probability")
            xlabel("D'(w vs b)")
            ylabel("Probability")
            formatSVMPlot(gca, gcf, 20)
        end
    end
end

%% Functions

function results = calculateDistances(X_subset, Y, imgs_shown, all_imgs, viable_units_bool, n_units_to_use)
    results = struct();
    jsd_epsilon = 0;
    catg2_ind1 = 261;
    
    done = false;
    n_attempts = 0;
    while (~done) && (n_attempts < 200)
        try
            % Select random subset of units
            units_to_use = randsample(find(viable_units_bool), n_units_to_use);

            X_train = X_subset(:, units_to_use);
            Y_train = Y(:);

            % Check condition number
            condition_num = cond(X_train(Y_train==1, :));
            if condition_num > 300
                n_attempts = n_attempts + 1;
                continue
            end

            % Get the pairwise dists
            all_dists = squareform(pdist(X_train, 'mahalanobis'));
            cat_inds = 1:sum(Y_train==1);
            dog_inds = (sum(Y_train==1)+1):(length(Y_train));
            cat_dists = all_dists(cat_inds, cat_inds);
            dog_dists = all_dists(dog_inds, dog_inds);
            dists_within_cats = upperTriuVals(cat_dists);
            dists_within_dogs = upperTriuVals(dog_dists);

            % Ensure the only zeros are on the diagonal
            assert(sum(sum(all_dists == 0)) == size(all_dists,1))
            assert(sum(sum(cat_dists == 0)) == size(cat_dists,1))
            assert(sum(sum(dog_dists == 0)) == size(dog_dists,1))
            
            % Filter out by stimulus
            dists_betw_catgs = all_dists(cat_inds, dog_inds);
            MD_within_all = [dists_within_cats; dists_within_dogs];

            % Store the full distributions (wont do this for bootstraps)
            results.within_cats = dists_within_cats;
            results.within_dogs = dists_within_dogs;
            results.betw_catgs = dists_betw_catgs;
            results.within_catgs = MD_within_all;
            
            results.within_cats_mean = mean(dists_within_cats);
            results.within_cats_median = median(dists_within_cats);
            results.within_cats_var = var(dists_within_cats);
            results.within_dogs_mean = mean(dists_within_dogs);
            results.within_dogs_median = median(dists_within_dogs);
            results.within_dogs_var = var(dists_within_dogs);
            results.betw_catgs_mean = mean(dists_betw_catgs);
            results.betw_catgs_median = median(dists_betw_catgs);
            results.betw_catg_variance = var(dists_betw_catgs);
            
            % Calculate additional metrics
            results.within_catgs_mean = mean(MD_within_all);
            results.within_catgs_median = median(MD_within_all);
            results.within_catg_variance = var(MD_within_all);

            % Calculate JSD
            [within_counts, within_edges] = histcounts(MD_within_all, 'Normalization', 'probability');
            [betw_counts, ~] = histcounts(dists_betw_catgs, within_edges, 'Normalization', 'probability');
            results.jsd_within_vs_betw = JSDiv(within_counts + jsd_epsilon, betw_counts + jsd_epsilon);

            % Calculate d-prime
            % computeCohen_d(x1, x2) --> if mean(x1) > mean(x2), then d' >
            % 0.
            results.dprime_within_vs_betw = computeCohen_d(MD_within_all, dists_betw_catgs);

            done = true;
        catch
            n_attempts = n_attempts + 1;
        end
    end
    
    if n_attempts == 200
        warning('Unable to find good set of units for this data');
        results = [];
    end
end


function results = calculateDistances_resample(X_subset, Y, imgs_shown, all_imgs, viable_units_bool, n_units_to_use)
    results = struct();
    jsd_epsilon = 1e-12;
    
    done = false;
    n_attempts = 0;
    while (~done) && (n_attempts < 200)
        try
            % Select random subset of units
            units_to_use = randsample(find(viable_units_bool), n_units_to_use);

            X_train = X_subset(:, units_to_use);
            Y_train = Y(:);

            % Check condition number
            condition_num = cond(X_train(Y_train==1, :));
            if condition_num > 300
                n_attempts = n_attempts + 1;
                continue
            end

            % Get the pairwise dists
            all_dists = squareform(pdist(X_train, 'mahalanobis'));

            % Ensure the only zeros are on the diagonal
            assert(sum(sum(all_dists == 0)) == size(all_dists,1))
            
            % Get results with resampling of trials
            cat_trials_to_use = randsample(1:sum(Y_train==1), sum(Y_train==1), true);  % ie with replacement
            dog_trials_to_use = randsample(((sum(Y_train==1)+1):length(Y_train)), sum(Y_train==2), true);  % ie with replacement
            cat_inds_to_use = nchoosek(cat_trials_to_use, 2);
            dog_inds_to_use = nchoosek(dog_trials_to_use, 2);
            
            dists_within_cats = zeros(size(cat_inds_to_use, 1), 1);
            for ii = 1:length(dists_within_cats)
               dists_within_cats(ii) = all_dists(cat_inds_to_use(ii,1), cat_inds_to_use(ii,2));
            end
            
            dists_within_dogs = zeros(size(dog_inds_to_use, 1), 1);
            for ii = 1:length(dists_within_dogs)
               dists_within_dogs(ii) = all_dists(dog_inds_to_use(ii,1), dog_inds_to_use(ii,2));
            end
            
            dists_betw_catgs = zeros(length(cat_trials_to_use) * length(dog_trials_to_use),1);
            ii = 1;
            for cat_tr = cat_trials_to_use
                for dog_tr = dog_trials_to_use
                    dists_betw_catgs(ii) = all_dists(cat_tr, dog_tr);
                    ii = ii + 1;
                end
            end
            
            
            results.within_cats_mean = mean(dists_within_cats);
            results.within_cats_median = median(dists_within_cats);
            results.within_cats_var = var(dists_within_cats);
            results.within_dogs_mean = mean(dists_within_dogs);
            results.within_dogs_median = median(dists_within_dogs);
            results.within_dogs_var = var(dists_within_dogs);
            results.betw_catgs_mean = mean(dists_betw_catgs);
            results.betw_catgs_median = median(dists_betw_catgs);
            results.betw_catg_variance = var(dists_betw_catgs);
            
            % Calculate additional metrics
            MD_within_all = [dists_within_cats; dists_within_dogs];
            results.within_catgs_mean = mean(MD_within_all);
            results.within_catgs_median = median(MD_within_all);
            results.within_catg_variance = var(MD_within_all);

            % Calculate JSD
            [within_counts, within_edges] = histcounts(MD_within_all, 'Normalization', 'probability');
            [betw_counts, ~] = histcounts(dists_betw_catgs, within_edges, 'Normalization', 'probability');
            results.jsd_within_vs_betw = JSDiv(within_counts + jsd_epsilon, betw_counts + jsd_epsilon);

            % Calculate d-prime
            % computeCohen_d(x1, x2) --> if mean(x1) > mean(x2), then d' >
            % 0.
            results.dprime_within_vs_betw = computeCohen_d(MD_within_all, dists_betw_catgs);

            done = true;
        catch
            n_attempts = n_attempts + 1;
        end
    end
    
    if n_attempts == 200
        warning('Unable to find good set of units for this data');
        results = [];
    end
end


function vals = upperTriuVals(x)
    vals = triu(x);
    vals = vals(vals > 0);
end

function pdist_inds = getPdistInds(n)
    idx1 = [];  % Vector to hold the first index of each pair
    idx2 = [];  % Vector to hold the second index of each pair

    for i = 1:n-1
        idx1 = [idx1; i*ones(n-i, 1)];  % Repeat i for 'n-i' times
        idx2 = [idx2; (i+1:n)'];        % Indices from i+1 to n
    end
    pdist_inds = [idx1 idx2];
end


function formatSVMPlot(ax, fig, fontsize)
if nargin == 2
    fontsize=28;
end
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize', fontsize, 'YGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end
