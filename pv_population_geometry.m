% population response geometry across learning

%% Load behavioral data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';

% Use this for main analyses
svmRecordPath = 'XMA2/Monkey_structs/SVM_Records.mat';
pv_path = 'XMA2/Monkey_structs';
behav_file = 'MaxMarta_xma2_behav_and_metaNI.mat';
spikeCountPath = 'XMA2/Spike_count_mats';

% Load behavioral data
load(fullfile(EXT_HD, pv_path, behav_file)) % behavior and neural summaries, but w/o spike times


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


%% Mahalanobis distances wrt distributions
% here we use MATLAB's mahal function to compute the distance of a set of
% images to a distribution (ie if we had 100 cat trials, calculate
% the mahal dists of a held-out set of 20 based on the covaraince of the
% remaining 80).

% the other way to do this would be to use the entire set of category trials
% except two, and calc pairwise mahal dists.

n_folds = 5;
fraction_min_n_units = 0.85;  % want to use the same num neurons pre vs. post, so we will grab some fraction of the minimum

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

                % Find index in X's 3rd dim for requested interval
                int_idx = find(cellfun(@(a) all(a == interval), rIntervals_original));

                % Subset to given interval
                X_subset = X_full(:, :, int_idx);  % resulting mat is trials x units
                Y_subset = Y;
                
                % Split data into train/test folds
                num_cues = 520; % how many cues total exist
                balancedKFoldList = zeros(1, num_cues); % for each cue being used, what fold is it in?
                balancedKFoldList(1:(num_cues/2)) = makeKFoldList(num_cues/2, n_folds);  % make folds separately for cats/dogs to assure even catg splits
                balancedKFoldList((num_cues/2 + 1):end) = makeKFoldList(num_cues/2, n_folds);

                MD_within_cats = [];
                MD_within_dogs = [];
                MD_cats_wrt_dogs = [];
                MD_dogs_wrt_cats = [];
                
                for k = 1:n_folds 
                    
                    
                    % Select random subset of units for this fold
                    units_to_use = randsample(find(viable_units_bool), n_units_to_use);
                    
                    % Select training imgs
                    trainingImgs = all_imgs(balancedKFoldList ~= k); 
                    idx = ismember(imgs_shown, trainingImgs); 
                    
                    X_train = X_subset(idx, units_to_use);
                    Y_train = Y(idx);
                    X_test = X_subset(~idx, units_to_use);
                    Y_test = Y(~idx);
                    
                    
                    % Call signature: MD = mahal(test samples, reference samples)
                    
                    % Within cats
                    MD_within_cats = [MD_within_cats; sqrt(mahal(X_test(Y_test == 1, :), X_train(Y_train == 1, :)))];
                    
                    % Within dogs
                    MD_within_dogs =[MD_within_dogs; sqrt(mahal(X_test(Y_test == 2, :), X_train(Y_train == 2, :)))];
                    
                    % Cats wrt dogs
                    MD_cats_wrt_dogs = [MD_cats_wrt_dogs; sqrt(mahal(X_test(Y_test == 1, :), X_train(Y_train == 2, :)))];
                    
                    % Dogs wrt cats
                    MD_dogs_wrt_cats = [MD_dogs_wrt_cats; sqrt(mahal(X_test(Y_test == 2, :), X_train(Y_train == 1, :)))];
                    
                end
                
                id = get_good_interval_name2(interval, array, "MD_within_cats");
                GEO(m).Sessions(sessn).(id) = MD_within_cats;
                
                id = get_good_interval_name2(interval, array, "MD_within_dogs");
                GEO(m).Sessions(sessn).(id) = MD_within_dogs;
                
                id = get_good_interval_name2(interval, array, "MD_cats_wrt_dogs");
                GEO(m).Sessions(sessn).(id) = MD_cats_wrt_dogs;
                
                id = get_good_interval_name2(interval, array, "MD_dogs_wrt_cats");
                GEO(m).Sessions(sessn).(id) = MD_dogs_wrt_cats;
                
            end
        end
    end
end


%% Plot mahalnobis results

f1 = figure('Position', [500, 700, 265, 400]);
hold on

f2 = figure('Position', [200, 300, 265, 400]);
hold on

outlier_MD_bounds = [5 25];

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    disp(Monkeys(m).Name)
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};
        
        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};

            for i = 1:length(rSessions)
                sessn = rSessions(i);
                
                within_cats_id = get_good_interval_name2(interval, array, "MD_within_cats");
                MD_within_cats = GEO(m).Sessions(sessn).(within_cats_id);
                MD_within_cats = MD_within_cats(~isnan(MD_within_cats));
                MD_within_cats = MD_within_cats((MD_within_cats > outlier_MD_bounds(1)) & (MD_within_cats  < outlier_MD_bounds(2)));
                
                within_dogs_id = get_good_interval_name2(interval, array, "MD_within_dogs");
                MD_within_dogs= GEO(m).Sessions(sessn).(within_dogs_id);
                MD_within_dogs = MD_within_dogs(~isnan(MD_within_dogs));
                MD_within_dogs = MD_within_dogs((MD_within_dogs > outlier_MD_bounds(1)) & (MD_within_dogs  < outlier_MD_bounds(2)));
                
                id = get_good_interval_name2(interval, array, "MD_cats_wrt_dogs");
                MD_cats_wrt_dogs = GEO(m).Sessions(sessn).(id);
                MD_cats_wrt_dogs = MD_cats_wrt_dogs(~isnan(MD_cats_wrt_dogs));
                MD_cats_wrt_dogs = MD_cats_wrt_dogs((MD_cats_wrt_dogs > outlier_MD_bounds(1)) & (MD_cats_wrt_dogs < outlier_MD_bounds(2)));
                
                id = get_good_interval_name2(interval, array, "MD_dogs_wrt_cats");
                MD_dogs_wrt_cats = GEO(m).Sessions(sessn).(id);
                MD_dogs_wrt_cats = MD_dogs_wrt_cats(~isnan(MD_dogs_wrt_cats));
                MD_dogs_wrt_cats = MD_dogs_wrt_cats((MD_dogs_wrt_cats > outlier_MD_bounds(1)) & (MD_dogs_wrt_cats < outlier_MD_bounds(2)));
                
%                 disp(size(MD_within_dogs))
                
                % Plot within catg dists in f1
                figure(f1)
                subplot(2,1,m)
                hold on
                errorbar(i, mean(MD_within_cats), std(MD_within_cats)/sqrt(length(MD_within_cats)),...
                    "o", 'LineWidth', 2, 'Color', 'r') 
                errorbar(i+0.25, mean(MD_within_dogs), std(MD_within_dogs)/sqrt(length(MD_within_dogs)),...
                    "o", 'LineWidth', 2, 'Color', 'b')
                
                % Plot within vs betw in f2
                all_within = [MD_within_dogs; MD_within_cats];
                all_betw = [MD_cats_wrt_dogs; MD_dogs_wrt_cats];
                GEO(m).Sessions(sessn).Mahal_within_vs_betw = {all_within, all_betw};
                
                figure(f2)
                subplot(2,1,m)
                hold on
                errorbar(i, mean(all_within), std(all_within)/sqrt(length(all_within)),...
                    "o", 'LineWidth', 2, 'Color', 'k') 
                errorbar(i+0.25, mean(all_betw), std(all_betw)/sqrt(length(all_betw)),...
                    "o", 'LineWidth', 2, 'Color', [0.5, 0.5, 0.5])
            end
            
            figure(f1)
%             title("Cat vs dog")
            xlim([0.5, 2.5])
            formatSVMPlot(gca, gcf);
            xlabel("Session")
            ylabel("Mahal dist.")
            
            figure(f2)
%             title("Within vs across")
            formatSVMPlot(gca, gcf);
            xlabel("Session")
            ylabel("Mahal dist.")
            xlim([0.5, 2.5])
            
            
            
            % Do some stats
            [h, p] = ttest2(GEO(m).Sessions(rSessions(1)).(within_cats_id), ...
                            GEO(m).Sessions(rSessions(1)).(within_dogs_id));
            fprintf("Monkey %s, cat vs dog PRE: p = %0.4g \n", Monkeys(m).Name, p)           
             
            [h, p] = ttest2(GEO(m).Sessions(rSessions(2)).(within_cats_id), ...
                             GEO(m).Sessions(rSessions(2)).(within_dogs_id));
            fprintf("Monkey %s, cat vs dog POST: p = %0.4g \n", Monkeys(m).Name, p)           
                        
            [h, p] = ttest2(GEO(m).Sessions(rSessions(1)).(within_cats_id), ...
                            GEO(m).Sessions(rSessions(2)).(within_cats_id));
            fprintf("Monkey %s, cat pre vs cat post: p = %0.4g \n", Monkeys(m).Name, p)
            
            [h, p] = ttest2(GEO(m).Sessions(rSessions(1)).(within_dogs_id), ...
                            GEO(m).Sessions(rSessions(2)).(within_dogs_id));
            fprintf("Monkey %s, dog pre vs dog post: p = %0.4g \n", Monkeys(m).Name, p)
            
            [h, p] = ttest2(GEO(m).Sessions(rSessions(1)).Mahal_within_vs_betw{1}, ...
                            GEO(m).Sessions(rSessions(1)).Mahal_within_vs_betw{2});
            fprintf("Monkey %s, within vs betw PRE: p = %0.4g \n", Monkeys(m).Name, p)
            
            [h, p] = ttest2(GEO(m).Sessions(rSessions(2)).Mahal_within_vs_betw{1}, ...
                            GEO(m).Sessions(rSessions(2)).Mahal_within_vs_betw{2});
            fprintf("Monkey %s, within vs betw POST: p = %0.4g \n", Monkeys(m).Name, p)
            
%             [h, p] = ttest2(GEO(m).Sessions(rSessions(1)).Mahal_within_vs_betw{1}, ...
%                             GEO(m).Sessions(rSessions(2)).Mahal_within_vs_betw{2});
%             fprintf("Monkey %s, within, PRE: p = %0.4g \n", Monkeys(m).Name, p)
            
            
             
        end
        
    end
end


%% Pairewise dist metrics take 2 (random unit subset folds)

n_folds = 5;
fraction_min_n_units = 0.5;  % want to use the same num neurons pre vs. post, so we will grab some fraction of the minimum

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
                
                % Split data into train/test folds
                num_cues = 520; % how many cues total exist
                balancedKFoldList = zeros(1, num_cues); % for each cue being used, what fold is it in?
                balancedKFoldList(1:(num_cues/2)) = makeKFoldList(num_cues/2, n_folds);  % make folds separately for cats/dogs to assure even catg splits
                balancedKFoldList((num_cues/2 + 1):end) = makeKFoldList(num_cues/2, n_folds);

                MD_within_cats = [];
                all_cat_img_pair_inds = [];
                
                MD_within_dogs = [];
                all_dog_img_pair_inds = [];
                
                MD_betw_catgs = [];
                
%                 MD_cats_wrt_dogs = [];
%                 MD_dogs_wrt_cats = [];
                
                for k = 1:n_folds 
%                     fprintf("k=%d \n", k)
                    
                    done = false;
                    n_attempts = 0;
                    while (~done) && (n_attempts < 100)
                        try
                            
                            % Select random subset of units for this fold
                            units_to_use = randsample(find(viable_units_bool), n_units_to_use);

                            % Select training imgs
                            trainingImgs = all_imgs(balancedKFoldList ~= k); 
                            
                            % Post-hoc: Exclude imgs 307/312 for Max
                            % because they are massive outliers. Check if
                            % this changes the results --> no.
                            if m == 2
                                trainingImgs = trainingImgs(~ismember(trainingImgs, [307, 312]));
                            end
                            
                            % Subset trials based on training images
                            idx = ismember(imgs_shown, trainingImgs); 
                            X_train = X_subset(idx, units_to_use);
                            Y_train = Y(idx);

                            % Some combos of units give poorly conditioned
                            % covariance matrices, try randomly until we
                            % find one that will calculate successfully.
%                             disp(cond(X_train(Y_train==1, :)))
                            condition_num = cond(X_train(Y_train==1, :));
                            if condition_num > 300
                                n_attempts = n_attempts + 1;
%                                 fprintf("%d \n ", n_attempts)
                                continue
                            end
                            
                            % Get the pairwise dists
                            all_dists = squareform(pdist(X_train, 'mahalanobis'));
                            cat_dists = squareform(pdist(X_train(Y_train==1, :), 'mahalanobis'));
                            dog_dists = squareform(pdist(X_train(Y_train==2, :), 'mahalanobis'));
                            
                            done = true;
                        catch
                            n_attempts = n_attempts + 1;
%                             fprintf("%d \n ", n_attempts)
                        end
                    end

                    % Ensure the only zeros are on the diagonal
                    assert(sum(sum(all_dists == 0)) == size(all_dists,1))
                    assert(sum(sum(cat_dists == 0)) == size(cat_dists,1))
                    assert(sum(sum(dog_dists == 0)) == size(dog_dists,1))

                    % Filter out by stimulus
                    dists_betw_catgs = all_dists(1:catg2_ind1, catg2_ind1:end);
                    dists_within_cats = upperTriuVals(cat_dists);
                    dists_within_dogs = upperTriuVals(dog_dists);
                    
                    % ID which dists are coming from which img pairs so we
                    % can check later. We do this by creating a vector that
                    % follows pdist's ordering (getPdistInds), and then
                    % using that to figure out which images correspond to
                    % which dists. It's only so many lines because you
                    % can't do f-ing multiple indexing ops in a row in
                    % MATLAB.
                    imgs_used = imgs_shown(idx);
                    imgs_used_cats = imgs_used(Y_train==1);
                    cat_img_pair_inds = getPdistInds(sum(Y_train==1));
                    cat_img_pair_inds = imgs_used_cats(cat_img_pair_inds);
                    imgs_used_dogs = imgs_used(Y_train==2);
                    dog_img_pair_inds = getPdistInds(sum(Y_train==2));
                    dog_img_pair_inds = imgs_used_dogs(dog_img_pair_inds);
                    assert(length(cat_img_pair_inds) == length(dists_within_cats))
                    assert(length(dog_img_pair_inds) == length(dists_within_dogs))
                    
                    MD_within_cats = [MD_within_cats; dists_within_cats];
                    all_cat_img_pair_inds = [all_cat_img_pair_inds; cat_img_pair_inds];
                    
                    MD_within_dogs = [MD_within_dogs; dists_within_dogs];
                    all_dog_img_pair_inds = [all_dog_img_pair_inds; dog_img_pair_inds];
                    
                    MD_betw_catgs = [MD_betw_catgs(:); dists_betw_catgs(:)];

                end
                id = get_good_interval_name2(interval, array, "m_pairwise_dists_betw_catgs");
                GEO(m).Sessions(sessn).(id) = MD_betw_catgs;
                
                id = get_good_interval_name2(interval, array, "m_pairwise_dists_within_cats");
                GEO(m).Sessions(sessn).(id) = MD_within_cats;
                
                id = get_good_interval_name2(interval, array, "m_pairwise_dists_img_pairs_within_cats");
                GEO(m).Sessions(sessn).(id) = all_cat_img_pair_inds;

                id = get_good_interval_name2(interval, array, "m_pairwise_dists_within_dogs");
                GEO(m).Sessions(sessn).(id) = MD_within_dogs;
                
                id = get_good_interval_name2(interval, array, "m_pairwise_dists_img_pairs_within_dogs");
                GEO(m).Sessions(sessn).(id) = all_dog_img_pair_inds;
            end
        end
    end
end

%% Plot pairwise dist results, within catg

figure("Position", [500 500 1200 800])
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
                
                cat_id = get_good_interval_name2(interval, array, "m_pairwise_dists_within_cats");
                m_dists_within_cats = GEO(m).Sessions(sessn).(cat_id);
                variances(i, 1) = var(m_dists_within_cats);
                cat_mask = (m_dists_within_cats > outlier_MD_bounds(1)) & (m_dists_within_cats < outlier_MD_bounds(2));
                no_outlier_variances(i,1) = var(m_dists_within_cats(cat_mask));
                
                dog_id = get_good_interval_name2(interval, array, "m_pairwise_dists_within_dogs");
                m_dists_within_dogs = GEO(m).Sessions(sessn).(dog_id);
                variances(i, 2) = var(m_dists_within_dogs);
                dog_mask = (m_dists_within_dogs > outlier_MD_bounds(1)) & (m_dists_within_dogs < outlier_MD_bounds(2));
                no_outlier_variances(i,2) = var(m_dists_within_dogs(dog_mask));
                
                all_within_dists = [m_dists_within_cats; m_dists_within_dogs];
                all_within_dists_masked = [m_dists_within_cats(cat_mask); m_dists_within_dogs(dog_mask)];
                
                if contains(Monkeys(m).Sessions(sessn).ShortName, "Pre")
                    pre_within_dists_masked = all_within_dists_masked;
                elseif contains(Monkeys(m).Sessions(sessn).ShortName, "Post")
                    post_within_dists_masked = all_within_dists_masked;
                end
                
                % Try removing imgs 307/312 which are outliers in mk 2 -->
                % makes no diff.
%                 if m == 2
%                     id = get_good_interval_name2(interval, array, "m_pairwise_dists_img_pairs_within_dogs");
%                     all_dog_img_pair_inds = GEO(m).Sessions(sessn).(id);
%                     outlier_bool = (any((all_dog_img_pair_inds == 312),2) | any((all_dog_img_pair_inds == 307),2));
%                     m_dists_within_dogs = m_dists_within_dogs(~outlier_bool);
%                 end
                
%                 subplot(1,2,1)
%                 hold on
%                 errorbar(i, mean(m_dists_within_cats(:)), std(m_dists_within_cats(:))/sqrt(length(m_dists_within_cats)),...
%                     "o", 'LineWidth', 2, 'Color', mlc(1)) 
%                 errorbar(i+0.25, mean(m_dists_within_dogs(:)), std(m_dists_within_dogs(:))/sqrt(length(m_dists_within_dogs)),...
%                     "o", 'LineWidth', 2, 'Color', mlc(2))
                

                % Plot cats/dogs separately --> they're the same
%                 subplot(2,1,1)
%                 hold on
%                 histogram(m_dists_within_cats, 'BinEdges', 5:0.1:35, 'Normalization', 'probability')
%                 title("Cats")
%                 set(gca, 'YScale', 'log')
%                 
%                 subplot(2,1,2)
%                 hold on
%                 histogram(m_dists_within_dogs, 'BinEdges', 5:0.1:35, 'Normalization', 'probability')
%                 title("Dogs")
%                 set(gca, 'YScale', 'log')

                % Concat all pairs
                subplot(2,1,m)
                hold on
%                 linear_bins = 0:0.1:35;
                linear_bins = 0:0.1:50;
                histogram(all_within_dists, 'BinEdges', linear_bins, 'Normalization', 'probability', 'FaceColor', mlc(i))
                
                % Draw mean + std
                mu = mean(all_within_dists_masked);
                sigma = std(all_within_dists_masked);
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
            
            
            [h, p] = vartest2(pre_within_dists_masked, post_within_dists_masked);
            fprintf("Monkey %s, int %d, vartest pre/within vs post/within: p=%0.4g \n", Monkeys(m).Name, iInt, p)
            disp([var(pre_within_dists_masked) var(post_within_dists_masked)]);
            
            [h, p] = ttest2(pre_within_dists_masked, post_within_dists_masked);
            fprintf("Monkey %s, int %d, ttest pre/within vs post/within p=%0.4g \n", Monkeys(m).Name, iInt, p)
            disp([mean(pre_within_dists_masked) mean(post_within_dists_masked)])
            
            histcounts_pre = histcounts(pre_within_dists_masked, 'BinEdges', linear_bins, 'Normalization', 'probability');
            histcounts_post = histcounts(post_within_dists_masked, 'BinEdges', linear_bins, 'Normalization', 'probability');
            jsd = JSDiv(histcounts_pre, histcounts_post);
            disp(jsd)
        end
    end
end

%% Plot pairwise dist results, between vs within catg

xlim_by_monk = {[6 20], [7, 22]};

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    disp(Monkeys(m).Name)
    
    figure
    hold on
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};
        
        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
            
            variances = zeros(2, 2);  % sessions x categories
            no_outlier_variances = zeros(2,2);
            for i = 1:length(rSessions)
                sessn = rSessions(i);
                
                cat_id = get_good_interval_name2(interval, array, "m_pairwise_dists_within_cats");
                m_dists_within_cats = GEO(m).Sessions(sessn).(cat_id);
                
                dog_id = get_good_interval_name2(interval, array, "m_pairwise_dists_within_dogs");
                m_dists_within_dogs = GEO(m).Sessions(sessn).(dog_id);
                
                betw_id = get_good_interval_name2(interval, array, "m_pairwise_dists_betw_catgs");
                m_dists_betw_catgs = GEO(m).Sessions(sessn).(betw_id);
                
                all_within = [m_dists_within_cats; m_dists_within_dogs];
                
                % Concat all pairs
                subplot(2,1,i)
                hold on
                
%                 linear_bins = 0:0.1:35;
                linear_bins = 0:0.1:50;
                histogram(all_within, 'BinEdges', linear_bins, 'Normalization', 'probability')
                histogram(m_dists_betw_catgs, 'BinEdges', linear_bins, 'Normalization', 'probability')
                
                [h, p] = ttest2(all_within, m_dists_betw_catgs);
                fprintf("Monkey %s, session %d, int %d, ttest within vs betw p=%0.4g \n", Monkeys(m).Name, sessn, iInt, p)
                disp([mean(all_within), mean(m_dists_betw_catgs)])
                
%                 xlim(xlim_by_monk{m})  % to show var redn. on linear-y plots
                title(sprintf("%s, within vs betw", Monkeys(m).Sessions(sessn).ShortName))
%                 legend(["Within", "Between"])
                xlabel("Mahal. dist.")
                ylabel("Prob.")
                set(gca, "YScale", "log")
                
                formatSVMPlot(gca, gcf)
            end            
        end
    end
end

%% Inspect where high pairwise dists are coming from

large_cutoffs_by_monk = [20 30];

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    disp(Monkeys(m).Name)
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};
        
        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
            
            figure
            hold on

            for i = 1:length(rSessions)
                sessn = rSessions(i);
                
                id = get_good_interval_name2(interval, array, "m_pairwise_dists_within_cats");
                m_dists_within_cats = GEO(m).Sessions(sessn).(id);
                
                id = get_good_interval_name2(interval, array, "m_pairwise_dists_within_dogs");
                m_dists_within_dogs= GEO(m).Sessions(sessn).(id);
                
                id = get_good_interval_name2(interval, array, "m_pairwise_dists_img_pairs_within_cats");
                all_cat_img_pair_inds = GEO(m).Sessions(sessn).(id);
                
                id = get_good_interval_name2(interval, array, "m_pairwise_dists_img_pairs_within_dogs");
                all_dog_img_pair_inds = GEO(m).Sessions(sessn).(id);
                
                subplot(2,1,1)
                hold on
                large_bool = m_dists_within_cats > large_cutoffs_by_monk(m);
                large_inds = all_cat_img_pair_inds(large_bool, :);
                histogram(large_inds(:), 'BinEdges', 1:260)
                title("Cats")
%                 ylim([0 200])
                
                subplot(2,1,2)
                hold on
                large_bool = m_dists_within_dogs > large_cutoffs_by_monk(m);
                large_inds = all_dog_img_pair_inds(large_bool, :);
                
                if m==2
                    outlier_bool = (any((large_inds == 312),2) | any((large_inds == 307),2));
                    large_inds = large_inds(~outlier_bool, :);
                end
                histogram(large_inds(:), 'BinEdges', 261:520)
                title("Dogs")
%                 ylim([0 200])
            end
            legend(["Pre", "Post"])
        end
    end
end
    
%% Functions
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

function kFoldList = makeKFoldList(len, k)
    subSizes = diff(round(linspace(0, len, k+1)));
    regions = repelem(1:k, subSizes);
    kFoldList = regions(randperm(len));
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