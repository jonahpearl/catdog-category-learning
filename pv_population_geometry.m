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
                    MD_within_cats = [MD_within_cats; mahal(X_test(Y_test == 1, :), X_train(Y_train == 1, :))];
                    
                    % Within dogs
                    MD_within_dogs =[MD_within_dogs; mahal(X_test(Y_test == 2, :), X_train(Y_train == 2, :))];
                    
                    % Cats wrt dogs
                    MD_cats_wrt_dogs = [MD_cats_wrt_dogs; mahal(X_test(Y_test == 1, :), X_train(Y_train == 2, :))];
                    
                    % Dogs wrt cats
                    MD_dogs_wrt_cats = [MD_dogs_wrt_cats; mahal(X_test(Y_test == 2, :), X_train(Y_train == 1, :))];
                    
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

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    disp(Monkeys(m).Name)
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};
        
        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
            
            figure('Position', [500, 700, 550, 250])
            hold on

            for i = 1:length(rSessions)
                sessn = rSessions(i);
                
                within_cats_id = get_good_interval_name2(interval, array, "MD_within_cats");
                MD_within_cats = GEO(m).Sessions(sessn).(within_cats_id);
                MD_within_cats = MD_within_cats(~isnan(MD_within_cats));
                
                within_dogs_id = get_good_interval_name2(interval, array, "MD_within_dogs");
                MD_within_dogs= GEO(m).Sessions(sessn).(within_dogs_id);
                MD_within_dogs = MD_within_dogs(~isnan(MD_within_dogs));
                
                id = get_good_interval_name2(interval, array, "MD_cats_wrt_dogs");
                MD_cats_wrt_dogs= GEO(m).Sessions(sessn).(id);
                MD_cats_wrt_dogs = MD_cats_wrt_dogs(~isnan(MD_cats_wrt_dogs));
                
                id = get_good_interval_name2(interval, array, "MD_dogs_wrt_cats");
                MD_dogs_wrt_cats= GEO(m).Sessions(sessn).(id);
                MD_dogs_wrt_cats = MD_dogs_wrt_cats(~isnan(MD_dogs_wrt_cats));
                
%                 disp(size(MD_within_dogs))
                
                subplot(1,2,1)
                hold on
                errorbar(i, mean(MD_within_cats), std(MD_within_cats)/sqrt(length(MD_within_cats)),...
                    "o", 'LineWidth', 2, 'Color', 'r') 
                errorbar(i+0.25, mean(MD_within_dogs), std(MD_within_dogs)/sqrt(length(MD_within_dogs)),...
                    "o", 'LineWidth', 2, 'Color', 'b')
                
                subplot(1,2,2)
                hold on
                errorbar(i, mean(MD_cats_wrt_dogs), std(MD_cats_wrt_dogs)/sqrt(length(MD_cats_wrt_dogs)),...
                    "o", 'LineWidth', 2, 'Color', 'r') 
                errorbar(i+0.25, mean(MD_dogs_wrt_cats), std(MD_dogs_wrt_cats)/sqrt(length(MD_dogs_wrt_cats)),...
                    "o", 'LineWidth', 2, 'Color', 'b') 
                
            end
            
            ax1 = subplot(1,2,1);
%             legend(["Cat", "Dog"])
            title("Within categories")
            xlim([0.5, 2.5])
%             formatSVMPlot(gca, gcf);
            
            ax2 = subplot(1,2,2);
%             legend(["Catgs wrt dogs", "Dogs wrt cats"]) 
            title("Across categories") 
            xlim([0.5, 2.5])
%             formatSVMPlot(gca, gcf);
            
            sgtitle(sprintf("Monk %s, interval %d to %d", Monkeys(m).Name, interval(1), interval(2)), 'Interpreter', 'none')
            xlabel("Session")
            ylabel("Mahalanobis dist.")
            
%             linkaxes([ax1, ax2])
%             if m == 1
%                 ylim([125, 155])
%             elseif m == 2
%                 ylim([170, 180])  
%             end
            
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
             
        end
        
    end
end


%% Pairewise dist metrics take 2

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
                fprintf("Monkey %s, session %d, int %d \n", Monkeys(m).Name, i, iInt)

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

%                             disp(cond(X_train(Y_train==1, :)))
                            condition_num = cond(X_train(Y_train==1, :));
                            if condition_num > 300
                                n_attempts = n_attempts + 1;
%                                 fprintf("%d \n ", n_attempts)
                                continue
                            end
                            
                            % Get the pairwise dists
                            cat_dists = squareform(pdist(X_train(Y_train==1, :), 'mahalanobis'));
                            dog_dists = squareform(pdist(X_train(Y_train==2, :), 'mahalanobis'));
                            
                            done = true;
                        catch
                            n_attempts = n_attempts + 1;
%                             fprintf("%d \n ", n_attempts)
                        end
                    end

                    % Ensure the only zeros are on the diagonal
                    assert(sum(sum(cat_dists == 0)) == size(cat_dists,1))
                    assert(sum(sum(dog_dists == 0)) == size(dog_dists,1))

                    % Filter out by stimulus
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

                end
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


%% Plot pairwise dist results

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    disp(Monkeys(m).Name)
    
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};
        
        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
            
            figure
            hold on

            variances = zeros(2, 2);  % sessions x categories
            for i = 1:length(rSessions)
                sessn = rSessions(i);
                
                id = get_good_interval_name2(interval, array, "m_pairwise_dists_within_cats");
                m_dists_within_cats = GEO(m).Sessions(sessn).(id);
                variances(i, 1) = var(m_dists_within_cats);
                
                id = get_good_interval_name2(interval, array, "m_pairwise_dists_within_dogs");
                m_dists_within_dogs = GEO(m).Sessions(sessn).(id);
                variances(i, 2) = var(m_dists_within_dogs);
                
                % Try removing imgs 307/312 which are outliers in mk 2
                if m == 2
                    id = get_good_interval_name2(interval, array, "m_pairwise_dists_img_pairs_within_dogs");
                    all_dog_img_pair_inds = GEO(m).Sessions(sessn).(id);
                    outlier_bool = (any((all_dog_img_pair_inds == 312),2) | any((all_dog_img_pair_inds == 307),2));
                    m_dists_within_dogs = m_dists_within_dogs(~outlier_bool);
                end
                
%                 subplot(1,2,1)
%                 hold on
%                 errorbar(i, mean(m_dists_within_cats(:)), std(m_dists_within_cats(:))/sqrt(length(m_dists_within_cats)),...
%                     "o", 'LineWidth', 2, 'Color', mlc(1)) 
%                 errorbar(i+0.25, mean(m_dists_within_dogs(:)), std(m_dists_within_dogs(:))/sqrt(length(m_dists_within_dogs)),...
%                     "o", 'LineWidth', 2, 'Color', mlc(2))
                subplot(2,1,1)
                hold on
                histogram(m_dists_within_cats, 'BinEdges', 0:0.5:50)
                title("Cats")
                set(gca, 'YScale', 'log')
                
                subplot(2,1,2)
                hold on
                histogram(m_dists_within_dogs, 'BinEdges', 0:0.5:50)
                title("Dogs")
                set(gca, 'YScale', 'log')
                
            end
            
            subplot(2,1,1)
            legend(["Pre", "Post"])
            sgtitle(sprintf("Monkey %s", Monkeys(m).Name), "Interpreter", "none")
            
            disp(variances)
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