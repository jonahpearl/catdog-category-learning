% abcat SVM on subsets of neur. pop. containing strongest catg encoding
% units

% now used for timecourses as well

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
short_name_regexp = '\w*cats\w*';

fullSVMPath = fullfile(EXT_HD, pv_path, 'SVM_results_%g.mat');

% Load behavioral data
load(fullfile(EXT_HD, pv_path, behav_file)) % behavior and neural summaries, but w/o spike times


%% If needed, get date strs and short names of sessions
marta_xls = fullfile(EXT_HD, 'RecordingMarta.xlsx');
max_xls = fullfile(EXT_HD, 'RecordingMax.xlsx');

for m = 1:length(Monkeys)
    date_strs = {Monkeys(m).Sessions.DateStr};
    if regexp(Monkeys(m).Name, 'Marta\w*')
        shortNames = get_short_names(date_strs, marta_xls, short_name_regexp);
    else
        shortNames = get_short_names(date_strs, max_xls, short_name_regexp);
    end
    for i = 1:length(Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).ShortName = shortNames{i};
    end
end

%% Set decoding parameters

%{
The goal here is two fold. 
1) to see how much of the SVM accuracy is coming from inidividual best
neurons.
2) to see how the contributions of those best neurons changed over
training. One might imagine that the top 5 neurons pre/post decode the
same, and it's just an increase in number of neurons post training that
add discriminative power via increasing SNR or via coding along new
dimensions (imgs); OR that the top 5 neurons post decode better than top 5
pre, suggesting individual units may have altered their tuning / new units
have come online with better tuning than observed before. My
"responsiveness" analysis suggests that the latter is the case. Let's see!

%}

% Params specific to this analysis. 

% num_best_units_to_use = 1:1:100;
% num_best_units_to_use = [1:1:50 52:2:200];
% num_best_units_to_use = [50 100 110 120 150 175];
% num_best_units_to_use = [10 25 50 100];
num_best_units_to_use = [100];  % defaulting to 100 for timecourse too now, because MATLAB's SPARSA solver does better at 100 than w more.

ranking_interval = {[175 275]};

% For this analysis, from a reviewer suggestion, we use either the 40 training
% images, or random 40-image subsets of the test set.
training_imgs = {1:20, 261:280};
test_imgs = {21:260, 281:520};

% General SVM Parameters
random_seed = 10; % for reproducibility 
rArrayLocs = {'te'};
% rArrayLocs = {'anterior', 'middle', 'posterior'};
rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)
% rSessionsByMonk = {1:9, 1:7};
% rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};  % Fig 2C
% rSessionsByMonk = {[1 2 3 5 6 7 9], [1 2 6 7]};  % Fig 2C with matched num trials (base 3-5 of Max have super low trial counts)

matchedTrialNumPrePost = true; % if true, match num trials (per catg) used in the classifier pre/post training
fraction_min_num_trials = 0.85;  % what fraction of the min number of trials to use for each classifier (otherwise you get no error bars for the session with the minimum number)

ignoreVal = 20; % if neuron has less than this num spikes, do not use it.
runShuffle = false; % run the shuffled condition?
    nShuffles = 5;
    
n_img_random_subsets = 10;
    
manuallyReduceIntervals = true; % test a subset of all the intervals for faster testing
%     manualIntervals = {[75 175], [175 275]}; % NB, also need to change variable fname_base to grab file containing desired intervals
    manualIntervals = {[175 275]}; % NB, also need to change variable fname_base to grab file containing desired intervals
    
    % for SVM timecourses (Fig 2B)
%     starts = -100:10:300;
%     manualIntervals = arrayfun(@(x) [x x+100], starts, 'UniformOutput', false);
%     ranking_interval = {[175 275]};  % interval to rank units wrt one-D SVM's

% path to spike counts. See pv_get_interval_spike_counts if you want to add
% more intervals.
step = 5; % Interval parameters (for loading the correct spike counts file)
width = 100;
X_fname_base = sprintf('%%s_allNeurons_step%d_wd%d.mat', step, width); % contains full time-course of spike counts (Fig 2b)
% X_fname_base = sprintf('%%s_allNeurons_variableBin_1.mat'); % contains 75-175, 175-225, 175-275, and 175-350.

% Note that folds are not generated purely randomly -- we impose
% constraints about not re-using the same image in train vs test ("abstract
% category") and about balancing the number of cats and dogs in each fold
% (so that the shuffle comes out right at 50 %).
kfold_num = 5; % k-fold cross validation. 

% Other Parameters
TE_LOCS = {'anterior', 'middle', 'posterior'};
catg2_ind1 = 261;

%% Create param struct for saving into record

if ~manuallyReduceIntervals
    ints_used = -1;
else
    ints_used = manualIntervals;
end

paramStruct = struct('RandomSeed', random_seed,...
    'RArrayLocs', {rArrayLocs}, ...
    'IgnoreVal', ignoreVal,...
    'RunShuffle', runShuffle,...
    'KFoldNum', kfold_num, ...
    'SessionsUsed', {rSessionsByMonk},...
    'IntervalsUsed', {ints_used},....
    'TrainingImgSet', 'random_subsets',...
    'TestingImgSet', 'random_subsets',...
    'NumBestUnits', {num_best_units_to_use},...
    'MatchedTrialNumPrePost', matchedTrialNumPrePost);

%% Record of ID nums so far

ID = 180761;  % random 20/20 vs actual 20/20, pre vs post, 175-275

% Load record
fullRecordPath = fullfile(EXT_HD, svmRecordPath);
load(fullRecordPath, 'Record');

%% Load pre-ranked single unit SVM data

ID = 323266;

load(sprintf(fullSVMPath, ID), 'data');
SVM = data;
clear data

%% Run SVMs

% To rank by single single-unit SVMs
ranking_field_name = 'Ranking_SingleUnitSVM_KFL';

if ~matchedTrialNumPrePost
    output_field_name_template = 'KFL_SingleUnitSVMRanking_sparsa_Top_%d';
elseif matchedTrialNumPrePost
    output_field_name_template = 'KFL_SngUntRank_sprs_T_%d_TrNMch_RndImgSub_%d';
end

% svm_solver = {'sgd', 'sparsa'};
svm_solver = {'sparsa'};

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    
    % Calculate num trials to use if matching
    if matchedTrialNumPrePost
        num_trials_per_catg_to_use = inf;
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            Y = SVM(m).Sessions(sessn).Session_Y_catg;
            min_catg_presentations = min([sum(Y==1) sum(Y==2)]);
            if min_catg_presentations < num_trials_per_catg_to_use
                num_trials_per_catg_to_use = min_catg_presentations;
            end
        end
    end
    num_trials_per_catg_to_use = round(num_trials_per_catg_to_use  * ...
                                       fraction_min_num_trials * ...  % actually sub-sample a little bit
                                       (nShuffles-1)/nShuffles * ...  % x fractional size of training set
                                       0.087);  % x 1/10 because we're doing 20/20 here instead of 260/260
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Get data from appropriate storage place
        % X_full is (trials) x (units) x (intervals)
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf(X_fname_base, MonkID);
        [X_full, rIntervals] = load_interval_data(fullfile(EXT_HD, spikeCountPath, fileName)); % X is spike counts, rIntervals is list of each interval
        Y = SVM(m).Sessions(sessn).Session_Y_catg; % list of categories for each trial (image) in X (either 1 or 2)
        imgs_shown = SVM(m).Sessions(sessn).Session_Y_imageID;  % list of imgs for each trial
        
        % For finding intervals in X_full
        rIntervals_original = rIntervals;
        
        for iLoc = 1:length(rArrayLocs)
            array = rArrayLocs{iLoc};
            
            for iInt = 1:length(manualIntervals)
                interval = manualIntervals{iInt};
                
                % Find index in X's 3rd dim for requested interval
                int_idx = find(cellfun(@(a) all(a == interval), rIntervals_original));
                
                % Get unit rankings for this array / interval
%                 field_name = get_good_interval_name2(interval, array, ranking_field_name);
                field_name = get_good_interval_name2(ranking_interval{1}, array, ranking_field_name);
                rankings = [SVM(m).Sessions(sessn).UnitInfo.(field_name)];
                
                for iUnitSubset = 1:length(num_best_units_to_use)
                    unit_subset = 1:(num_best_units_to_use(iUnitSubset));
                    if length(unit_subset) > length(SVM(m).Sessions(sessn).UnitInfo)
                        continue
                    end
                    
                    units_to_use = find(ismember(rankings, unit_subset));
                    
                    % Subset to N best units, and given interval
                    X_subset = X_full(:, units_to_use, int_idx);
                    Y_subset = Y;
                    
                    % Goal is to train/test on different image sets here,
                    % so we will create different "img sets" to be passed
                    % in.
                    for iImgSubset = 1:(n_img_random_subsets+1)
                        
                        % Start by using the 20/20 set
                        if iImgSubset == 1
                            img_set = [training_imgs{1} training_imgs{2}];
                        else
                            img_set = [randsample(test_imgs{1}, 20) randsample(test_imgs{2}, 20)];
                        end
                        
                        % Pre-allocate the storage vector for the SVMs.
                        kflValues = zeros(kfold_num,1);
                        if runShuffle
                            kflValues_SHUFFLE = zeros(kfold_num, nShuffles);
                        end
                        iK = 1;
                        
                        % Create list of which imgs will be tested in each
                        % cv fold
                        num_cues = length(img_set); % how many cues total exist
                        balancedKFoldList = zeros(1, num_cues); % for each cue being used, what fold is it in?
                        balancedKFoldList(1:(num_cues/2)) = makeKFoldList(num_cues/2, kfold_num);  % make folds separately for cats/dogs to assure even catg splits
                        balancedKFoldList((num_cues/2 + 1):end) = makeKFoldList(num_cues/2, kfold_num);
                        
                        % Finally, run the SVM with these folds as the
                        % cross-validation folds.
                        for k = 1:kfold_num 

                            % Which images to use as training in this fold.
                            % IMGSUBSET edit
                            trainingImgs = img_set(balancedKFoldList ~= k); 
                            
                            % Which trials to use as training, based on the
                            % imgs.
                            idx = ismember(imgs_shown, trainingImgs); 
                            
                            % Further match total number of cat/dog images
                            % if requested.
                            if matchedTrialNumPrePost
                                cat_inds = find(ismember(imgs_shown, trainingImgs) & (Y_subset == 1));
                                dog_inds = find(ismember(imgs_shown, trainingImgs) & (Y_subset == 2));
                                cat_inds_to_keep = randsample(cat_inds, num_trials_per_catg_to_use);
                                dog_inds_to_keep = randsample(dog_inds, num_trials_per_catg_to_use);
                                idx(cat_inds(~ismember(cat_inds, cat_inds_to_keep))) = 0;
                                idx(dog_inds(~ismember(dog_inds, dog_inds_to_keep))) = 0;
                            end
                            
                            assert(sum(Y_subset(idx)==1) == sum(Y_subset(idx)==2))
                            disp(sum(Y_subset(idx)==1))

                            % Run the model.
                            model = fitclinear(X_subset(idx,:),...
                                Y_subset(idx),...
                                'Regularization', 'lasso',...
                                'Solver', svm_solver);

                            % Get prediction error.
                            preds = predict(model, X_subset(~idx,:));
                            err = sum(preds ~= Y_subset(~idx)) / sum(~idx);
                            kflValues(iK) = err;
                            
                            if runShuffle
                                for iShuff = 1:nShuffles
                                    Y_subset_SHUFFLE = Y_subset(randperm(length(Y_subset)));
                                    model_SHUFFLE = fitclinear(X_subset(idx,:), Y_subset_SHUFFLE(idx),'Regularization', 'lasso', 'Solver', svm_solver);
                                    preds_SHUFFLE = predict(model_SHUFFLE, X_subset(~idx,:));
                                    err_SHUFFLE = sum(preds_SHUFFLE ~= Y_subset_SHUFFLE(~idx)) / sum(~idx);
                                    kflValues_SHUFFLE(iK, iShuff) = err_SHUFFLE;
                                end
                            end

                            % Iterate.
                            iK = iK+1; 
                        end
                        % Store data. "RndImgSub_0" will be the 20/20 set
                        % and "RndImgSub_1" to "RndImgSub_N" will be the
                        % random subsets.
                        kflID = get_good_interval_name2(interval, array,...
                            sprintf(output_field_name_template, num_best_units_to_use(iUnitSubset), iImgSubset-1));
                        SVM(m).Sessions(sessn).(kflID) = kflValues;
                        disp(mean(kflValues))
                        
                        if runShuffle
                            % needed to write "SHUFF" instead here for some
                            % versions because output exceeded the max field
                            % name length lol.
                            kflID_S = get_good_interval_name2(interval, array, sprintf(strcat(output_field_name_template, '_SHUFF'), num_best_units_to_use(iUnitSubset)));
                            SVM(m).Sessions(sessn).(kflID_S) = kflValues_SHUFFLE;
                        end
                    end
                    % Report progress.
                    fprintf('Done with %s loc %s interval # %d, %d best units, session %d \n', Monkeys(m).Name, array, iInt, num_best_units_to_use(iUnitSubset), sessn)
                end
            end
        end
    end
end

%% Plot a given subset/interval over sessions (Fig 2C)

interval_to_plot = [175, 275];
n_top_units = 100;
base_kfl_name = 'KFL_SngUntRank_sprs_T_%d_TrNMch_RndImgSub_%d';
loc = "te";

mkYLims = {[0.5 0.8], [0.5 0.65]};  % fig 2c

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    subplot(2,1,m)
    hold on

    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Pull the data back out
        random_subset_kfls = [];
        for iImgSubset = 1:(n_img_random_subsets+1)
            kfl_name = sprintf(base_kfl_name, n_top_units, iImgSubset-1);
            kflID = get_good_interval_name2(interval_to_plot, loc, kfl_name);
            if iImgSubset == 1
                training_kfl = 1 - SVM(m).Sessions(sessn).(kflID);
            else
                random_subset_kfls = [random_subset_kfls; 1 - SVM(m).Sessions(sessn).(kflID)];
            end
        end
        
        % Plot real 20/20 data
        errorbar(i, mean(training_kfl), std(training_kfl)/sqrt(length(training_kfl)), ...
            'o', 'LineWidth', 2, 'Color', mlc(i), ...
            'MarkerSize', 0.1)
        
        % Add random subsets alongside
        errorbar(i + 0.25, mean(random_subset_kfls), std(random_subset_kfls)/sqrt(length(random_subset_kfls)), ...
            'o', 'LineWidth', 2, 'Color', mlc(i), 'HandleVisibility', "off", ...
            'MarkerSize', 0.1)
        
        disp([mean(training_kfl) mean(random_subset_kfls)])
    end
    
    formatSVMPlot(gca, gcf, 16)
    title(sprintf('%s', Monkeys(m).Name), 'Interpreter', 'none')
    xticks([1 1.25 2 2.25])
    xticklabels(["20/20" "Rnd" "20/20" "Rnd"])
    legend(["Pre", "Post"])

    % Save the plots
%     pause(0.5)
%     saveas(gcf, fullfile(figureSavePath, sprintf('pv_SVM_Timecourse_%s_%s_%g', Monkeys(m).Name, sigID, ID)), 'epsc')
end

%% =============== %%

%% Run SVMs, subtracting more good units each time

% To rank by single single-unit SVMs
ranking_field_name = 'Ranking_SingleUnitSVM_KFL';
output_field_name_template = 'SingleUnitRemovalSVMs_KFL_Without_Top_%d';

% svm_solver = {'sgd', 'sparsa'};
svm_solver = {'sparsa'};

num_best_units_to_remove = [1:50 52:2:200];

% To rank by GLM coeff
% ranking_field_name = 'Ranking_GLM_coeff';
% output_field_name_template = 'KFL_GLMRanking_Top_%d';

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Get data from appropriate storage place
        % X_full is (trials) x (units) x (intervals)
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf(X_fname_base, MonkID);
        [X_full, rIntervals] = load_interval_data(fullfile(EXT_HD, spikeCountPath, fileName)); % X is spike counts, rIntervals is list of each interval
        Y = SVM(m).Sessions(sessn).Session_Y_catg; % list of categories for each trial (image) in X (either 1 or 2)
        imgs_shown = SVM(m).Sessions(sessn).Session_Y_imageID;  % list of imgs for each trial
        
        % For finding intervals in X_full
        rIntervals_original = rIntervals;
        
        for iLoc = 1:length(rArrayLocs)
            array = rArrayLocs{iLoc};
            
            for iInt = 1:length(manualIntervals)
                interval = manualIntervals{iInt};
                
                % Find index in X's 3rd dim for requested interval
                int_idx = find(cellfun(@(a) all(a == interval), rIntervals_original));
                
                % Get unit rankings for this array / interval
                field_name = get_good_interval_name2(interval, array, ranking_field_name);
                rankings = [SVM(m).Sessions(sessn).UnitInfo.(field_name)];
                
                for iUnitSubset = 1:length(num_best_units_to_remove)
                    
                    % get output field name, skip if already present
                    kflID = get_good_interval_name2(interval, array, sprintf(output_field_name_template, num_best_units_to_remove(iUnitSubset)));
                    if ismember(kflID, fields(SVM(m).Sessions(sessn))) && ~isempty(SVM(m).Sessions(sessn).(kflID))
                        continue
                    end
                    
                    % start from 100 best units, because > 100 units and
                    % the SVM solver starts to do weird things.
                    initial_units_to_use = find(ismember(rankings, 1:100));
                    
                    % remove the n best units
                    rankings_to_remove = 1:(num_best_units_to_remove(iUnitSubset));
                    units_to_remove = find(ismember(rankings, rankings_to_remove));
                    units_to_use = initial_units_to_use(~ismember(initial_units_to_use, units_to_remove));
                    
                    % Pre-allocate the storage vector for the SVMs.
                    kflValues = zeros(kfold_num,1);
                    iK = 1;
                    
                    % Subset to N best units, and given interval
                    X_subset = X_full(:, units_to_use, int_idx);
                    Y_subset = Y;
                    
                    % IMGSUBSET edit
                    % if we're training and testing on the same img set,
                    % run cross validation as normal.
                    % If we're training and testing on different img sets,
                    % it's much easier, just train and test on the full
                    % data matrices since there's already no overlap
                    % between them.
                    
                    if svm_train_set_idx == svm_test_set_idx
                        % Create balanced image folds, for "abstract category"
                        % decoding.
                        num_cues = length(unique(img_sets{svm_train_set_idx})); % how many cues total exist
                        balancedKFoldList = zeros(1, num_cues); % for each cue being used, what fold is it in?
                        balancedKFoldList(1:(num_cues/2)) = makeKFoldList(num_cues/2, kfold_num);  % make folds separately for cats/dogs to assure even catg splits
                        balancedKFoldList((num_cues/2 + 1):end) = makeKFoldList(num_cues/2, kfold_num);
                        
                        % Finally, run the SVM with these folds as the
                        % cross-validation folds.
                        for k = 1:kfold_num 

                            % Which images to use as training in this fold.
                            % IMGSUBSET edit
                            trainingImgs = img_sets{svm_train_set_idx}(balancedKFoldList ~= k); 

                            
                            % Which trials to use as training, based on the
                            % imgs.
                            idx = ismember(imgs_shown, trainingImgs); 

                            % Run the model.
                            model = fitclinear(X_subset(idx,:),...
                                Y_subset(idx),...
                                'Regularization', 'lasso',...
                                'Solver', svm_solver);

                            % super slow and doesn't seem to help
%                             model = fitcsvm(X_subset(idx,:),...
%                             Y_subset(idx));


                            % Get prediction error.
                            preds = predict(model, X_subset(~idx,:));
                            err = sum(preds ~= Y_subset(~idx)) / sum(~idx);
                            kflValues(iK) = err;

                            % Iterate.
                            iK = iK+1;
                        end
                        
                    elseif svm_train_set_idx ~= svm_test_set_idx
                        % IMGSUBSETS edit
                        % Which trials to use as training and testing, based
                        % on the design of this analysis.
                        svm_train_bool = ismember(imgID_Subset, img_sets{svm_train_set_idx}); % on each trial, should we use this data for the SVM?
                        svm_test_bool = ismember(imgID_Subset, img_sets{svm_test_set_idx});
                        
                        % Here, our training and test data matrices are
                        % different by design, so we can't really
                        % cross-validate per se. We'll just run the full
                        % model kfold_num times in order to generate some
                        % error bars wrt the random seed.
                        for k = 1:kfold_num 
                            % Run the model.
                            model = fitclinear(X_subset(svm_train_bool,:),...
                                Y_subset(svm_train_bool),...
                                'Regularization', 'lasso',...
                                'Solver', svm_solver);
                            % Get prediction error.
                            preds = predict(model, X_subset(svm_test_bool,:));
                            err = sum(preds ~= Y_subset(svm_test_bool)) / sum(svm_test_bool);
                            kflValues(iK) = err;
                            % Iterate.
                            iK = iK+1;
                        end
                    end
                    
                    % Store data.
                    SVM(m).Sessions(sessn).(kflID) = kflValues;
    %                 if runShuffle
    %                     kflID_S = get_good_interval_name2(interval, loc, 'KFL_SHUFFLE');
    %                     SVM(m).Sessions(sessn).(kflID_S) = kflValues_SHUFFLE;
    %                 end
    
                    % Report progress.
                    fprintf('Done with %s loc %s interval # %d, removed %d best units (%d remaining), session %d \n', SVM(m).Name, array, iInt, num_best_units_to_remove(iUnitSubset), length(units_to_use), sessn)
                end
            end
        end
    end
end

%% Plot performance over sessions of removing units

num_best_units_to_remove = [1:50 52:2:200];

array = 'te';
% base_kfl_name = 'KFL_GLMRanking_Top_%d';
base_kfl_name = 'SingleUnitRemovalSVMs_KFL_Without_Top_%d';
fig_path = '/Users/jonahpearl/Documents/BJR group/Catdog paper/Feb 2023 addtl figs';

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};

            % pre alloc vector to hold mean + std for accuracies
            means = nan(length(rSessions), length(num_best_units_to_remove));
            stds = nan(length(rSessions), length(num_best_units_to_remove));

            for i = 1:length(rSessions)
                sessn = rSessions(i);
                for iUnitSubset = 1:length(num_best_units_to_remove)
                    kflID = get_good_interval_name2(interval, array, sprintf(base_kfl_name, num_best_units_to_remove(iUnitSubset)));
                    if ~ismember(kflID, fields(SVM(m).Sessions(sessn)))
                        continue
                    end
                    kfls = SVM(m).Sessions(sessn).(kflID); 

                    means(i, iUnitSubset) = mean(1 - kfls);
                    stds (i, iUnitSubset) = std(1 - kfls);
                end
            end

            % plot the data
            figure
            hold on
            for i = 1:length(rSessions)
                errorbar(num_best_units_to_remove, means(i,:), stds(i,:), 'DisplayName', Monkeys(m).Sessions(rSessions(i)).ShortName)
            end
            set(gca, 'XScale', 'log')
            xlabel('Num best units removed')
            ylabel('Accuracy')
            title(sprintf('%s, interval %d - %d', Monkeys(m).Name, interval(1), interval(2)), 'Interpreter', 'none')
            legend
%             saveas(gcf, fullfile(fig_path, sprintf('%s_%d_%d.pdf', Monkeys(m).Name, interval(1), interval(2))))
        end
    end
end

%% Curve fitting for ^^
% looks like the data are bi-sigmoidal lol. We'll see if i can verify this
% rigorously but first let's figure out curve fitting...

array = 'te';
% base_kfl_name = 'KFL_GLMRanking_Top_%d';
% base_kfl_name  = 'KFL_SingleUnitSVMRanking_Top_%d';
base_kfl_name = 'SingleUnitRemovalSVMs_KFL_Without_Top_%d';
fig_path = '/Users/jonahpearl/Documents/BJR group/Catdog paper/Feb 2023 addtl figs';

% ok this is totally contrived, let's just do a single sigmoid
% bisigmoid = fittype('a1/(1 + exp(-b1*(x-c1))) + a2/(1 + exp(-b2*(x-c2))) + d',...
%     'dependent', 'y', 'independent', 'x',...
%     'coefficients', {'a1', 'b1', 'c1', 'a2', 'b2', 'c2', 'd'});
num_best_units_to_remove = [1:50 52:2:96];

sigmoid = fittype('a1/(1 + exp(-b1*(x-c1))) + d',...
    'dependent', 'y', 'independent', 'x',...
    'coefficients', {'a1', 'b1', 'c1', 'd'});


for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};

            % pre alloc vector to hold mean + std for accuracies
            means = nan(length(rSessions), length(num_best_units_to_remove));
            stds = nan(length(rSessions), length(num_best_units_to_remove));

            for i = 1:length(rSessions)
                sessn = rSessions(i);
                for iUnitSubset = 1:length(num_best_units_to_remove)
                    kflID = get_good_interval_name2(interval, array, sprintf(base_kfl_name, num_best_units_to_remove(iUnitSubset)));
                    if ~ismember(kflID, fields(SVM(m).Sessions(sessn)))
                        break
                    end
                    kfls = SVM(m).Sessions(sessn).(kflID); 

                    means(i, iUnitSubset) = mean(1 - kfls);
                    stds (i, iUnitSubset) = std(1 - kfls);
                end
            end

            figure
            hold on
            xvals = log10(num_best_units_to_remove);
            yl = [0.48 max(means(:)+stds(:))+0.02];
            for i = 1:length(rSessions)
                f = fit(xvals', means(i,:)', sigmoid, 'start', [.05, 5, 1, 0.5]);
                errorbar(xvals, means(i,:), stds(i,:),...
                    'DisplayName', Monkeys(m).Sessions(rSessions(i)).ShortName,...
                    'Color', mlc(i))
                plot(xvals, f(xvals), '--', 'Color', mlc(i) + 0.1, 'HandleVisibility', 'off')
                plot([f.c1 f.c1], yl, '-', 'Color', mlc(i) + 0.1, 'HandleVisibility', 'off', 'LineWidth', 1.25)
                midpoint_slope = 0.25 * f.a1 * f.b1;  % take dy/dx of A/(1+exp(-b(x-c))), eval at x=c --> dy/dx = 1/4*A*b.
                fprintf('Monkey %s, session %s: amp %g, midpoint %g, b %g, midpoint slope %g, offset %g \n',...
                    Monkeys(m).Name,...
                    Monkeys(m).Sessions(rSessions(i)).ShortName,...
                        round(f.a1,3),...
                        round(f.c1,3),...
                        round(f.b1,3),...
                        round(midpoint_slope,3),...
                        round(f.d, 3))
            end
            set(gca, 'XScale', 'log')
            xlabel('Log10 num best units used')
            ylabel('Accuracy')
            title(sprintf('%s, interval %d - %d', Monkeys(m).Name, interval(1), interval(2)), 'Interpreter', 'none')
            legend
            
        end
    end
end

%% =============== %%

%% Run SVMs across arrays (supp fig 3)

% To rank by single single-unit SVMs
ranking_field_name = 'Ranking_SingleUnitSVM_KFL';
output_field_name_template = 'KFL_SingleUnitSVMRanking_sparsa_Top_%d';

% svm_solver = {'sgd', 'sparsa'};
svm_solver = {'sparsa'};

rArrayLocs = {'anterior', 'middle', 'posterior'};
rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)

runShuffle = false; % run the shuffled condition? obv will slow things down ~(nShuffles)-fold!
nShuffles = 5;

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Get data from appropriate storage place
        % X_full is (trials) x (units) x (intervals)
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf(X_fname_base, MonkID);
        [X_full, rIntervals] = load_interval_data(fullfile(EXT_HD, spikeCountPath, fileName)); % X is spike counts, rIntervals is list of each interval
        Y = SVM(m).Sessions(sessn).Session_Y_catg; % list of categories for each trial (image) in X (either 1 or 2)
        imgs_shown = SVM(m).Sessions(sessn).Session_Y_imageID;  % list of imgs for each trial
        
        % For finding intervals in X_full
        rIntervals_original = rIntervals;
        
        for iLoc = 1:length(rArrayLocs)
            array = rArrayLocs{iLoc};
            
            for iInt = 1:length(manualIntervals)
                interval = manualIntervals{iInt};
                
                % Find index in X's 3rd dim for requested interval
                int_idx = find(cellfun(@(a) all(a == interval), rIntervals_original));
                
                % Get unit rankings for this array / interval
%                 field_name = get_good_interval_name2(interval, array, ranking_field_name);
                field_name = get_good_interval_name2(ranking_interval{1}, array, ranking_field_name);
                rankings = [SVM(m).Sessions(sessn).UnitInfo.(field_name)];
                
                for iUnitSubset = 1:length(num_best_units_to_use)
                    unit_subset = 1:(num_best_units_to_use(iUnitSubset));
                    if length(unit_subset) > length(SVM(m).Sessions(sessn).UnitInfo)
                        continue
                    end
                    
                    units_to_use = find(ismember(rankings, unit_subset));
                    
                    % Pre-allocate the storage vector for the SVMs.
                    kflValues = zeros(kfold_num,1);
                    if runShuffle
                        kflValues_SHUFFLE = zeros(kfold_num, nShuffles);
                    end
                    iK = 1;
                    
                    % Subset to N best units, and given interval
                    X_subset = X_full(:, units_to_use, int_idx);
                    Y_subset = Y;
                    
                    % IMGSUBSET edit
                    % if we're training and testing on the same img set,
                    % run cross validation as normal.
                    % If we're training and testing on different img sets,
                    % it's much easier, just train and test on the full
                    % data matrices since there's already no overlap
                    % between them.
                    
                    if svm_train_set_idx == svm_test_set_idx
                        % Create balanced image folds, for "abstract category"
                        % decoding.
                        num_cues = length(unique(img_sets{svm_train_set_idx})); % how many cues total exist
                        balancedKFoldList = zeros(1, num_cues); % for each cue being used, what fold is it in?
                        balancedKFoldList(1:(num_cues/2)) = makeKFoldList(num_cues/2, kfold_num);  % make folds separately for cats/dogs to assure even catg splits
                        balancedKFoldList((num_cues/2 + 1):end) = makeKFoldList(num_cues/2, kfold_num);
                        
                        % Finally, run the SVM with these folds as the
                        % cross-validation folds.
                        for k = 1:kfold_num 

                            % Which images to use as training in this fold.
                            % IMGSUBSET edit
                            trainingImgs = img_sets{svm_train_set_idx}(balancedKFoldList ~= k); 

                            
                            % Which trials to use as training, based on the
                            % imgs.
                            idx = ismember(imgs_shown, trainingImgs); 

                            % Run the model.
                            model = fitclinear(X_subset(idx,:),...
                                Y_subset(idx),...
                                'Regularization', 'lasso',...
                                'Solver', svm_solver);

                            % super slow and doesn't seem to help
%                             model = fitcsvm(X_subset(idx,:),...
%                             Y_subset(idx));


                            % Get prediction error.
                            preds = predict(model, X_subset(~idx,:));
                            err = sum(preds ~= Y_subset(~idx)) / sum(~idx);
                            kflValues(iK) = err;
                            
                            if runShuffle
                                for iShuff = 1:nShuffles
                                    Y_subset_SHUFFLE = Y_subset(randperm(length(Y_subset)));
                                    model_SHUFFLE = fitclinear(X_subset(idx,:), Y_subset_SHUFFLE(idx),'Regularization', 'lasso', 'Solver', svm_solver);
                                    preds_SHUFFLE = predict(model_SHUFFLE, X_subset(~idx,:));
                                    err_SHUFFLE = sum(preds_SHUFFLE ~= Y_subset_SHUFFLE(~idx)) / sum(~idx);
                                    kflValues_SHUFFLE(iK, iShuff) = err_SHUFFLE;
                                end
                            end

                            % Iterate.
                            iK = iK+1;
                        end
                        
                    elseif svm_train_set_idx ~= svm_test_set_idx
                        % IMGSUBSETS edit
                        % Which trials to use as training and testing, based
                        % on the design of this analysis.
                        svm_train_bool = ismember(imgID_Subset, img_sets{svm_train_set_idx}); % on each trial, should we use this data for the SVM?
                        svm_test_bool = ismember(imgID_Subset, img_sets{svm_test_set_idx});
                        
                        % Here, our training and test data matrices are
                        % different by design, so we can't really
                        % cross-validate per se. We'll just run the full
                        % model kfold_num times in order to generate some
                        % error bars wrt the random seed.
                        for k = 1:kfold_num 
                            % Run the model.
                            model = fitclinear(X_subset(svm_train_bool,:),...
                                Y_subset(svm_train_bool),...
                                'Regularization', 'lasso',...
                                'Solver', svm_solver);
                            % Get prediction error.
                            preds = predict(model, X_subset(svm_test_bool,:));
                            err = sum(preds ~= Y_subset(svm_test_bool)) / sum(svm_test_bool);
                            kflValues(iK) = err;
                            % Iterate.
                            iK = iK+1;
                        end
                    end
                    
                    % Store data.
                    kflID = get_good_interval_name2(interval, array, sprintf(output_field_name_template, num_best_units_to_use(iUnitSubset)));
                    SVM(m).Sessions(sessn).(kflID) = kflValues;
                    if runShuffle
                        kflID_S = get_good_interval_name2(interval, loc, sprintf(strcat(output_field_name_template, '_SHUFFLE'), num_best_units_to_use(iUnitSubset)));
                        SVM(m).Sessions(sessn).(kflID_S) = kflValues_SHUFFLE;
                    end
    
                    % Report progress.
                    fprintf('Done with %s loc %s interval # %d, %d best units, session %d \n', Monkeys(m).Name, array, iInt, num_best_units_to_use(iUnitSubset), sessn)
                end
            end
        end
    end
end

%% Plot timecourse across arrays 

% Plotting params
plot_alpha = 0.4; % transparency of sem fill
mkYLims = {[0.45 0.85], [0.45 0.65]};

rArrayLocs = {'anterior', 'middle', 'posterior'};
top_n_units = 100;  % top 100 for fig 2B
base_kfl_name = sprintf('KFL_SingleUnitSVMRanking_sparsa_Top_%d', top_n_units);  % top 10, 25, 50, or 100

% f = figure2('Position', [400 400 1000 600])
figure('Position', [400 400 500 300])
% tiledlayout(length(Monkeys), length(rArrayLocs))
for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    % Pre-allocate vectors for plotting
    accMeans = zeros(length(rSessions), length(manualIntervals), length(rArrayLocs));
    accSems = zeros(length(rSessions), length(manualIntervals), length(rArrayLocs));
    
    % Prepare figure
%     figure2('Position', [400 400 1000 600])
%     subplot(2,1,m)
    
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        % New subplot for each array location
%         subplot(2, ceil(length(rArrayLocs)/2), iLoc)
        nexttile
        hold on

        % Big TE with smaller arrays underneath
%         if strcmp(loc, 'te')
%             subplot(2, 3,[1 2 3])
%             hold on
%         else
%             subplot(2,3,iLoc+2)
%             hold on
%         end
        
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            for iInt = 1:length(manualIntervals)
                interval = manualIntervals{iInt};
                
                % Get field names
                kflID = get_good_interval_name2(interval, loc, base_kfl_name);
                
                % Get data
                kfls = SVM(m).Sessions(sessn).(kflID);
                kfls = kfls(:); % reshape into 1 x numel
                accMeans(i, iInt, iLoc) = mean(1 - kfls);
                accSems(i, iInt, iLoc) = std(1 - kfls) / numel(kfls);
            end
            
            
            % Simplify variable names
            meanVector = accMeans(i, :, iLoc);
            semVector = accSems(i, :, iLoc);
            
            % Plot a line for each session with filled sem.
            % mlc() is a function to get the RGB vals for the default
            % MATLAB colors.
            plot(starts, meanVector, '-', 'LineWidth', 2, 'Color', mlc(i),...
                'DisplayName', Monkeys(m).Sessions(sessn).ShortName)
            fill([starts fliplr(starts)],...
                [meanVector + semVector, fliplr(meanVector - semVector)], mlc(i),...
                'FaceAlpha', plot_alpha,'linestyle','none', ...
                'HandleVisibility', 'off');
            
            % Plot signf inds
            sigID = sprintf('ClustPermSignfInds_%s_Sessions_%d_vs_%d_Top%d', loc, rSessions(1), rSessions(2), top_n_units);
            try 
                if ~ isempty(Monkeys(m).(sigID))
                    plot(starts(Monkeys(m).(sigID)), mkYLims{m}(2)-0.02, 'ko', 'MarkerFaceColor', 'k')
                end
            catch
                warning('No field found for cluster permutation statistics')
            end
            
            % Add labels
            if iLoc == 1
                xlabel('Time from cue on (ms)')
                ylabel('SVM accuracy')
%                 legend('Location', 'eastoutside')
            end
            
            % Make graphs have the same y-axes within each monkey
            ylim(mkYLims{m})
            
            % Detailed labels if desired
            %             ylabel('SVM abcat accuracy (mean +/- SEM)')
            title(sprintf('%s', loc), 'Interpreter', 'none')
            
            
            % Make the plot look nice
            formatSVMPlot(gca, gcf, 20)
        end
    end
    % More detailed labels if desired.
    %         sgtitle(sprintf('%s', Monkeys(m).Name), 'Interpreter', 'none')

    % Save the plots
%     pause(0.5)
%     saveas(gcf, fullfile(figureSavePath, sprintf('pv_SVM_Timecourse_%s_%s_%g', Monkeys(m).Name, sigID, ID)), 'epsc')
end

%% =============== %%
%% Save the data

% Remove unecessary fields
for m = 1:length(SVM)
%     SVM(m).Sessions = rmfield(SVM(m).Sessions, {'CueInfo', 'TrialInfo', 'TimesBT', 'CodesBT', 'Fixn_err_TC'});
    
    % (storing these is redundant with param struct, but makes me feel
    % safer)
    
    % if we use some metric other than GLM coeffs to rank best units,
    % should still be ok as long as the field name differs.
    
    % Store intervals used
    SVM(m).RIntervals = rIntervals;
    
    % Store num best units used
    SVM(m).Num_best_units = num_best_units_to_use;
     
end


% Save the data and add a row to the Record
fullSVMPath = fullfile(EXT_HD, pv_path, 'SVM_results_%g.mat');
fullRecordPath = fullfile(EXT_HD, svmRecordPath);
save_SVM_data(SVM, paramStruct, fullSVMPath, fullRecordPath);
fprintf('Data saved. \n')

%% Functions

function tStats = getrealtstats(Data, m, rSessions, loc, rIntervals, kfl_base_name)
    tStats = zeros(1, length(rIntervals));
    for iInt = 1:length(rIntervals)
        interval = rIntervals{iInt};

        % Get field names
        kflID = get_good_interval_name2(interval, loc, kfl_base_name);

        % Run t-tests
        s1 = Data(m).Sessions(rSessions(1)).(kflID);
        s2 = Data(m).Sessions(rSessions(2)).(kflID);
        [~,~,~,stats] = ttest2(s1, s2, 'Tail', 'both');
        tStats(iInt) = stats.tstat;
    end
end

function tStats_SHUFFLED = getshuffledtstats(Data, m, rSessions, loc, rIntervals, kfl_base_name)

    % First calculate number of shuffles
    kflID = get_good_interval_name2(rIntervals{1}, loc, kfl_base_name);
    exampleKFLs = Data(m).Sessions(rSessions(1)).(kflID);
    nShuffles = size(exampleKFLs,2);  % SVM shuffles have been pre-run; now just extract the results
    
    tStats_SHUFFLED = zeros(nShuffles, length(rIntervals));
    for iInt = 1:length(rIntervals)
        interval = rIntervals{iInt};
        
        % Get field names
        kflID = get_good_interval_name2(interval, loc, kfl_base_name);

        % Get data
        s1 = Data(m).Sessions(rSessions(1)).(kflID);
        s2 = Data(m).Sessions(rSessions(2)).(kflID);
        
        for iShuff = 1:nShuffles
            [~,~,~,stats] = ttest2(s1(:,iShuff), s2(:, iShuff), 'Tail', 'both');
            tStats_SHUFFLED(iShuff, iInt) = stats.tstat;
        end
    end
end

function clusterStats = calcShuffClusterStats(tStats, nClust)
    
    % Run the same algorithm as the real data
    clusterStats = calculateTClusterStats(tStats);
    
    % make same length as real data
    if length(clusterStats) < nClust
        clusterStats(length(clusterStats):nClust) = 0;
    elseif length(clusterStats) > nClust
        clusterStats(nClust+1:end) = []; % sort is descending, so this removes smallest values
    end
end

function [clusterStats, intsInCluster] = calculateTClusterStats(tStats)
    starts = [1 find(diff(sign(tStats)))+1];
    ends = [starts(2:end)-1 length(tStats)];
    clusters = cell(1, length(starts));
    intsInCluster = cell(1, length(starts)); % inds in rIntervals that belong to each cluster, for plotting later
    for iClust = 1:length(clusters)
        clusters{iClust} = tStats(starts(iClust):ends(iClust));
        intsInCluster{iClust} = starts(iClust):ends(iClust);
    end
    [clusterStats, idx] = sort(abs(cellfun(@sum, clusters)), 'descend');
    intsInCluster = intsInCluster(idx); % sort these to correspond to sorted cluster list
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