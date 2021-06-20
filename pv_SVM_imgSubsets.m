
% SVMs on different image subsets

%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
svmRecordPath = 'XMA2/Monkey_structs/SVM_Records.mat';
pv_path = 'XMA2/Monkey_structs';

% Load behavioral data
load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_behav_and_metaNI.mat')) % behavior and neural summaries, but w/o spike times

%% Set parameters

%{
The goal here is to see how much of the SVM accuracy is coming from the
20/20 training image set, which the monkeys seem to have memorized. Since
we're already using abcat SVMs, it's unlikel that the entire accuracy is
coming from the 20/20 set, but it's probably contributing more than its
proportion.

% Key changes of the SVM code for this analysis are marked with IMGSUBSET in
comments.

To test this:
-- train and test the SVMs only on the 240/240 img set
-- train and test only on 20/20? Though this will probably not go well with
abcat since the trial nums are tiny.
-- train the SVMs on the 20/20 set, and test on the 240/240 set
-- vice versa


%}

% Params specific to this img subset analysis
% "train" --> 20/20 set. "test" --> 240/240 set.
img_sets = {[1:20 261:280], [21:260 281:520]}; 
img_set_names = {'train', 'test'};
svm_train_set_idx = 1; % 1 for 20/20 training set, 2 for 240/240 testing set
svm_test_set_idx = 1;

% SVM Parameters
random_seed = 10; % for reproducibility 
rArrayLocs = {'te', 'anterior', 'middle', 'posterior'}; % relevant subsets
% rArrayLocs = {'te'};
% rSessionsByMonk = {[7 9], [6 7]};
rSessionsByMonk = {1:9, 1:7};
ignoreVal = 20; % if neuron has less than this num spikes, do not use it.
runShuffle = false; % run the shuffled condition?
    nShuffles = 5;
matchInputSizes = false; % make input matrices all same size (control) ?
    % These parameters only used if matching input sizes. If false, these
    % parameters are automatically set to NaN, NaN, and 1.
    proportion_of_min_pop_size = 0.75; 
    proportion_of_min_trials = 0.75;
    nRandomSubsamples = 50;
manuallyReduceIntervals = true; % test a subset of all the intervals for faster testing
    manualIntervals = {[175 275], [175 350]}; % NB, also need to change variable fname_base to grab file containing desired intervals
kfold_num = 5; % k-fold cross validation. 
% Note that folds are not generated purely randomly -- we impose
% constraints about not re-using the same image in train vs test ("abstract
% category") and about balancing the number of cats and dogs in each fold
% (so that the shuffle comes out right at 50 %).

% Interval parameters (for loading the correct spike counts file)
step = 5;
width = 100;


% Other Parameters
spikeCountPath = 'XMA2/Spike_count_mats';
TE_LOCS = {'anterior', 'middle', 'posterior'};
catg2_ind1 = 261;

% fname_base = sprintf('%%s_allNeurons_step%d_wd%d.mat', step, width); % double %% escapes the first %s
fname_base = sprintf('%%s_allNeurons_variableBin_1.mat'); % contains 75-175, 175-225, 175-275, and 175-350.

%% Create param struct for saving into record

if ~matchInputSizes
    proportion_of_min_pop_size = NaN; 
    proportion_of_min_trials = NaN;
    nRandomSubsamples = 1;
end

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
    'MatchInputBool', matchInputSizes,...
    'MatchInputParams', [proportion_of_min_pop_size, proportion_of_min_trials, nRandomSubsamples],...
    'SessionsUsed', {rSessionsByMonk},...
    'IntervalsUsed', {ints_used},....
    'TrainingImgSet', img_set_names{svm_train_set_idx},...
    'TestingImgSet', img_set_names{svm_test_set_idx});

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
        Monkeys(m).Sessions(sessn).Session_Y_catg = Y;
        
        Y = [];
        for j = 1:length(Monkeys(m).Sessions(sessn).CueInfo)
            Y = vertcat(Y, repelem(Monkeys(m).Sessions(sessn).CueInfo(j).CueID, Monkeys(m).Sessions(sessn).CueInfo(j).NumApp)');
        end
        if length(Y) ~= length(Monkeys(m).Sessions(sessn).Session_Y_catg)
            error('Error compiling session Y''s ')
        end
        Monkeys(m).Sessions(sessn).Session_Y_imageID = Y;
    end
end

%% Get bools of units to use for each array loc
% Don't worry about matching input sizes yet

for m = 1:length(Monkeys)
    for i = 1:length(Monkeys(m).Sessions)
       
        % Catch units with too few spikes
        enoughSpikesBool = [Monkeys(m).Sessions(i).UnitInfo.SpikeNum] > ignoreVal;
        
        % Look at area
        subBool = cell(1, length(rArrayLocs));
        for iLoc = 1:length(rArrayLocs)
            if strcmp(rArrayLocs{iLoc}, 'te')
                subBool{iLoc} = ismember({Monkeys(m).Sessions(i).UnitInfo.Location}, TE_LOCS);
            else
                subBool{iLoc} = ismember({Monkeys(m).Sessions(i).UnitInfo.Location}, rArrayLocs{iLoc});
            end
            
            % Add enoughSpikes bools
            subBool{iLoc} = subBool{iLoc} & enoughSpikesBool;
        end
        
        % Finalize
        Monkeys(m).Sessions(i).SVMBools = subBool;
    end
end

%% Run SVMs

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    
    % Match input sizes if requested
    if matchInputSizes
        error('match input sizes not yet implemented for img subsets analysis') % IMGSUBSET edit
        [numUnitsVec, numTrials] = getMatchedInputSizes(Monkeys(m).Sessions, ...
            rSessions, rArrayLocs, ...
            proportion_of_min_pop_size, proportion_of_min_trials);
    else 
        nRandomSubsamples = 1; % ie, no subsampling at all
    end
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Get data from appropriate storage place
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf(fname_base, MonkID);
        [X_full, rIntervals] = load_interval_data(fullfile(EXT_HD, spikeCountPath, fileName)); % X is spike counts, rIntervals is list of each interval
        Y = Monkeys(m).Sessions(sessn).Session_Y_catg; % list of categories for each image in X (1 or 2)
        unitBools = Monkeys(m).Sessions(sessn).SVMBools; % cell array of boolean for units to use for each loc
        
        % For finding intervals in X_full
        rIntervals_original = rIntervals;
        
        % Update intervals if requested
        if manuallyReduceIntervals
            rIntervals = manualIntervals;
        end
        
        
        for iLoc = 1:length(rArrayLocs)
            loc = rArrayLocs{iLoc};
            
            % Get boolean for UnitInfo of units in this location
            locBools = unitBools{iLoc};
            
            % Get minimum-adjusted number of units across sessions for this array location.
            if matchInputSizes
                numUnits = numUnitsVec(iLoc); 
            end
            
            % Catch subsets with no data
            if sum(locBools ) == 0
                fprintf('No neurons')
                % do other stuff, store nans
            end
            
            % Take the correct units (not subsetted yet)
            X = X_full(:, locBools, :);
            
            for iInt = 1:length(rIntervals)
                interval = rIntervals{iInt};
                
                % Find index in X's 3rd dim for requested interval
                int_idx = find(cellfun(@(a) all(a == interval), rIntervals_original));
                
                % Pre-allocate the storage vector for the SVMs.
                kflValues = zeros(nRandomSubsamples*kfold_num,1);
                iK = 1;
                if runShuffle
                    kflValues_SHUFFLE = zeros(nRandomSubsamples*kfold_num, nShuffles);
                end
                
                
                for j = 1:nRandomSubsamples
                    
                    % Create subsets of data to use on this loc
                    % iteration. (Or, if not matching, then only have one
                    % loc, and just rename the vars for ease.)
                    if matchInputSizes
                        trial_idx = randsample(size(X,1), numTrials);
                        unit_idx = randsample(size(X,2), numUnits);
                        X_subset = X(trial_idx, unit_idx, int_idx); % data for only those trials and a random loc of neurons
                        Y_subset = Y(trial_idx); % ground truth category for those trials
                        imgID_Subset = Monkeys(m).Sessions(sessn).Session_Y_imageID(trial_idx); % image IDs for those trials
                    else
                        X_subset = X(:,:, int_idx);
                        Y_subset = Y;
                        imgID_Subset = Monkeys(m).Sessions(sessn).Session_Y_imageID;
                    end
                    
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
                        balancedKFoldList(1:(num_cues/2)) = makeKFoldList(num_cues/2, kfold_num);
                        balancedKFoldList((num_cues/2 + 1):end) = makeKFoldList(num_cues/2, kfold_num);
                        
                        % Finally, run the SVM with these folds as the
                        % cross-validation folds.
                        for k = 1:kfold_num 

                            % Which images to use as training in this fold.
                            % IMGSUBSET edit
                            trainingImgs = img_sets{svm_train_set_idx}(balancedKFoldList ~= k); 

                            % Asserts to make sure we did it right
                            assert(length(imgID_Subset) == size(X,1)) % Number of trials should match.
                            assert(sum(ismember(trainingImgs, img_sets{3-svm_train_set_idx}))==0) % Should not be any imgs from the other set
                            
                            % Which trials to use as training, based on the
                            % imgs.
                            idx = ismember(imgID_Subset, trainingImgs); 

                            % Run the model.
                            model = fitclinear(X_subset(idx,:), Y_subset(idx),'Regularization', 'lasso');

                            % Get prediction error.
                            preds = predict(model, X_subset(~idx,:));
                            err = sum(preds ~= Y_subset(~idx)) / sum(~idx);
                            kflValues(iK) = err;

                            % Run with shuffled values if requested
                            if runShuffle
                                for iShuff = 1:nShuffles
                                    Y_subset_SHUFFLE = Y_subset(randperm(length(Y_subset)));
                                    model_SHUFFLE = fitclinear(X_subset(idx,:), Y_subset_SHUFFLE(idx),'Regularization', 'lasso');
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
                            model = fitclinear(X_subset(svm_train_bool,:), Y_subset(svm_train_bool),'Regularization', 'lasso');
                            % Get prediction error.
                            preds = predict(model, X_subset(svm_test_bool,:));
                            err = sum(preds ~= Y_subset(svm_test_bool)) / sum(svm_test_bool);
                            kflValues(iK) = err;
                            % Iterate.
                            iK = iK+1;
                        end
                    end
                end
                
                % Store data.
                kflID = get_good_interval_name2(interval, loc, 'KFL');
                Monkeys(m).Sessions(sessn).(kflID) = kflValues;
                if runShuffle
                    kflID_S = get_good_interval_name2(interval, loc, 'KFL_SHUFFLE');
                    Monkeys(m).Sessions(sessn).(kflID_S) = kflValues_SHUFFLE;
                end
                
                % Report progress.
                fprintf('Done with %s loc %s interval # %d session %d \n', Monkeys(m).Name, loc, iInt, sessn)
            end
        end
    end
end

%% Save the data

% Remove unecessary fields
for m = 1:length(Monkeys)
    Monkeys(m).Sessions = rmfield(Monkeys(m).Sessions, {'UnitInfo', 'CueInfo', 'TrialInfo', 'TimesBT', 'CodesBT', 'Fixn_err_TC'});
end

% Store intervals used
Monkeys(m).RIntervals = rIntervals;

% Save the data and add a row to the Record
fullSVMPath = fullfile(EXT_HD, pv_path, 'SVM_results_%g.mat');
fullRecordPath = fullfile(EXT_HD, svmRecordPath);
save_SVM_data(Monkeys, paramStruct, fullSVMPath, fullRecordPath);
fprintf('Data saved. \n')

%% Functions

function kFoldList = makeKFoldList(len, k)
    subSizes = diff(round(linspace(0, len, k+1)));
    regions = repelem(1:k, subSizes);
    kFoldList = regions(randperm(len));
end

function [numUnits, numTrials] = getMatchedInputSizes(Sessions, rSessions, rArrayLocs, U, T)
trialNums = zeros(length(rSessions), 1);
unitNums = zeros(length(rSessions), length(rArrayLocs));
for i = 1:length(rSessions)
    sessn = rSessions(i);
    trialNums(i) = sum([Sessions(sessn).CueInfo.NumApp]);
    unitNums(i, :) = cellfun(@sum, Sessions(sessn).SVMBools);
    numTrials = ceil(min(trialNums)*T); % a single number
    numUnits = ceil(min(unitNums)*U); % length(rArrayLocs) numbers, takes U*min of each loc across sessions
end
end
