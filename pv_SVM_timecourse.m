% run SVMs on each interval of data from the passive viewing tasks

%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';

% Load behavioral data
load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_behav_and_metaNI.mat')) % behavior and neural summaries, but w/o spike times

%% Run once -- add unit info (except actual spike count lists) to behav
% 
% % Load the data, broken up into sessions for storage
% Data = load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_ni.mat'));
% 
% % Stitch the sessions back into one big structure
% [status, MU] = stitch_monkeyStruct_from_parts(Data);
% clear Data
% 
% for m = 1:length(Monkeys)
%     for i = 1:length(Monkeys(m).Sessions)
%         Monkeys(m).Sessions(i).UnitInfo = MU(m).Sessions(i).UnitInfo;
%         for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
%             Monkeys(m).Sessions(i).UnitInfo(j).SpikeNum = length(Monkeys(m).Sessions(i).UnitInfo(j).Spike_times);
%         end
%         Monkeys(m).Sessions(i).UnitInfo = rmfield(Monkeys(m).Sessions(i).UnitInfo, 'Spike_times');
%     end
% end
% clear MU

%% Set parameters

% SVM Parameters
random_seed = 10; % for reproducibility 
rSubsets = {'te', 'anterior', 'middle', 'posterior'}; % relevant subsets
rSessions = {[7 9], [6 7]};
ignoreVal = 20; % if neuron has less than this num spikes, do not use it.
run_shuffle = false; % run the shuffled condition?
matchInputSizes = false; % make input matrices all same size (control) ?
    % These parameters only used if matching input sizes
    proportion_of_min_pop_size = 0.75; 
    proportion_of_min_trials = 0.75;
    n_subsets = 20;
kfold_num = 5; % k-fold cross validation. 
% Note that folds are not generated purely randomly -- we impose
% constraints about not re-using the same image in train vs test ("abstract
% category") and about balancing the number of cats and dogs in each fold
% (so that the shuffle comes out right at 50 %).

% Other Parameters
spikeCountPath = 'XMA2/Spike_count_mats';
TE_LOCS = {'anterior', 'middle', 'posterior'};
catg2_ind1 = 261;

%% Create param struct for saving into record
paramStruct = struct('RandomSeed', random_seed, 'RSubsets', {rSubsets}, ...
    'IgnoreVal', ignoreVal, 'RunShuffle', run_shuffle, 'KFoldNum', kfold_num, ...
    'MatchInputBool', matchInputSizes,...
    'MatchInputParams', [proportion_of_min_pop_size, proportion_of_min_trials, n_subsets],...
    'SessionsUsed', {rSessions});


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

%% Get bools of units to use for each subset
% Don't worry about matching input sizes yet

for m = 1:length(Monkeys)
    for i = 1:length(Monkeys(m).Sessions)
       
        % Catch units with too few spikes
        enoughSpikesBool = [Monkeys(m).Sessions(i).UnitInfo.SpikeNum] > ignoreVal;
        
        % Look at area
        subBool = cell(1, length(rSubsets));
        for iSub = 1:length(rSubsets)
            if strcmp(rSubsets{iSub}, 'te')
                subBool{iSub} = ismember({Monkeys(m).Sessions(i).UnitInfo.Location}, TE_LOCS);
            else
                subBool{iSub} = ismember({Monkeys(m).Sessions(i).UnitInfo.Location}, rSubsets{iSub});
            end
            
            % Add enough spikes data
            subBool{iSub} = subBool{iSub} & enoughSpikesBool;
        end
        
        % Finalize
        Monkeys(m).Sessions(i).SVMBools = subBool;
    end
end


%% Run SVMs

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = 1:length(Monkeys(m).Sessions);
    
    % Match input sizes if requested
    [numUnits, numTrials] = getMatchedInputSizes(Monkeys(m).Sessions, ...
        rSessions, rSubsets, ...
        proportion_of_min_pop_size, proportion_of_min_trials);
                
    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Get data from appropriate storage place
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf('%s_allNeurons.mat', MonkID);
        [X, rIntervals] = load_interval_data(fullfile(EXT_HD, spikeCountPath, fileName));
        Y = Monkeys(m).Sessions(sessn).Session_Y_catg;
        
        for iSub = 1:length(rSubsets)
            subset = rSubsets{iSub};
            unitBools = Monkeys(m).Sessions(i).SVMBools;
            
            for p = 1:length(rIntervals)
                interval = rIntervals{p};

                
                % avoid empty data
                if isempty(X)
                    fprintf('No neurons for session %d, interval %d to %d, subset %s\n',...
                        sessn, interval(1), interval(2), neuron_subset)
                    fprintf('Skipping and setting KFL to NaN\n')
                    kfl_id = get_good_interval_name2(interval, neuron_subset, 'KFL_abcat');
                    Monkeys(m).Sessions(sessn).(kfl_id) = repelem(NaN, n_subsets)';
                    continue
                end

                % now select not only random subsets of neurons
                % to use, but random and non-overlapping subsets of images
                % to use as train vs test data.
                kfl_vec = zeros(n_subsets*kfold_num,1);
                iK = 1;
                for j = 1:n_subsets
    %                 trial_idx = randsample(size(X,1), num_trials_to_use); % which trials to use
    %                 X_subset = X(trial_idx,randsample(size(X,2),nn)); % data for only those trials and a random subset of neurons
    %                 Y_subset = Y(trial_idx); % ground truth category for those trials
    %                 img_ids = Monkeys(m).Sessions(sessn).Session_Img_IDs(trial_idx); % image IDs for those trials
                    X_subset = X(:,randsample(size(X,2),nn)); % data for all trials and a random subset of neurons
                    Y_subset = Y(:); % ground truth category for those trials
                    img_ids = Monkeys(m).Sessions(sessn).Session_Y_imageID(:); % (1xnum trials) of image IDs 
                    num_cues = length(unique(Monkeys(m).Sessions(sessn).Session_Y_imageID)); % how many cues total exist
                    
                    % get balanced image subsets
                    balanced_kfold_list = zeros(1, num_cues);
                    balanced_kfold_list(1:catg2_ind1-1) = crossvalind('Kfold', num_cues/2, kfold_num);
                    balanced_kfold_list(catg2_ind1:end) = crossvalind('Kfold', num_cues/2, kfold_num);
                    kfold_list = balanced_kfold_list;
                    
                    for k = 1:kfold_num % use each fold as held-out testing set and average the prediction error
                        train_imgs = find(kfold_list ~= k); % vector of cue IDs to be used as training cues
                        idx = ismember(img_ids, train_imgs); % boolean vector of size (1x num trials) of which trials to be used as training
                        model = fitclinear(X_subset(idx,:), Y_subset(idx),'Regularization', 'lasso');
                        preds = predict(model, X_subset(~idx,:));
                        err = sum(preds ~= Y_subset(~idx)) / sum(~idx);
                        kfl_vec(iK) = err;
                        iK = iK+1;
                    end
                end
                % store results
                kfl_id = get_good_interval_name2(interval, neuron_subset, 'KFL_abcat');
                Monkeys(m).Sessions(sessn).(kfl_id) = kfl_vec;
                
                if ~run_shuffle
                    fprintf('Done with %s subset %s interval # %d session %d (no shuffle) \n', Monkeys(m).Name, neuron_subset, p, sessn)
                    continue
                end
                
                % ditto but shuffled
                kfl_vec_shuffled = zeros(n_subsets*kfold_num,1);
                iK = 1;
                for j = 1:n_subsets
    %                 trial_idx = randsample(size(X,1), num_trials_to_use); % which trials to use
    %                 X_subset = X(trial_idx,randsample(size(X,2),nn)); % data for only those trials and a random subset of neurons
    %                 Y_subset = Y(trial_idx); % ground truth category for those trials
    %                 img_ids = Monkeys(m).Sessions(sessn).Session_Img_IDs(trial_idx); % image IDs for those trials
                    X_subset = X(:,randsample(size(X,2),nn)); % data for all trials and a random subset of neurons
                    Y_subset = Y(:); % ground truth category for those trials
                    img_ids = Monkeys(m).Sessions(sessn).Session_Y_imageID(:); % (1xnum trials) of image IDs 

                    % === SHUFFLE ===
                    shuff_inds = randperm(length(Y_subset));
                    Y_subset = Y_subset(shuff_inds);
                    img_ids = img_ids(shuff_inds);
                    
                    num_cues = length(unique(Monkeys(m).Sessions(sessn).Session_Y_imageID)); % how many cues total exist
                    
                    % OLD WAY for kfold list
%                     kfold_list = crossvalind('Kfold', num_cues, kfold_num); % not actualyl using this for crossvalidating. vector saying what fold each cue is in, regardless of trial_idx
                    
                    % NEW WAY for kfold list
                    balanced_kfold_list = zeros(1, num_cues);
                    balanced_kfold_list(1:catg2_ind1-1) = crossvalind('Kfold', num_cues/2, kfold_num);
                    balanced_kfold_list(catg2_ind1:end) = crossvalind('Kfold', num_cues/2, kfold_num);
                    kfold_list = balanced_kfold_list;
                    
                    for k = 1:kfold_num % use each fold as held-out testing set and average the prediction error
                        train_imgs = find(kfold_list ~= k); % vector of cue IDs to be used as training cues
                        idx = ismember(img_ids, train_imgs); % boolean vector of size (1x num trials) of which trials to be used as training
                        model = fitclinear(X_subset(idx,:), Y_subset(idx),'Regularization', 'lasso');
                        preds = predict(model, X_subset(~idx,:));
                        err = sum(preds ~= Y_subset(~idx)) / sum(~idx);
                        kfl_vec_shuffled(iK) = err;
                        iK = iK + 1;
                    end
                end
                
                % store shuffled results
                kfl_shuff_id = get_good_interval_name2(interval, neuron_subset, 'KFL_abcat_SHUFFLED');
                Monkeys(m).Sessions(sessn).(kfl_shuff_id) = kfl_vec_shuffled;
                fprintf('Done with %s subset %s interval # %d session %d \n', Monkeys(m).Name, neuron_subset, p, sessn)
            end
        end
    end
end

%% Functions

function [numUnits, numTrials] = getMatchedInputSizes(Sessions, rSessions, rSubsets, U, T)
trialNums = zeros(length(rSessions), 1);
unitNums = zeros(length(rSessions), length(rSubsets));
for i = 1:length(rSessions)
    sessn = rSessions(i);
    trialNums(i) = sum([Sessions(sessn).CueInfo.NumApp]);
    unitNums(i, :) = cellfun(@sum, Sessions(sessn).SVMBools);
    numTrials = ceil(min(trialNums)*T); % a single number
    numUnits = ceil(min(unitNums)*U); % length(rSubsets) numbers, takes U*min of each subset across sessions
end
end
