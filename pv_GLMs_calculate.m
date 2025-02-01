% script for running GLM analyses

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

%% Load record if desired

% ID = 633648; % 75-175, 175-275, 275-375
% ID = 748804; % all days, same 3 ints
% ID = 224797; % same as above but with GLM coeffs

% ==== revisions ====
% ID: 21110 -- pre / post, 175-275, matched trial nums. no shuffle.

% to look at what has already been run
fullRecordPath = fullfile(EXT_HD, glmRecordPath);
load(fullRecordPath)


%% Set parameters

% Model predictors
include_trial_metadata = false; % include loc in trial and trial num in regression?
if include_trial_metadata
    numCoeffRows = 7; % in the model with trial metadata, relevant outputs are in rows 2-7
else
    numCoeffRows = NaN; % otherwise, it's just intercept and slope, and we want slope, ie row 2
end

% GLM Parameters
random_seed = 10; % for reproducibility 
% rSessionsByMonk = {[7 9], [6 7]};
rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};
ignoreVal = 20; % if neuron has less than this num spikes, do not use it.
runShuffle = false; % run the shuffled condition?
    nShuffles = 100;

matchedTrialNumPrePost = true; % if true, match num trials (per catg) used in the classifier pre/post training
fraction_min_num_trials = 0.95;  % what fraction of the min number of trials to use for each classifier (otherwise you get no error bars for the session with the minimum number)

manual_intervals = {[175 275]};

% Interval parameters
step = 5;
width = 100;

% Other Parameters
spikeCountPath = 'XMA2/Spike_count_mats';
TE_LOCS = {'anterior', 'middle', 'posterior'};
catg2_ind1 = 261;
loc_in_trials_fname_base = '%s_imgLocInTrials.mat';
trial_num_fname_base = '%s_trialNum.mat';

%% Create param struct for saving into record

paramStruct = struct('RandomSeed', random_seed, ...
    'IgnoreVal', ignoreVal,...
    'RunShuffle', runShuffle, ...
    'SessionsUsed', {rSessionsByMonk},...
    'IncludeTrialMetadata', include_trial_metadata);

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

%% Run GLMs

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    
   % Calculate num trials to use if matching
    if matchedTrialNumPrePost
        num_trials_per_catg_to_use = inf;
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            Y = Monkeys(m).Sessions(sessn).Session_Y_catg;
            min_catg_presentations = min([sum(Y==1) sum(Y==2)]);
            if min_catg_presentations < num_trials_per_catg_to_use
                num_trials_per_catg_to_use = min_catg_presentations;
            end
        end
    end
    num_trials_per_catg_to_use = round(fraction_min_num_trials * num_trials_per_catg_to_use);
    disp(num_trials_per_catg_to_use)
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Get spike count data from appropriate storage place
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf('%s_allNeurons_step%d_wd%d.mat', MonkID, step, width);
        [X_full, rIntervals] = load_interval_data(fullfile(EXT_HD, spikeCountPath, fileName)); % X is spike counts, rIntervals is list of each interval
        Y = Monkeys(m).Sessions(sessn).Session_Y_catg; % list of categories for each image in X (1 or 2)
        
        % Get other regressors we might use
        fileName = sprintf(loc_in_trials_fname_base, MonkID);
        loc_in_trial_mat = load(fullfile(EXT_HD, spikeCountPath, fileName), 'X');
        loc_in_trial_mat = loc_in_trial_mat.X;
        fileName = sprintf(trial_num_fname_base, MonkID);
        trial_num_mat = load(fullfile(EXT_HD, spikeCountPath, fileName), 'X');
        trial_num_mat = trial_num_mat.X;
        
        % Manually reduce intervals tested
        warning('Reducing intervals tested...comment out these lines to run all intervals')
        rIntervals_original = rIntervals; % for finding idx in 3rd dim of X_full
        
        % Pre-allocate the storage vectors
        if include_trial_metadata
            pVals = zeros(size(X_full,2), length(manual_intervals ), numCoeffRows);
            coeffs = zeros(size(X_full,2), length(manual_intervals ), numCoeffRows);
        else
            pVals = zeros(size(X_full,2), length(manual_intervals ));
            coeffs = zeros(size(X_full,2), length(manual_intervals ));
        end
        
        if runShuffle && include_trial_metadata
            error('shuffle not implemented for trial metadata yet')
        elseif runShuffle
            pVals_SHUFFLE = zeros(size(X_full,2), length(manual_intervals ), nShuffles);
            coeffs_SHUFF = zeros(size(X_full,2), length(manual_intervals ));
        end
                
        for iUnit = 1:size(X_full,2)
            
            %Catch units with too few spikes
            if Monkeys(m).Sessions(sessn).UnitInfo(iUnit).SpikeNum < ignoreVal
                pVals(iUnit, :) = NaN;
                if runShuffle
                    pVals_SHUFFLE(iUnit, :, :) = NaN;
                end
                continue
            end
            
            for iInt = 1:length(manual_intervals)
                
                % Find the desired interval in the data
                interval = manual_intervals {iInt};
                interval_idx = find(cellfun(@(a) all(a == interval), rIntervals_original));
                
                % Match num trials across sessions
                trial_idx = true(length(Y),1);
                if matchedTrialNumPrePost
                    cat_inds = find(Y == 1);
                    dog_inds = find(Y == 2);
                    cat_inds_to_keep = randsample(cat_inds, num_trials_per_catg_to_use);
                    dog_inds_to_keep = randsample(dog_inds, num_trials_per_catg_to_use);
                    trial_idx(cat_inds(~ismember(cat_inds, cat_inds_to_keep))) = 0;
                    trial_idx(dog_inds(~ismember(dog_inds, dog_inds_to_keep))) = 0;
                end
                assert(sum(Y(trial_idx)==1) == sum(Y(trial_idx)==2))

                if include_trial_metadata
                    t = table(X_full(trial_idx, iUnit, interval_idx), loc_in_trial_mat, trial_num_mat, Y(trial_idx),...
                        'VariableNames', {'spikeCount', 'locInTrial', 'trialNum', 'imgCatg'});
                    formula = 'spikeCount ~ locInTrial * trialNum * imgCatg - locInTrial:trialNum:imgCatg';
                    glm = fitglm(t, formula, 'Distribution', 'poisson');
                    pVals(iUnit, iInt, :) = glm.Coefficients{:,4}; % 4th col is pval
                    coeffs(iUnit, iInt, :) = glm.Coefficients{:,1}; % 1st col is estimate
                else
                    t = table(X_full(trial_idx, iUnit, interval_idx), Y(trial_idx), ....
                        'VariableNames', {'spikeCount', 'imgCatg'});
                    glm = fitglm(t, 'spikeCount ~ imgCatg', 'Distribution', 'poisson');
                    pVals(iUnit, iInt) = glm.Coefficients{2,4}; % 2,4 is ind for pval of slope
                    coeffs(iUnit, iInt) = glm.Coefficients{2,1}; % estimate of slope
                end
                
                % Run with shuffled values if requested
%                 if runShuffle
%                     for iShuff = 1:nShuffles
%                         t = table(X_full(:, iUnit), Y(randperm(length(Y))));
%                         glm = fitglm(t, 'SC ~ catg', 'Distribution', 'poisson');
%                         pVals_SHUFFLE(iUnit, iInt, iShuff) = glm.Coefficients{2,4};
%                         coeffs(iUnit, iInt, iShuff) = glm.Coefficients{2,1}; % estimate of slope
%                     end
%                 end
            end
            
            % Report progress.
            fprintf('Done with %s unit # %d session %d \n', Monkeys(m).Name, iUnit, sessn)
        end 
        
        % Store data.
        Monkeys(m).Sessions(sessn).GLM_Pvals = pVals;
        Monkeys(m).Sessions(sessn).GLM_coeffs = coeffs;
        Monkeys(m).Sessions(sessn).GLM_intervals = manual_intervals;
        
        
        if runShuffle
            Monkeys(m).Sessions(sessn).GLM_Pvals_SHUFFLE = pVals_SHUFFLE;
            Monkeys(m).Sessions(sessn).GLM_coeffs_SHUFFLE = coeffs_SHUFFLE;
        end

    end
end

%% Save the data

% Remove unecessary fields
for m = 1:length(Monkeys)
    Monkeys(m).Sessions = rmfield(Monkeys(m).Sessions, {'UnitInfo', 'CueInfo', 'TrialInfo', 'TimesBT', 'CodesBT', 'Fixn_err_TC'});
end

% Save the data and add a row to the Record
fullGLMPath = fullfile(EXT_HD, pv_path, 'GLM_results_%g.mat');
fullRecordPath = fullfile(EXT_HD, glmRecordPath);
save_SVM_data(Monkeys, paramStruct, fullGLMPath, fullRecordPath);
fprintf('File saved \n')



