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

% Use this for car/truck
% svmRecordPath = 'XMA2/Monkey_structs/SVM_Records.mat';
% pv_path = 'XMA2/Monkey_structs';
% spikeCountPath = 'XMA2/Spike_count_mats_CARTRUCK';  % car truck skip
% behav_file = 'MaxMarta_xma2_CARTRUCK_ni.mat';
% short_name_regexp = '\w*cars\w*';

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
% num_best_units_to_use = [1:2:50 52:4:200];
% num_best_units_to_use = [50 100 110 120 150 175];
% num_best_units_to_use = [10 25 50 100];
num_best_units_to_use = [100];  % defaulting to 100 for timecourse too now, because MATLAB's SPARSA solver does better at 100 than w more.

ranking_interval = {[175 275]};

% Might as well just test each image set...
% "train" --> 20/20 set. "test" --> 240/240 set.
img_sets = {[1:20 261:280], [21:260 281:520], [1:260, 261:520]}; % {[cats_set1 dogs_set1], [cats_set2...]}
img_set_names = {'train', 'test', 'all'};

% Pick here which img sets to train/test on. Will need to re-run this
% script for each desired combo. (Too many for-loops deep already...easier
% to do this manually.)
svm_train_set_idx = 3; % 1 for 20/20 training set, 2 for 240/240 testing set, 3 for all
svm_test_set_idx = 3;


% General SVM Parameters
random_seed = 10; % for reproducibility 
rArrayLocs = {'te'};
% rArrayLocs = {'anterior', 'middle', 'posterior'};
rSessionsByMonk = {[7 9], [6 7]};  % (Fig 4A / 4B)
% rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};  % Fig 4C
% rSessionsByMonk = {[], 1:7};  % tmp hotfix to just run max
% rSessionsByMonk = {[1 2 3 5 6 7 9], [1 2 6 7]};  % Fig 4C with matched num trials (base 3-5 of Max have super low trial counts)
% rSessionsByMonk = {[10 11 12], [8 9 10]};  % Fig ??, cat/dog sessions bookending the car/truck training

% rSessionsByMonk = {[], [4 5]};  % car/truck learning
% rSessionsByMonk = {1:7, 1:5};  % car/truck baselines

matchedTrialNumPrePost = true; % if true, match num trials (per catg) used in the classifier pre/post training
fraction_min_num_trials = 0.85;  % what fraction of the min number of trials to use for each classifier (otherwise you get no error bars for the session with the minimum number)

ignoreVal = 20; % if neuron has less than this num spikes, do not use it.
runShuffle = false; % run the shuffled condition? need to run it for the timecourse clustering stats
    nShuffles = 20;
    
manuallyReduceIntervals = true; % test a subset of all the intervals for faster testing
%     manualIntervals = {[75 175], [175 275]}; % NB, also need to change variable fname_base to grab file containing desired intervals
    manualIntervals = {[175 275]}; % NB, also need to change variable fname_base to grab file containing desired intervals
    X_fname_base = sprintf('%%s_allNeurons_variableBin_1.mat'); % contains 75-175, 175-225, 175-275, and 175-350.
    
    % for SVM timecourses
%     starts = -100:10:300;
%     step = 5; % Interval parameters (for loading the correct spike counts file)
%     width = 100;
%     manualIntervals = arrayfun(@(x) [x x+100], starts, 'UniformOutput', false);
%     ranking_interval = {[175 275]};  % interval to rank units wrt one-D SVM's
%     X_fname_base = sprintf('%%s_allNeurons_step%d_wd%d.mat', step, width); % contains full time-course of spike counts 


% Note that folds are not generated purely randomly -- we impose
% constraints about not re-using the same image in train vs test ("abstract
% category") and about balancing the number of cats and dogs in each fold
% (so that the shuffle comes out right at 50 %).
kfold_num = 5; % k-fold cross validation. 

% Other Parameters
TE_LOCS = {'anterior', 'middle', 'posterior'};
catg2_ind1 = 261;

%% Count total num appearances of each image

for m = 1:length(Monkeys)
    img_counts = [];
    rSessions = rSessionsByMonk{m};
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        if contains(Monkeys(m).Sessions(sessn).ShortName, "UNI")
            continue
        end
        img_counts = [Monkeys(m).Sessions(sessn).CueInfo.NumApp];
        disp(img_counts(1:10))
        fprintf("Monkey %s, session %s, %0.2f +/- %0.2f pres / img \n", ...
            Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName, ...
            mean(img_counts), std(img_counts))
    
    end
    
end

%% XXX Load GLM data

% Load record

glmRecordPath = 'XMA2/Monkey_structs/GLM_Records.mat';
pv_path = 'XMA2/Monkey_structs';
fullGLMPath = fullfile(EXT_HD, pv_path, 'GLM_results_%g.mat');


% ID = 633648; % 75-175, 175-275, 275-375
% ID = 748804; % all days, same 3 ints
ID = 224797; % same as above but with GLM coeffs

load(sprintf(fullGLMPath, ID), 'data');
GLM = data;
clear data

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
    'TrainingImgSet', img_set_names{svm_train_set_idx},...
    'TestingImgSet', img_set_names{svm_test_set_idx},...
    'NumBestUnits', {num_best_units_to_use},...
    'MatchedTrialNumPrePost', matchedTrialNumPrePost);

%% Record of ID nums so far

% ID: 771321 -- 3:3:100, {[75 175], [175 275]}, pre/post, ranked wrt GLM coeffs.
% Later thought: it's nice that the analysis still works with GLM coeffs,
% but the "best" neurons wrt the GLM coeff rankings aren't very good
% compraed to the "best" ones found wrt the one-D SVMs. So this probably
% just doesn't capture what we're interested in and we should use SVM
% rankings instead.

% DONT USE
% ID: 748804 -- 1:1:100 and 1:5:300 {[175 275]}, pre/post, ranked wrt single-unit
% SVMs. Really wonky looking results. Oops scratcht that -- did kfl instead
% of acc! lol.
% DONT USE

% ID: 498507 -- 1:1:200, {[175 275]}, pre/post. but with weird switch betw sgd and
% sparsa SVM solvers.

% ID: 224797 -- num_best_units_to_use = [1:1:50 52:2:200]; {[175 275]}, pre/post.
% Still weird sga / sparsa switch, for now seemingly nothing to be done
% about it.

% ID: 198063 -- ^^, but also has analysis where we remove the best units
% one at a time.

% ID: 921995 -- full timecourse (step=10 ms), shuff x5, pre/post, [10 25 50
% 100] unit subsets. (Fig 4A)

% ID: 316514 -- all sessions, [10 25 50 100], {[175 275]}, no shuffle. (Fig
% 4B)

% ID: 723146 -- full timecourse for separate arrays (pre/post, 100 best
% units).

% ID: 257183 -- same as above (723146) but with shuffles.


% *** trial num matched for revisions below here ***

% ID: 331368 all sessions, [100] best units, 175-275, no shuff, matched
% num cats/dogs. (for Fig 4C)

% ID: 291617 -- same as above ^, but with 240/240 images.

% ID: 994359 -- pre/post, [100] best units, full timecourse, no shuff, matched
% nums. For Fig 4A.

% ID: 297691 --same as above ^, but with 20x shuffles.

% ID: 66973 -- same as above ^, but with 240/240 images.

% XXXX ID: 812755 -- same as above, but with 20x shuffles and 240/240.

% ID: 958258 -- pre/post, [100] best units, full timecourse, arrays with
% shuffle.

% ID: 760531 -- pre/post 1:100 best units, matched nums. Fig 4E

% ID: 163919 -- cat/dog sessions that bookend the car/truck training.
% Similar to Fig 4A, [100] best units, 175-275, matched num trials. Odd
% results — Marta's decoding acc falls to zero?

% ===============
% CARTRUCK data
% ID: 326003 -- car/truck PF's. Marta post (#7) is bad due to array damage, it seems.

% ID: 643100 -- car/truck PF's with baselines days, matched num trials.

% ID: 377431 -- car/truck timecourse, pre/post, matched, shuffle.
% ===============

% Load record
fullRecordPath = fullfile(EXT_HD, svmRecordPath);
load(fullRecordPath, 'Record');

%% (skip this if analyzing anew) Find desired SVM data

ID = 812755;

fNames = fields(Record);
row = find([Record.ID] == ID);
for f = 1:length(fNames)
    s = Record(row).(fNames{f});
    field = fNames{f};
    if strcmp(field, 'MatchInputParams')
        fprintf('%s: %s \n', field, mat2str(s))
    elseif strcmp(field, 'SessionsUsed')
        fprintf('%s: %s--%s, %s--%s \n', field, Monkeys(1).Name, mat2str(s{1}), Monkeys(2).Name, mat2str(s{2}))
    elseif isa(s, 'cell')
        fprintf('%s: %s \n', field, join(cellfun(@string ,s), ', '))
    else
        fprintf('%s: %s \n', field, string(s))
    end
    
end

%% (skip this if analyzing anew) (but useful to load in pre-ranked units) Load SVM data
% ID = 771321;
% ID = 498507;
% ID = 198063;
% ID = 224797;

% ID = 921995;  % Fig 4A (timecourse)
% ID = 316514;  % Fig 4B (all sessions)
% ID = 198063;  % Fig 4C (adding units)
% ID = 257183; % full timecourse for separate arrays (pre/post, 100 best units).

% revisions
ID = 291617; % 240/240
% ID = 207437;
% ID = 994359;  % timecourse
% ID = 331368;  % Fig 4B (all sessions), matched
% ID = 760531;  % increasing num units
% ID = 958258;  % timecourse w array shuffles


% ID = 326003;  % car-truck, not matched tr nums
% ID = 643100;  % car-truck, baseline + matched

load(sprintf(fullSVMPath, ID), 'data');
SVM = data;
clear data

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
        SVM(m).Sessions(sessn).Session_Y_catg = Y;
        
        Y = [];
        for j = 1:length(Monkeys(m).Sessions(sessn).CueInfo)
            Y = vertcat(Y, repelem(Monkeys(m).Sessions(sessn).CueInfo(j).CueID, Monkeys(m).Sessions(sessn).CueInfo(j).NumApp)');
        end
        if length(Y) ~= length(SVM(m).Sessions(sessn).Session_Y_catg)
            error('Error compiling session Y''s ')
        end
        SVM(m).Sessions(sessn).Session_Y_imageID = Y;
    end
end

%% (don't use) Rank units based on GLM coeffs

% The GLM coefficients are tempting as a metric for category coding, but
% they're actually kind of hard to interpret. You can have high coeffs due
% to outliers; and lots of images being a little different or a few images
% being a lot different could end up with the same coefficient; so it's not
% actually informative wrt category decodability.

% Monkeys(m).Sessions(i).GLM_coeffs is (num units) x (num intervals used)
for m = 1:length(Monkeys)
    sessions_to_use = rSessionsByMonk{m};
    
    for i = 1:length(sessions_to_use)
        session = sessions_to_use(i);
        
        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};
            glm_intervals_calcd = GLM(m).Sessions(session).GLM_intervals;
            idx = find(cellfun(@(a) all(a == interval), glm_intervals_calcd));
            if isempty(idx)
                error('No GLM coeffs found for interval %d - %d', interval(1), interval(2));
            end
            
            % glm coeffs should exist for all units in the data
            coeffs = GLM(m).Sessions(session).GLM_coeffs(:, idx);
            assert (length(coeffs)==length(Monkeys(m).Sessions(session).UnitInfo))
            
            
            % sort each array subset separately for convenience
            for iLoc = 1:length(rArrayLocs)
                array = rArrayLocs{iLoc};
                if strcmp(rArrayLocs{iLoc}, 'te')
                    units_to_use_bool = ismember({Monkeys(m).Sessions(session).UnitInfo.Location}, TE_LOCS);
                    units_to_use_list = find(units_to_use_bool);
                else
                    units_to_use_bool = ismember({Monkeys(m).Sessions(session).UnitInfo.Location}, array);
                    units_to_use_list = find(units_to_use_bool);
                end
                
                field_name = get_good_interval_name2(interval, array, 'Ranking_GLM_coeff');
                
                % rank all units in descending order, ie ranking(i) gives
                % the ranking of unit i, and find(ranking==1) gives the
                % unit number of the unit with the highest abs coeff.
                ranking = argsort(flipud(argsort(abs(coeffs))));  % argsorting twice is quick way to rank!
                
                % only include units that we're using
                ranking_only_included_units = ranking(units_to_use_bool);
                
                % re-rank
                ranking_only_included_units = argsort(argsort(ranking_only_included_units));
                
                result = zeros(size(coeffs,1),1);
                result(units_to_use_bool) = ranking_only_included_units;
                result(~units_to_use_bool) = NaN;
                
                % sanity check
                ind_of_highest_coeff = find(result==1);
                assert(max(abs(coeffs(units_to_use_bool))) == abs(coeffs(ind_of_highest_coeff)));
                
                for j = 1:length(result)
                    SVM(m).Sessions(session).UnitInfo(j).(field_name) = result(j);
                end
            end
        end
    end
end

%% (run once, then save output) Calculate single-unit SVM performance

% does not match num trials across sessions, since doesn't matter for
% ranking within a session.

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);

         % TMP hotfix:
%         if ~isempty(SVM(m).Sessions(sessn).UnitInfo)
%             continue
%         end

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
            
%             for iInt = 1:length(manualIntervals)
            for iInt = 1:length(ranking_interval)
%                 interval = manualIntervals{iInt};
                interval = ranking_interval{iInt};
                
                % Find index in X's 3rd dim for requested interval
                int_idx = find(cellfun(@(a) all(a == interval), rIntervals_original));
                
                for iUnit = 1:length(Monkeys(m).Sessions(sessn).UnitInfo)
                    units_to_use = [iUnit];
                    
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
                            model = fitclinear(X_subset(idx,:), Y_subset(idx),'Regularization', 'lasso');

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
                            model = fitclinear(X_subset(svm_train_bool,:), Y_subset(svm_train_bool),'Regularization', 'lasso');
                            % Get prediction error.
                            preds = predict(model, X_subset(svm_test_bool,:));
                            err = sum(preds ~= Y_subset(svm_test_bool)) / sum(svm_test_bool);
                            kflValues(iK) = err;
                            % Iterate.
                            iK = iK+1;
                        end
                    end
                    
                    % Store data.
                    kflID = get_good_interval_name2(interval, array, 'KFL_SingleUnitSVMs');
                    SVM(m).Sessions(sessn).UnitInfo(iUnit).(kflID) = kflValues;
    
                    % Report progress.
                    fprintf('Done with %s loc %s interval # %d, unit %d, session %d \n', Monkeys(m).Name, array, iInt, iUnit, sessn)
                end
            end
        end
    end
end

%% Rank units based on single-unit SVM performance

% rArrayLocs = {'anterior', 'middle', 'posterior',};
% rArrayLocs = {'te'};

% Monkeys(m).Sessions(i).GLM_coeffs is (num units) x (num intervals used)
for m = 1:length(Monkeys)
    sessions_to_use = rSessionsByMonk{m};
    
    for i = 1:length(sessions_to_use)
        session = sessions_to_use(i);
        
%         for iInt = 1:length(manualIntervals)
        for iInt = 1:length(ranking_interval)
%             interval = manualIntervals{iInt};
            interval = ranking_interval{iInt};
            
            kflID = get_good_interval_name2(interval, 'te', 'KFL_SingleUnitSVMs');
            singleunit_svm_kfls = [SVM(m).Sessions(session).UnitInfo.(kflID)];
            singleunit_mean_svm_kfls = mean(singleunit_svm_kfls);
            
            % sort each array subset separately for convenience
            for iLoc = 1:length(rArrayLocs)
                array = rArrayLocs{iLoc};
                if strcmp(rArrayLocs{iLoc}, 'te')
                    units_to_use_bool = ismember({Monkeys(m).Sessions(session).UnitInfo.Location}, TE_LOCS);
                    units_to_use_list = find(units_to_use_bool);
                else
                    units_to_use_bool = ismember({Monkeys(m).Sessions(session).UnitInfo.Location}, array);
                    units_to_use_list = find(units_to_use_bool);
                end
                
                field_name = get_good_interval_name2(interval, array, 'Ranking_SingleUnitSVM_KFL');
                
                % rank all units in descending order, ie ranking(i) gives
                % the ranking of unit i, and find(ranking==1) gives the
                % unit number of the unit with the highest KFL.
                ranking = argsort(fliplr(argsort(1 - singleunit_mean_svm_kfls)));  % argsorting twice is quick way to rank!
                
                % only include units that we're using
                ranking_only_included_units = ranking(units_to_use_bool);
                
                % re-rank
                ranking_only_included_units = argsort(argsort(ranking_only_included_units));
                
                result = zeros(size(singleunit_mean_svm_kfls,1),1);
                result(units_to_use_bool) = ranking_only_included_units;
                result(~units_to_use_bool) = NaN;
                
                % sanity check
                ind_of_highest = find(result == 1);
                assert(max(1-singleunit_mean_svm_kfls(units_to_use_bool)) == (1-singleunit_mean_svm_kfls(ind_of_highest)));
                
                for j = 1:length(result)
                    SVM(m).Sessions(session).UnitInfo(j).(field_name) = result(j);
                end
            end
        end
    end
end

%% (run once, then save output) Run SVMs, adding more units each time

% To rank by single single-unit SVMs
ranking_field_name = 'Ranking_SingleUnitSVM_KFL';

if ~matchedTrialNumPrePost
    output_field_name_template = 'KFL_SingleUnitSVMRanking_sparsa_Top_%d';
elseif matchedTrialNumPrePost
%     output_field_name_template = 'KFL_SingleUnitSVMRanking_sparsa_Top_%d_TrNMatch';
    output_field_name_template = 'KFL_SUnSVMRk_sprs_T_%d_TrNMch';  % had to shorten for individ arrays and shuffles lol
end

% svm_solver = {'sgd', 'sparsa'};
% svm_solver = {'sparsa'};
svm_solver = 'sparsa';

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
        num_trials_per_catg_to_use = round(fraction_min_num_trials * num_trials_per_catg_to_use * (kfold_num-1)/kfold_num);
    end
    
    
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
                    
                    % TMP hotfix
%                     kflID = get_good_interval_name2(interval, array, sprintf(output_field_name_template, num_best_units_to_use(iUnitSubset)));
%                     if ~isempty(SVM(m).Sessions(sessn).(kflID))
%                         continue
%                     end
                    
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
                            
                            % Further match total number of cat/dog images
                            % if requested.
                            if matchedTrialNumPrePost
                                cat_inds = find(ismember(imgs_shown, trainingImgs) & (Y_subset == 1));
                                dog_inds = find(ismember(imgs_shown, trainingImgs) & (Y_subset == 2));
                                cat_inds_to_keep = randsample(cat_inds, num_trials_per_catg_to_use);
                                dog_inds_to_keep = randsample(dog_inds, num_trials_per_catg_to_use);
                                idx(cat_inds(~ismember(cat_inds, cat_inds_to_keep))) = 0;
                                idx(dog_inds(~ismember(dog_inds, dog_inds_to_keep))) = 0;
                                assert(sum(Y_subset(idx)==1) == sum(Y_subset(idx)==2))
                            end
%                             disp(sum(Y_subset(idx)==1))

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
%                     disp(mean(kflValues))
                    if runShuffle
                        % needed to write "SHUFF" instead here for some
                        % versions because output exceeded the max field
                        % name length lol.
                        kflID_S = get_good_interval_name2(interval, array, sprintf(strcat(output_field_name_template, '_SHUFF'), num_best_units_to_use(iUnitSubset)));
                        SVM(m).Sessions(sessn).(kflID_S) = kflValues_SHUFFLE;
                    end
    
                    % Report progress.
                    fprintf('Done with %s loc %s interval # %d, %d best units, session %d \n', Monkeys(m).Name, array, iInt, num_best_units_to_use(iUnitSubset), sessn)
                end
            end
        end
    end
end

%% Plot performance wrt adding units (Fig 4E)

array = 'te';
rArrayLocs = {'te'};

% base_kfl_name = 'KFL_GLMRanking_Top_%d';
% base_kfl_name = 'KFL_SingleUnitSVMRanking_Top_%d';
% base_kfl_name = 'KFL_SingleUnitSVMRanking_sparsa_Top_%d';  % 2D
% base_kfl_name = 'KFL_SingleUnitSVMRanking_sparsa_Top_%d_TrNMatch';  % Fig 4E
base_kfl_name = 'KFL_SUnSVMRk_sprs_T_%d_TrNMch';

intervals_to_plot = {[175 275]};

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(intervals_to_plot)
            interval = intervals_to_plot{iInt};

            % pre alloc vector to hold mean + std for accuracies
            means = nan(length(rSessions), length(num_best_units_to_use));
            stds = nan(length(rSessions), length(num_best_units_to_use));
            
            
            for i = 1:length(rSessions)
                sessn = rSessions(i);
%                 disp(sessn)
                
                for iUnitSubset = 1:length(num_best_units_to_use)
                    unit_subset = 1:(num_best_units_to_use(iUnitSubset));
                    kflID = get_good_interval_name2(interval, array, sprintf(base_kfl_name, num_best_units_to_use(iUnitSubset)));
                    if ~ismember(kflID, fields(SVM(m).Sessions(sessn)))
                        warning('No KFL field found -- did you use wrong field name template?')
                        continue
                    end
                    kfls = SVM(m).Sessions(sessn).(kflID); 
%                     disp(kfls)

                    means(i, iUnitSubset) = mean(1 - kfls);
                    stds (i, iUnitSubset) = std(1 - kfls);
                end
            end

            % plot the data
            figure
            hold on
            for i = 1:length(rSessions)
                errorbar(num_best_units_to_use, means(i,:), stds(i,:), 'DisplayName', Monkeys(m).Sessions(rSessions(i)).ShortName)
            end
            set(gca, 'XScale', 'log')
            xlabel('Num best units used')
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

% ok this is totally contrived and the fit is very finicky wrt init 
% conditions, let's just do a single sigmoid.
% bisigmoid = fittype('a1/(1 + exp(-b1*(x-c1))) + a2/(1 + exp(-b2*(x-c2))) + d',...
%     'dependent', 'y', 'independent', 'x',...
%     'coefficients', {'a1', 'b1', 'c1', 'a2', 'b2', 'c2', 'd'});

sigmoid = fittype('a1/(1 + exp(-b1*(x-c1))) + d',...
    'dependent', 'y', 'independent', 'x',...
    'coefficients', {'a1', 'b1', 'c1', 'd'});

% manualIntervals = {[175 275]};
% rSessionsByMonk = {[7 9], [6 7]};  
ylims_by_monk = {[0.6 0.8], [0.54, 0.65]};
lower_n_units_val_by_monk = [3 3];

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        
        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};

            
            % Pre alloc vector to hold mean + std for accuracies
            all_means = nan(length(rSessions), length(num_best_units_to_use));
            all_stds = nan(length(rSessions), length(num_best_units_to_use));

            % Retrieve the data
            for i = 1:length(rSessions)
                sessn = rSessions(i);
                for iUnitSubset = 1:length(num_best_units_to_use)
                    unit_subset = 1:(num_best_units_to_use(iUnitSubset));
                    kflID = get_good_interval_name2(interval, array, sprintf(base_kfl_name, num_best_units_to_use(iUnitSubset)));
                    if ~ismember(kflID, fields(SVM(m).Sessions(sessn)))
                        break
                    end
                    kfls = SVM(m).Sessions(sessn).(kflID); 

                    all_means(i, iUnitSubset) = mean(1 - kfls);
                    all_stds (i, iUnitSubset) = std(1 - kfls);
                end
            end
            
            % For simplicity, for now, truncate at 100 units, because
            % that's where the SVM solver gets weird. Thankfully the curves
            % are flat by there anyways.
            idx_to_use = find((num_best_units_to_use <= 100) & (num_best_units_to_use >= lower_n_units_val_by_monk(m)));
            
            % Prep values + pre-alloc for fitting
            means = all_means(:, idx_to_use);
            stds = all_stds(:, idx_to_use);
            xvals = log10(num_best_units_to_use(idx_to_use));
            coeff_means = cell(length(rSessions));
            confint_results = cell(length(rSessions));
            
            % For transparency, show the few data points at begining that
            % don't fit nicely into the sigmoid curve.
            idx_excluded_lower = find((num_best_units_to_use < lower_n_units_val_by_monk(m)));
            excluded_xvals_lower = log10(num_best_units_to_use(idx_excluded_lower));
            excluded_means = all_means(:, idx_excluded_lower);
            
            figure("Position", [400 400 360 300])
            hold on
            yl = [0.48 max(means(:)+stds(:))+0.02];
            
            for i = 1:length(rSessions)
                sessn = rSessions(i);
                
                % Get sigmoid fit
                f = fit(xvals', means(i,:)', sigmoid, 'start', [.05, 5, 1, 0.5]);
                
                % Get params + param confidence intervals
                coeff_means{i} = coeffvalues(f);
                confint_results{i} = confint(f);
                
                % Add effective midpoint slope as a param
                midpoint_slope_mean = 0.25 * f.a1 * f.b1;  % take dy/dx of A/(1+exp(-b(x-c))), eval at x=c --> dy/dx = 1/4*A*b.
                midpoint_slope_confint = 0.25.*confint_results{i}(:,1).*confint_results{i}(:,2);
                coeff_means{i} = [coeff_means{i} midpoint_slope_mean];
                confint_results{i} = [confint_results{i} midpoint_slope_confint];
                
                if contains(Monkeys(m).Sessions(sessn).ShortName, 'Pre')
                    pre_sessn = i;
                elseif contains(Monkeys(m).Sessions(sessn).ShortName, 'Post')
                    post_sessn = i;
                end
                
                % Plot underlying raw data + error
%                 errorbar(xvals, means(i,:), stds(i,:),...
%                     'DisplayName', Monkeys(m).Sessions(rSessions(i)).ShortName,...
%                     'Color', mlc(i))
                scatter(xvals, means(i,:), 10, mlc(i), 'filled', ...
                    'DisplayName', Monkeys(m).Sessions(rSessions(i)).ShortName)
                
                % Show single excluded value (n= 2 units) for Marta
%                 scatter(excluded_xvals_lower, excluded_means(i,:), 20, mlc(i),...
%                     'HandleVisibility', 'off');
                
                % Plot fit 95% conf intervals
                eb = predint(f, xvals, 0.95, 'functional', 'on');  % 'on' means simultaneously across all params
                patch([xvals fliplr(xvals)], [eb(:,1)', fliplr(eb(:,2)')], mlc(i), 'FaceAlpha', 0.25, 'HandleVisibility', 'off'); 
                
                plot(xvals, f(xvals), '--', 'Color', mlc(i) + 0.1, 'HandleVisibility', 'off')
                plot([f.c1 f.c1], yl, '-', 'Color', mlc(i) + 0.1, 'HandleVisibility', 'off', 'LineWidth', 1.25)
                fprintf('Monkey %s, session %s: amp %g, midpoint %g, b %g, midpoint slope %g, offset %g \n',...
                    Monkeys(m).Name,...
                    Monkeys(m).Sessions(rSessions(i)).ShortName,...
                        round(f.a1,3),...
                        round(10^f.c1,3),...
                        round(f.b1,3),...
                        round(midpoint_slope_mean,3),...
                        round(f.d, 3))
            end
%             axis equal
            set(gca, 'XScale', 'log')
            xlabel('Num best units used')
            ylabel('Accuracy')
            ylim(ylims_by_monk{m})
            xticks(log10([3 10 25 100]))
            xticklabels([3 10 25 100])
            title(sprintf('%s, interval %d - %d', Monkeys(m).Name, interval(1), interval(2)), 'Interpreter', 'none')
%             legend
            formatSVMPlot(gca, gcf, 16)
            
            % Compare pre/post coeff conf ints
            n_coeffs = 5;
            coeff_names = {'Amp', 'Shape', 'Midpoint', 'Intercept', 'Midpt Slope'};
            figure("Position", [700 700 600 150])
            for iCoeff = 1:n_coeffs
                subplot(1,n_coeffs,iCoeff)
                hold on
                pre_int = confint_results{pre_sessn}(:,iCoeff);
                pre_mean = coeff_means{pre_sessn}(iCoeff);
                post_int = confint_results{post_sessn}(:,iCoeff);
                post_mean = coeff_means{post_sessn}(iCoeff);
                errorbar(0, pre_mean, pre_int(1) - pre_mean, pre_int(2) - pre_mean, 'color', mlc(1), 'LineWidth', 2)
                errorbar(1, post_mean, post_int(1) - post_mean, post_int(2) - post_mean, 'color', mlc(2), 'LineWidth', 2)
                title(coeff_names{iCoeff})
                formatSVMPlot(gca, gcf, 16)
                xlim([-0.5, 1.5])
                xticks([])
            end
            sgtitle(sprintf('%s, interval %d - %d', Monkeys(m).Name, interval(1), interval(2)), 'Interpreter', 'none')
            
        end
    end
end

%% (before plotting timecourse) Run cluster-based permutation statistics for a given neural subset (ie num units = 50) (fast) (shuffle vs real, pre vs post)
% See: https://www.nature.com/articles/s41593-018-0148-7, "Cluster-based
% permutation procedure" in methods section.

% This section of code automatically looks at shuffled data -- do not 
% include in list of locs.
% rArrayLocs_stats = {'te'};
% rArrayLocs_stats = {'anterior', 'middle', 'posterior'};
top_n_units = 100;

% Statistics parameters
cluster_alpha = 0.05;


% kfl_base_name = sprintf('KFL_SingleUnitSVMRanking_sparsa_Top_%d', top_n_units);
kfl_base_name = sprintf('KFL_SUnSVMRk_sprs_T_%d_TrNMch', top_n_units);  % for individ arrays

% 1. run ttests at each time point
% 2. Generate clusters and cluster statistics (real and shuffled)
% 3. Compare ranked statistics from shuffle to real to find pvals (mult
% comp correction)

% for m = 1:length(Monkeys)
for m = 2
    rSessions = rSessionsByMonk{m};
    if numel(rSessions) ~= 2
        error('Shuffle permutations code expects two sessions to compare. \n Were you trying to do a line plot of one interval across days? See below.')
    end
    
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};

        % Get list of real tstats
        tStats = getrealtstats(SVM, m, rSessions, loc, manualIntervals, kfl_base_name); % 1 x num intervals
        
        % Get lists of tstats from shuffled data
        try
            shuff_name = "_SHUFFLE";
            tStats_SHUFFLED = getshuffledtstats(SVM, m, rSessions, loc, manualIntervals, strcat(kfl_base_name, shuff_name)); % num shuffle permutations x num intervals
        catch
            shuff_name = "_SHUFF";  % had to do this for some of them due to field name length limits!
            tStats_SHUFFLED = getshuffledtstats(SVM, m, rSessions, loc, manualIntervals, strcat(kfl_base_name, shuff_name)); % num shuffle permutations x num intervals
        end
        
        % Find real clusters
        % clusterStats: 1 x num clusters vector (num clusters calculated in
        % function). Is sorted.
        % intsInCluster: of the I intervals in rIntervals, which ones are
        % in each cluster? Sorted to correspond to clusterStats.
        
        [clusterStats, intsInCluster] = calculateTClusterStats(tStats); 
        
        % Find shuffled clusters
        nShuff = size(tStats_SHUFFLED,1);
        if nShuff == 1
            error('Cannot run permutation tests with only one subset of shuffles!')
        end
        nClust = size(clusterStats,2);
        clusterStats_SHUFFLED = zeros(nClust, nShuff);
        for iShuff = 1:nShuff
            clusterStats_SHUFFLED(:,iShuff) = calcShuffClusterStats(tStats_SHUFFLED(iShuff,:), nClust);
        end
        
        
        % Rank match to create p-values. Ie, ask, how many times is the
        % largest true cluster t-stat greater than the shuffled cluster
        % t-stats? How many times is the *second-largest* real t-stat >
        % *second-largest* shuffled t-stat? Etc...
        signf_clusters = [];
        for iRank = 1:length(clusterStats)
            true_val = clusterStats(iRank);
            shuffled_vals = clusterStats_SHUFFLED(iRank,:);  % 1 x nshuffles
            pval = (sum(true_val <= shuffled_vals) / numel(shuffled_vals));
            fprintf("clust %d, pval %0.4g \n", iRank, pval)
            signf = (pval < (cluster_alpha / (length(clusterStats)*2)));
            if signf
                signf_clusters = [signf_clusters iRank];
            end
        end
        
        % Of the intervals in rIntervals, which are significant?
        signf_ints = [intsInCluster{signf_clusters}];
        sigID = sprintf('ClustPermSignfInds_%s_Sessions_%d_vs_%d_Top%d', loc, rSessions(1), rSessions(2), top_n_units);
        sigID2 = sprintf('ClustPermSignfInds_byClust_%s_Sessions_%d_vs_%d_Top%d', loc, rSessions(1), rSessions(2), top_n_units);
        Monkeys(m).(sigID) = signf_ints;
        tmp = intsInCluster(signf_clusters);
        cluster_sorting = argsort(cellfun(@(v) v(1), tmp));
        Monkeys(m).(sigID2) = tmp(cluster_sorting);
    end
end

%% (before same) Check when timecourse first signf above shuffle
above_chance_alpha = 0.05;
% kfl_base_name = sprintf('KFL_SingleUnitSVMRanking_sparsa_Top_%d', top_n_units);
kfl_base_name = sprintf('KFL_SUnSVMRk_sprs_T_%d_TrNMch', top_n_units);  % for individ arrays

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        session = rSessions(i);
        for iLoc = 1:length(rArrayLocs)
            loc = rArrayLocs{iLoc};

            pvals = zeros(length(manualIntervals),1);
            for iInt = 1:length(manualIntervals)
                interval = manualIntervals{iInt};

                % Get real data
                kflID = get_good_interval_name2(interval, loc, kfl_base_name);
                accs = 1 - SVM(m).Sessions(session).(kflID);
                
                % Get shuffled data
                try
                    shuff_name = "_SHUFFLE";
                    shuffID = get_good_interval_name2(interval, loc, strcat(kfl_base_name, shuff_name));
                    accs_SHUFFLED = 1 - SVM(m).Sessions(session).(shuffID);
                catch
                    shuff_name = "_SHUFF";  % had to do this for some of them due to field name length limits!
                    shuffID = get_good_interval_name2(interval, loc, strcat(kfl_base_name, shuff_name));
                    accs_SHUFFLED = 1 - SVM(m).Sessions(session).(shuffID);
                end
                
                % Run t-test
                [~,pval] = ttest(accs_SHUFFLED(:), mean(accs), 'Tail', 'left');  % test accs > 0.50
                pvals(iInt) = pval;
            end
        end
        
        figure
        plot(pvals < 0.05)
        title(Monkeys(m).Sessions(session).ShortName)
        
        sigID = sprintf('FirstIntAbvChance_%s_Top%d', loc, top_n_units);
        Monkeys(m).Sessions(session).(sigID) = iInt;
    end
end

%% Plot timecourse for a given neural subset (Fig 4A/B)

% Plotting params
plot_alpha = 0.4; % transparency of sem fill
% mkYLims = {[0.45 0.85], [0.45 0.65]};
mkYLims = {[0.47 0.8], [0.47 0.8]};

% rArrayLocs = {'te'};
% rArrayLocs = {'anterior', 'middle', 'posterior'};

top_n_units = 100;  % top 100 for fig 2B
% base_kfl_name = sprintf('KFL_SingleUnitSVMRanking_sparsa_Top_%d', top_n_units);  % top 10, 25, 50, or 100
base_kfl_name = sprintf('KFL_SUnSVMRk_sprs_T_%d_TrNMch', top_n_units);  % for individ arrays

figure2('Position', [400 400 1000 600])
% tiledlayout(length(Monkeys), length(rArrayLocs))
for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    % Pre-allocate vectors for plotting
    accMeans = zeros(length(rSessions), length(manualIntervals), length(rArrayLocs));
    accSems = zeros(length(rSessions), length(manualIntervals), length(rArrayLocs));
    
    % Prepare figure
%     figure2('Position', [400 400 1000 600])
%     figure('Position', [580, 640, 1125, 355])
%     figure
    subplot(2,1,m)
    
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        % New subplot for each array location
%         subplot(2, ceil(length(rArrayLocs)/2), iLoc)
%         subplot(1, 3, iLoc)
%         nexttile
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
            plot(starts, meanVector, '-', 'LineWidth', 1, 'Color', mlc(i),...
                'DisplayName', Monkeys(m).Sessions(sessn).ShortName)
            fill([starts fliplr(starts)],...
                [meanVector + semVector, fliplr(meanVector - semVector)], mlc(i),...
                'FaceAlpha', plot_alpha,'linestyle','none', ...
                'HandleVisibility', 'off');
            
            % Plot signf inds
            sigID = sprintf('ClustPermSignfInds_%s_Sessions_%d_vs_%d_Top%d', loc, rSessions(1), rSessions(2), top_n_units);
            sigID2 = sprintf('ClustPermSignfInds_byClust_%s_Sessions_%d_vs_%d_Top%d', loc, rSessions(1), rSessions(2), top_n_units);
            try 
%                 if ~ isempty(Monkeys(m).(sigID))
%                     plot(starts(Monkeys(m).(sigID)), mkYLims{m}(2)-0.02, 'ko', 'MarkerFaceColor', 'k')
%                 end
                if ~isempty(Monkeys(m).(sigID2))
                    % Do some shenanigans to create well-spaced grayscale
                    % colors for the different t-test clusters
                    nclust = length(Monkeys(m).(sigID2));
                    grays = gray(10);
                    color_idxs = round(linspace(0,10,nclust+1))+1;
                    for iClust = 1:nclust
                        plot(starts(Monkeys(m).(sigID2){iClust}), mkYLims{m}(2)-0.02,...
                            'ko', 'MarkerFaceColor', grays(color_idxs(iClust), :),...
                            'MarkerEdgeColor', 'none')
                    end
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
            yticks([0.5 0.65 0.8])
            xlim([-100 300])
%             if m == 1
%                 yticks([0.5 0.8])
%             end
            
            % Detailed labels if desired
            %             ylabel('SVM abcat accuracy (mean +/- SEM)')
            title(sprintf('%s', loc), 'Interpreter', 'none')
            
            
            % Make the plot look nice
            formatSVMPlot(gca, gcf)
        end
    end
    % More detailed labels if desired.
    %         sgtitle(sprintf('%s', Monkeys(m).Name), 'Interpreter', 'none')

    % Save the plots
%     pause(0.5)
%     saveas(gcf, fullfile(figureSavePath, sprintf('pv_SVM_Timecourse_%s_%s_%g', Monkeys(m).Name, sigID, ID)), 'epsc')
end

%% Plot a given subset/interval over sessions (Fig 4C)

% rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};  % Fig 4C

top_n_units = 100;
interval_to_plot = [175 275];
% base_kfl_name = sprintf('KFL_SingleUnitSVMRanking_sparsa_Top_%d', top_n_units);  % top 10, 25, 50, or 100
base_kfl_name = sprintf('KFL_SingleUnitSVMRanking_sparsa_Top_%d_TrNMatch', top_n_units);
% base_kfl_name = sprintf('KFL_SUnSVMRk_sprs_T_%d_TrNMch', top_n_units);

figure2('Position', [400 400 1000 600])
% mkYLims = {[0.5 0.8], [0.5 0.65]};  % fig 4B, bottom at 0.5
% mkYLims = {[0.65 0.8], [0.57 0.63]};  % fig 4B, tight lims
mkYLims = {[0.5 0.75], [0.5 0.7]};  % car truck

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    % Pre-allocate vectors for plotting
    accMeans = zeros(length(rSessions), length(rArrayLocs));
    accSems = zeros(length(rSessions), length(rArrayLocs));
    
    % Prepare figure
%     figure2('Position', [400 400 1000 600])
    subplot(2,1,m)
    hold on
    
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        % New subplot for each array location
%         subplot(2, ceil(length(rArrayLocs)/2), iLoc)
%         nexttile
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
            
            % Get field names
            kflID = get_good_interval_name2(interval_to_plot, loc, base_kfl_name);

            % Get data
            kfls = SVM(m).Sessions(sessn).(kflID);
            kfls = kfls(:); % reshape into 1 x numel
            accMeans(i, iLoc) = mean(1 - kfls);
            accSems(i, iLoc) = std(1 - kfls) / numel(kfls);
        end
        
        % Plot a line for each array
%         errorbar(accMeans(:, iLoc), accSems(:,iLoc), '-o',...
%             'LineWidth', 2, 'Color', mlc(iLoc),...
%             'DisplayName', loc)

        disp(accMeans(:, iLoc))
        bar(accMeans(:, iLoc))
        errorbar(accMeans(:, iLoc), accSems(:,iLoc), 'o',...
            'LineWidth', 2, 'Color', 'k',...
            'MarkerSize', 0.1)
        

        % Add labels
        xticks(1:length(Monkeys(m).Sessions(rSessions)))
        xticklabels({Monkeys(m).Sessions(rSessions).ShortName})
        xtickangle(90)
        ylabel("Cross-val. SVM accuracy")

        % Make graphs have the same y-axes within each monkey
        ylim(mkYLims{m})

        % Make the plot look nice
        formatSVMPlot(gca, gcf, 16)
    end
    
    % More detailed labels if desired.
%     sgtitle(sprintf('%s', Monkeys(m).Name), 'Interpreter', 'none')

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
rSessionsByMonk = {[7 9], [6 7]};

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
% base_kfl_name = sprintf('KFL_SingleUnitSVMRanking_sparsa_Top_%d', top_n_units);  % top 10, 25, 50, or 100
base_kfl_name = sprintf('KFL_SUnSVMRk_sprs_T_%d_TrNMch', top_n_units);  % for individ arrays

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
            % Plot signf inds
            sigID2 = sprintf('ClustPermSignfInds_byClust_%s_Sessions_%d_vs_%d_Top%d', loc, rSessions(1), rSessions(2), top_n_units);
            try 
%                 if ~ isempty(Monkeys(m).(sigID))
%                     plot(starts(Monkeys(m).(sigID)), mkYLims{m}(2)-0.02, 'ko', 'MarkerFaceColor', 'k')
%                 end
                if ~isempty(Monkeys(m).(sigID2))
                    % Do some shenanigans to create well-spaced grayscale
                    % colors for the different t-test clusters
                    nclust = length(Monkeys(m).(sigID2));
                    grays = gray(10);
                    color_idxs = round(linspace(0,10,nclust+1))+1;
                    for iClust = 1:nclust
                        plot(starts(Monkeys(m).(sigID2){iClust}), mkYLims{m}(2)-0.02,...
                            'ko', 'MarkerFaceColor', grays(color_idxs(iClust), :),...
                            'MarkerEdgeColor', 'none')
                    end
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