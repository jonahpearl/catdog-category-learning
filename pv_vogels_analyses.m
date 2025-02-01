% repeat some analyses from Vogels 1999

%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';
spikeCountPath = 'XMA2/Spike_count_mats';

% Load behavioral data
load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_behav_and_metaNI.mat')) % behavior and neural summaries, but w/o spike times
TE_LOCS = {'anterior', 'middle', 'posterior'};

% Interval parameters for loading spike count mat
step = 5;
width = 100;

%% Add 'area' label to unitinfo, get short session names

% add Area labels
for m = 1:length(Monkeys)
    for i = 1:length(Monkeys(m).Sessions)
        for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
            switch Monkeys(m).Sessions(i).UnitInfo(j).Location
                case 'posterior'
                    Monkeys(m).Sessions(i).UnitInfo(j).Area = 'te';
                case 'middle'
                    Monkeys(m).Sessions(i).UnitInfo(j).Area = 'te';
                case 'anterior'
                    Monkeys(m).Sessions(i).UnitInfo(j).Area = 'te';
                case 'teo'
                    Monkeys(m).Sessions(i).UnitInfo(j).Area = 'teo';
            end
        end
    end
end

marta_xls = fullfile(EXT_HD, 'RecordingMarta.xlsx');
max_xls = fullfile(EXT_HD, 'RecordingMax.xlsx');
verification_string = 'cats';

for m = 1:length(Monkeys)
    date_strs = {Monkeys(m).Sessions.DateStr};
    if regexp(Monkeys(m).Name, 'Marta\w*')
        short_names = get_short_names(date_strs, marta_xls, sprintf('\\w*%s\\w*', verification_string)); % need double \\ to get single \ in output
    else
        short_names = get_short_names(date_strs, max_xls, sprintf('\\w*%s\\w*', verification_string));
    end
    for i = 1:length(Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).ShortName = short_names{i};
    end
end

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

%% Run Nst analysis (with and without category)

% Params
thresh_percent = 30;
rSessionsByMonk = {[7 9], [6 7]};
% rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};
ignoreVal = 20;
baseline_norm = false;
runShuffle = false;
nShuffles = 100;

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Get data from appropriate storage place
        % X_full is matrix of spike counts of dims (trials) x (units) x
        % (intervals), rIntervals is list of each spike-count interval.
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf('%s_allNeurons_step%d_wd%d.mat', MonkID, step, width);
        [X_full, rIntervals] = load_interval_data(fullfile(EXT_HD, spikeCountPath, fileName)); 
        Y_imgIDs = Monkeys(m).Sessions(sessn).Session_Y_imageID; % list of image ids for each presentation
        img_by_catg = [ones(260,1) 2*ones(260,1)]; % list of categories for each image (260 cats + 260 dogs)
        
        % Manually reduce intervals tested
        warning('Reducing intervals tested...comment out these lines to run all intervals')
        rIntervals_original = rIntervals; % for finding idx in 3rd dim of X_full
        rIntervals = {[75 175], [175 275], [275 375]}; % for controlling loops
        
        % Pre-allocate the storage vectors
        nst_all = zeros(size(X_full,2), length(rIntervals));
        nst_catdog = zeros(size(X_full,2), length(rIntervals), 2);
        
        if runShuffle
            nst_all_SHUFFLE = zeros(size(X_full,2), length(rIntervals), nShuffles);
            nst_catdog_SHUFFLE = zeros(size(X_full,2), length(rIntervals), 2, nShuffles);
        end
                
        for iUnit = 1:size(X_full,2)
            
            % Catch units with too few spikes
            if Monkeys(m).Sessions(sessn).UnitInfo(iUnit).SpikeNum < ignoreVal
                sparseness(iUnit, :) = NaN;
                if runShuffle
                    nst_all_SHUFFLE(iUnit, :, :) = NaN;
                    nst_catdog_SHUFFLE(iUnit, :, :, :) = NaN;
                end
                continue
            end
            
            % Run the analysis
            for iInt = 1:length(rIntervals)
                interval = rIntervals{iInt};
                idx = find(cellfun(@(a) all(a == interval), rIntervals_original));
                spike_counts = X_full(:, iUnit, idx);
                stim_means = splitapply(@mean, spike_counts, Y_imgIDs); % god bless MATLAB becoming more like python
                
                
                % Share response threshold across cats and dogs, because we
                % want to know about overall tuning related to cats vs
                % dogs.
                % (There's also a monkey struct where I did this the wrong
                % way, with a max for dogs and a max for cats
                % "MaxMarta_xma2_Nst_maxesSplit.mat")
                response_threshold = 0.333 * max(stim_means);
                
                % For all stims
                nst_all(iUnit, iInt) = sum(stim_means >= response_threshold);
                
                % For cats
                cat_response_threshold = 0.333 * max(stim_means(img_by_catg==1));
%                 nst_catdog(iUnit, iInt, 1) = sum(stim_means(img_by_catg==1) >= cat_response_threshold);
                nst_catdog(iUnit, iInt, 1) = sum(stim_means(img_by_catg==1) >= response_threshold);
                
                % For dogs
                dog_response_threshold = 0.333 * max(stim_means(img_by_catg==2));
%                 nst_catdog(iUnit, iInt, 2) = sum(stim_means(img_by_catg==2) >= dog_response_threshold);
                nst_catdog(iUnit, iInt, 2) = sum(stim_means(img_by_catg==2) >= response_threshold);
                
                
                % Run with shuffled values if requested
                if runShuffle
                    for iShuff = 1:nShuffles
                        stim_means = splitapply(@mean, spike_counts, Y_imgIDs(randperm(length(Y_imgIDs)))); % god bless MATLAB becoming more like python
                        
                        % For all stims
                        response_threshold = 0.333 * max(stim_means);
                        nst_all_SHUFFLE(iUnit, iInt) = sum(stim_means >= response_threshold);

                        % For cats
                        cat_response_threshold = 0.333 * max(stim_means(img_by_catg==1));
                        nst_catdog_SHUFFLE(iUnit, iInt, 1) = sum(stim_means(img_by_catg==1) >= cat_response_threshold);

                        % For dogs
                        dog_response_threshold = 0.333 * max(stim_means(img_by_catg==2));
                        nst_catdog_SHUFFLE(iUnit, iInt, 2) = sum(stim_means(img_by_catg==2) >= dog_response_threshold);
                    end
                end
            end
            
            % Report progress.
            fprintf('Done with %s unit # %d session %d \n', Monkeys(m).Name, iUnit, sessn)
        end 
        
        % Store data.
        Monkeys(m).Sessions(sessn).Nst = struct('Nst_all', nst_all, 'Nst_catdog', nst_catdog);
        
        if runShuffle
            Monkeys(m).Sessions(sessn).Nst_SHUFFLE = struct('Nst_all_SHUFFLE', nst_all_SHUFFLE, 'Nst_catdog_SHUFFLE', nst_catdog_SHUFFLE);
        end

    end
end

%% Save Nst analyses
% Remove unecessary fields
for m = 1:length(Monkeys)
    try
        Monkeys(m).Sessions = rmfield(Monkeys(m).Sessions, {'UnitInfo', 'CueInfo', 'TrialInfo', 'TimesBT', 'CodesBT', 'Fixn_err_TC'});
    catch
    end
end

save(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_Nst_sharedMax.mat'), 'Monkeys')

%% Run sparseness analysis (pretty fast)

% Calculates sparseness of average stimulus responses (spike counts) per
% unit. That is to say, for each unit, construct a 1 x (num images) mean
% response vector, then calculate that vector's sparseness.

% This follows Vogels 1999, tho NB Vogels uses per-trial baseline-subtracted spike
% counts, these spike counts aren't adjusted at all.
% Close to 0: very sparse. 1: totally uniform.

% Params
% rSessionsByMonk = {[7 9], [6 7]};
rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};
ignoreVal = 20;
baseline_norm = false;
runShuffle = false;
nShuffles = 10;

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Get data from appropriate storage place
        % X_full is matrix of spike counts of dims (trials) x (units) x
        % (intervals), rIntervals is list of each spike-count interval.
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf('%s_allNeurons_step%d_wd%d.mat', MonkID, step, width);
        [X_full, rIntervals] = load_interval_data(fullfile(EXT_HD, spikeCountPath, fileName)); 
        Y_imgIDs = Monkeys(m).Sessions(sessn).Session_Y_imageID; % list of image ids for each presentation
        
        % Manually reduce intervals tested
        warning('Reducing intervals tested...comment out these lines to run all intervals')
        rIntervals_original = rIntervals; % for finding idx in 3rd dim of X_full
        rIntervals = {[75 175], [175 275]}; % for controlling loops
%         rIntervals = {[175 275]}; % for controlling loops
        
        % Pre-allocate the storage vectors
        sparseness = zeros(size(X_full,2), length(rIntervals));
        dog_sparseness = zeros(size(X_full,2), length(rIntervals));
        cat_sparseness = zeros(size(X_full,2), length(rIntervals));
        
        if runShuffle
            sparseness_SHUFFLE = zeros(size(X_full,2), length(rIntervals), nShuffles);
        end
                
        for iUnit = 1:size(X_full,2)
            
            % Catch units with too few spikes
            if Monkeys(m).Sessions(sessn).UnitInfo(iUnit).SpikeNum < ignoreVal
                sparseness(iUnit, :) = NaN;
                if runShuffle
                    sparseness_SHUFFLE(iUnit, :, :) = NaN;
                end
                continue
            end
            
            % Run the analysis
            for iInt = 1:length(rIntervals)
                interval = rIntervals{iInt};
                idx = find(cellfun(@(a) all(a == interval), rIntervals_original));
                spike_counts = X_full(:, iUnit, idx);
                stim_means = splitapply(@mean, spike_counts, Y_imgIDs); % god bless MATLAB becoming more like python
                sparseness(iUnit, iInt) = calculate_sparseness(stim_means);
                
                % Run cat/dog specific sparsenesses
                cat_sparseness(iUnit, iInt) = calculate_sparseness(stim_means(1:260));
                dog_sparseness(iUnit, iInt) = calculate_sparseness(stim_means(261:end));
                
                
                % Run with shuffled values if requested
                if runShuffle
                    for iShuff = 1:nShuffles
                        stim_means = splitapply(@mean, spike_counts, Y_imgIDs(randperm(length(Y_imgIDs)))); % god bless MATLAB becoming more like python
                        sparseness_SHUFFLE(iUnit, iInt, iShuff) = calculate_sparseness(stim_means);
                    end
                end
            end
            
            % Report progress.
            fprintf('Done with %s unit # %d session %d \n', Monkeys(m).Name, iUnit, sessn)
        end 
        
        % Store data.
        Monkeys(m).Sessions(sessn).Unit_sparseness = sparseness;
        Monkeys(m).Sessions(sessn).Unit_sparseness_cat = cat_sparseness;
        Monkeys(m).Sessions(sessn).Unit_sparseness_dog = dog_sparseness;
        
        if runShuffle
            Monkeys(m).Sessions(sessn).Unit_sparseness_SHUFFLE = sparseness_SHUFFLE;
        end

    end
end

%% Save sparseness data

% Remove unecessary fields
for m = 1:length(Monkeys)
    try
        Monkeys(m).Sessions = rmfield(Monkeys(m).Sessions, {'UnitInfo', 'CueInfo', 'TrialInfo', 'TimesBT', 'CodesBT', 'Fixn_err_TC'});
    catch
    end
end

save(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_sparseness.mat'), 'Monkeys')

%% Run mean-responsive sc analysis

% Params
% rSessionsByMonk = {[7 9], [6 7]};
rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};
ignoreVal = 20;
vr_alpha = 0.05;

% Load visual responsiveness data
fname = 'MaxMarta_VR_img_TTest_jun2021.mat'; % changed ttest2 to ttest (because they're paired!)
Data = load(fullfile(EXT_HD, pv_path, fname));
[status, VR_data] = stitch_monkeyStruct_from_parts(Data);
clear Data

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Get data from appropriate storage place
        % X_full is matrix of spike counts of dims (trials) x (units) x
        % (intervals), rIntervals is list of each spike-count interval.
        MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
        fileName = sprintf('%s_allNeurons_step%d_wd%d.mat', MonkID, step, width);
        [X_full, rIntervals] = load_interval_data(fullfile(EXT_HD, spikeCountPath, fileName)); 
        Y_imgIDs = Monkeys(m).Sessions(sessn).Session_Y_imageID; % list of image ids for each presentation
        img_by_catg = [ones(260,1) 2*ones(260,1)]; % list of categories for each image (260 cats + 260 dogs)
        
        % Manually reduce intervals tested
        warning('Reducing intervals tested...comment out these lines to run all intervals')
        rIntervals_original = rIntervals; % for finding idx in 3rd dim of X_full
        rIntervals = {[175 275]}; % for controlling loops
        baseline_intervals = {[-150, -50]};
        
        % Pre-allocate the storage vectors
        mean_resp_cats = cell(size(X_full,2), length(rIntervals)); % mean response to all VR images (variable number) for each unit
        mean_resp_dogs = cell(size(X_full,2), length(rIntervals)); % also look at top 5 resps to contorl for sparse ffects
                
        for iUnit = 1:size(X_full,2)
            
            % Catch units with too few spikes
            if Monkeys(m).Sessions(sessn).UnitInfo(iUnit).SpikeNum < ignoreVal
                continue
            end
            
            % Run the analysis
            for iInt = 1:length(rIntervals)
                interval = rIntervals{iInt};
                baseline_int = baseline_intervals{iInt};
                
                % Get spike counts for this unit, for this interval
                idx = find(cellfun(@(a) all(a == interval), rIntervals_original));
                spike_counts = X_full(:, iUnit, idx);
                
                % Get mean spike count for each image
                stim_means = splitapply(@mean, spike_counts, Y_imgIDs); % god bless MATLAB becoming more like python
                
                % Find images for which unit was vis. responsive
                tid = get_good_interval_name2(interval, 'full', 'VisResp_test_img_TTEST'); % with a t-test
                bid = get_good_interval_name2(baseline_int, '', '');
                vr_id = strcat(tid,bid);
                data_mat = VR_data(m).Sessions(sessn).(vr_id); % stimuli x units
                cat_idx = data_mat(1:260, iUnit) < vr_alpha;
                dog_idx = data_mat(261:520, iUnit) < vr_alpha;
                
                % Store data
                mean_resp_cats{iUnit, iInt} = stim_means(cat_idx);
                mean_resp_dogs{iUnit, iInt} = stim_means(dog_idx);
                
            end
            
            % Report progress.
            fprintf('Done with %s unit # %d session %d \n', Monkeys(m).Name, iUnit, sessn)
        end 
        
        % Store data.
        Monkeys(m).Sessions(sessn).Mean_VR_responses = struct('Cats', mean_resp_cats, 'Dogs', mean_resp_dogs);
    end
end

%% Save mean-responsive data
% Remove unecessary fields
for m = 1:length(Monkeys)
    try
        Monkeys(m).Sessions = rmfield(Monkeys(m).Sessions, {'UnitInfo', 'CueInfo', 'TrialInfo', 'TimesBT', 'CodesBT', 'Fixn_err_TC'});
    catch
    end
end

save(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_mean_VR_responsive.mat'), 'Monkeys')

%% Load sparseness data

data = load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_sparseness.mat'), 'Monkeys');
for m = 1:length(data.Monkeys)
    for i = 1:length(data.Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).Unit_sparseness = data.Monkeys(m).Sessions(i).Unit_sparseness;
        Monkeys(m).Sessions(i).Unit_sparseness_dog = data.Monkeys(m).Sessions(i).Unit_sparseness_dog;
        Monkeys(m).Sessions(i).Unit_sparseness_cat = data.Monkeys(m).Sessions(i).Unit_sparseness_cat;
%         Monkeys(m).Sessions(i).Unit_sparseness_SHUFFLE = data.Monkeys(m).Sessions(i).Unit_sparseness_SHUFFLE;
    end
end
clear data

%% Load Nst data

data = load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_Nst_sharedMax.mat'), 'Monkeys');
for m = 1:length(data.Monkeys)
    for i = 1:length(data.Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).Nst = data.Monkeys(m).Sessions(i).Nst;
    end
end
clear data

%% Load mean-responsive data

data = load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_mean_VR_responsive.mat'), 'Monkeys');
for m = 1:length(data.Monkeys)
    for i = 1:length(data.Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).Mean_VR_responses = data.Monkeys(m).Sessions(i).Mean_VR_responses;
    end
end
clear data

%% Plotting params
% rSessionsByMonk = {[7 9] [6 7]}; % Fig 2!
rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};
% rSessionsByMonk = {[1 6 7 9], [1 5 6 7]};

% Choose arrays. Treat shuffle as a separate loc, will be easier.
% rArrayLocs = {'te', 'SHUFFLE_te'}; 
rArrayLocs = {'te'};

% rArrayLocs = {'te', 'anterior', 'middle', 'posterior', 'SHUFFLE_te', 'SHUFFLE_anterior', 'SHUFFLE_middle', 'SHUFFLE_posterior'};
% rArrayLocs = {'anterior', 'middle', 'posterior', 'SHUFFLE_anterior', 'SHUFFLE_middle', 'SHUFFLE_posterior'}; 
% rArrayLocs = {'te', 'anterior', 'middle', 'posterior'}; 
% rArrayLocs = {'anterior', 'middle', 'posterior'}; 

%% Plot mean-responsive data

% goal: plot cdfs of stim mean sc's, with lines colored /thickness by how
% many stims that neuron responds to.

% Unclear how to quantify this...maybe as sigmoidal or expontential fits,
% and look at the population distbn of the means of those fits? Yikes like
% three layers deep (visual resp + mean stim response + fit distbn of mean
% stim responses)

colors = cbrewer('seq', 'YlOrRd', 10);


for m = 1:length(Monkeys)
figure
tiledlayout(2,2)

rSessions = rSessionsByMonk{m};
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        mvr = Monkeys(m).Sessions(sessn).Mean_VR_responses;
        cats = {mvr.Cats};
        dogs = {mvr.Dogs};
        
        % Determine color scaling
        min_num_resp = 0;
        max_num_resp = max(max(cellfun(@length, cats)), max(cellfun(@length, dogs)));
        cround = @(num) 1+round(9*(num - min_num_resp)/(max_num_resp - min_num_resp));
        
        nexttile
        hold on
        
        for j = 1:length(cats)
            if isempty(cats{j})
                continue
            end
            [f,x] = ecdf(cats{j});
            plot(x,f, 'Color', colors(cround(length(cats{j})),:))
        end
        
        nexttile
        hold on
        for j = 1:length(dogs)
            if isempty(dogs{j})
                continue
            end
            [f,x] = ecdf(dogs{j});
            plot(x,f, 'Color', colors(cround(length(dogs{j})),:))
        end
        
    end
end

%% Plot Nst all (without catg)

ranksum_alpha = 0.05;

close all
figure(1)
figure(2)
% rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};
rSessionsByMonk = {[7 9] [6 7]};

for m = 1:length(Monkeys)

rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            % Get indices of units in the pValues matrix to use in calculating
            % proprotion of units with signf GLMs
            if strcmp(loc, 'te')
                units = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, loc);
            else
                units = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Location}, loc);
            end
            
            
            % Get data, store for stats
            nst = Monkeys(m).Sessions(sessn).Nst.Nst_all(units, 2)/520;
            if regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre.*')
               pre =  nst; 
               color = mlc(1);
               lw = 1.5;
            elseif regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post.*')
                post = nst;
                color = mlc(2);
                lw = 1.5;
            else
                color = [0.4 0.4 0.4];
                lw = 0.75;
            end
            
            % Plot data
            
            % One subplot per session
%             figure(1)
%             subplot(2, length(rSessions), sub2ind([length(rSessions),2], i, m))
%             hold on
%             title(sprintf('%s, session %s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName))
%             ecdf(nst)
%             legend(rArrayLocs)
%             xlim([0 1])
            
            % One subplot per array, comparing sessions
            figure(2)
            subplot(2, length(rArrayLocs), sub2ind([length(rArrayLocs), 2], iLoc, m))
            hold on
            [f,x] = ecdf(nst);
%             plot(x,f, 'DisplayName', Monkeys(m).Sessions(sessn).ShortName, 'Color', color, 'LineWidth', lw)
            plot(x,f, 'Color', color, 'LineWidth', lw)
            title(sprintf('%s, %s', Monkeys(m).Name, loc))
%             legend({'Pre', 'Post'})
            xlim([0 1])
        end
        
        figure(2)
        formatPlot(gca, gcf)
%         legend
        
        % Stats test pre vs post: test for broadening of tuning curves
        [pval, h] = ranksum(pre, post, 'alpha', ranksum_alpha, 'tail', 'left');
        fprintf('%s, %s, Nst pre vs post (ranksum): h = %d, p = %0.2f \n',...
            Monkeys(m).Name, loc, h, pval)
          
    end
end

%% Scatter Nst cat vs Nst dog
% figure(3)
% figure(4)
rSessionsByMonk = {[7 9], [6 7]};
rArrayLocs = {'te'};
nst_iInt_to_plot = 2;  % intervals are [75 175], [175 275], [275 375]
 
figure('Position', [400 400 400 400])
for m = 1:length(Monkeys)

rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            % Get indices of units in the pValues matrix to use in calculating
            % proprotion of units with signf GLMs
            if strcmp(loc, 'te')
                units = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, loc);
            else
                units = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Location}, loc);
            end
            
            if regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre.*')  
               lstyle = '--';
            elseif regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post.*')
                lstyle = '-';
            else
                lstyle = '-';
            end
            
            % Get data
            cats = Monkeys(m).Sessions(sessn).Nst.Nst_catdog(units, nst_iInt_to_plot, 1)/260;
            dogs = Monkeys(m).Sessions(sessn).Nst.Nst_catdog(units, nst_iInt_to_plot, 2)/260;
            overall = Monkeys(m).Sessions(sessn).Nst.Nst_all(units, nst_iInt_to_plot)/520;
            
            
            subplot(2,1,m)
            hold on
            [f,x] = ecdf(cats);
            plot(x,f, 'Color', 'red', 'LineStyle', lstyle)
            [f,x] = ecdf(dogs);
            plot(x,f, 'Color', 'blue','LineStyle', lstyle)
            
            [f,x] = ecdf(overall);
            plot(x,f, 'Color', 'k', 'LineStyle', lstyle)
            
            
            title(sprintf('%s, %s', Monkeys(m).Name, loc));
%             legend({'Pre', 'Post'})
            xlim([0 1]);
            formatPlot(gca, gcf);
            axis square
            
            % Plot data
%             figure(3)
%             subplot(2,2, sub2ind([2,2], i, m))
%             hold on
%             title(sprintf('%s, session %s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName))
%             scatter(dogs, cats, 'filled')
%             legend(rArrayLocs)
%             xlim([0 1])
%             ylim([0 1])
%             
%             figure(4)
%             subplot(2,3, sub2ind([3,2], iLoc, m))
%             hold on
%             scatter(dogs, cats, 'filled')
%             title(sprintf('%s, %s', Monkeys(m).Name, loc))
%             legend({'Pre', 'Post'})
%             xlim([0 1])
%             ylim([0 1])
        end
    end
end

%% Major axis regression on Nst data
%%%% Run major-axis regression to assess slopes, store residuals
% MAR alpha is not dynamic for figures -- if you want to use a different
% MAR alpha (for slope CI's), come back here and re-run the analysis first.

mar_alpha = 0.05;
alpha_string_scale_factor = 100;
interval_idx = 2; % {[75 175], [175 275], [275 375]};
nst_intervals_used = {[75 175], [175 275], [275 375]};
make_plots = false;
rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};

for m = 1:length(Monkeys)

rSessions = rSessionsByMonk{m};
        
    for i = 1:length(rSessions)
        sessn = rSessions(i);

        % Get data
        cats = Monkeys(m).Sessions(sessn).Nst.Nst_catdog(:, interval_idx, 1)/520;
        dogs = Monkeys(m).Sessions(sessn).Nst.Nst_catdog(:, interval_idx, 2)/520;

        % Calculate correlation
        x = corrcoef(cats, dogs);

        % Get angle of best-fit line with major axis regression
        % method. I checked this function against a reference, it's
        % correct.
        [coeffs,coeff_ints,~,~] = maregress(dogs, cats, mar_alpha);
        theta = rad2deg(atan(coeffs(2)));
        theta_bounds = rad2deg(atan(coeff_ints(2,:)));
        intercept = coeffs(1);
        intercept_bounds = coeff_ints(1,:);

        % Print for easy reference
        interval = nst_intervals_used{interval_idx};
        fprintf('%s, %s, all units, interval %d to %d, correlation %0.2f, slope %0.2f (%0.2f to %0.2f)\n', ...
            Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName,...
            interval(1), interval(2),...
            x(1,2), tan(theta/180*pi), tan(theta_bounds/180*pi))

        % Calculate distance from points to regression line:
        % D = abs(ax0+by0+c) / sqrt(a^2+b^2), where point is
        % (x0,y0), and line is ax+by+c = 0.
        % We have y=mx+c form, so 0=coeffs(2)*x-y+coeffs(1) --> a=coeffs(2), b = -1, and c =
        % coeffs(1).
        perpendicular_residuals = abs(coeffs(2)*dogs + -1*cats + coeffs(1)) / sqrt(coeffs(2)^2 + 1);

        % Store data
        propn_id = get_good_interval_name2(interval, 'full', 'NstMAR_CatDogPropns');
        slope_id = get_good_interval_name2(interval, 'full', 'NstMAR_Slope');
        intercept_id = get_good_interval_name2(interval, 'full', 'NstMAR_Intercept');
        resid_id = get_good_interval_name2(interval, 'full', 'NstMAR_Resids');
        rval_id = get_good_interval_name2(interval, 'full', 'NstMAR_RValue');
        Monkeys(m).Sessions(sessn).(propn_id) = [cats; dogs]; % 2 x units
        Monkeys(m).Sessions(sessn).(resid_id) = perpendicular_residuals; % 1 x units
        Monkeys(m).Sessions(sessn).(slope_id) = [theta_bounds(1) theta theta_bounds(2)]; % 1 x 3
        Monkeys(m).Sessions(sessn).(intercept_id) = [intercept_bounds(1) intercept intercept_bounds(2)]; % 1 x 3
        Monkeys(m).Sessions(sessn).(rval_id) = x(1,2); % pearson correlation coeff 
    end
end

%% Plot selected MAR results with scatterplots
xtickcolors = cbrewer('qual', 'Dark2', 3);
epoch_colors = xtickcolors([3 2 1], :);
% interval_to_plot = [175 350];
interval_to_plot = [175 275];
vr_alpha = 0.05;
alpha_string_scale_factor = 100;
area_to_plot = 'te';
color_by_waveform_class = false;
colors = cbrewer('qual', 'Set1', 5);
f1 = figure('Position', [400 400 860 800]);
hold on
f2 = figure('Position', [400 400 400 800]);
% figure('Position', [400 400 1600 500])
hold on

for m = 1:length(Monkeys)
    sessions_to_use = rSessionsByMonk{m};
    diffs = cell(1,length(sessions_to_use));
    resids = cell(1,length(sessions_to_use));
    
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        
        % Get data
        propn_id = get_good_interval_name2(interval, 'full', 'NstMAR_CatDogPropns');
        t = Monkeys(m).Sessions(sessn).(propn_id);
        cat_props = t(1,:);
        dog_props = t(2,:);
        diffs{i} = dog_props - cat_props;
        slope_id = get_good_interval_name2(interval, 'full', 'NstMAR_Slope');
        slope = tan(deg2rad(Monkeys(m).Sessions(sessn).(slope_id)(2)));
        intercept_id = get_good_interval_name2(interval, 'full', 'NstMAR_Intercept');
        intercept = Monkeys(m).Sessions(sessn).(intercept_id)(2);
        xvals = 0:0.01:1;
        resid_id = get_good_interval_name2(interval, 'full', 'NstMAR_Resids');
        perpendicular_resids = Monkeys(m).Sessions(sessn).(resid_id);
        resids{i} = perpendicular_resids;
        
        % ID session name and set some params accordingly.
        switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
            case 'Base'
                col = epoch_colors(1,:);
                thk = 1;
            case 'Pre'
                col = epoch_colors(2,:);
                thk = 2;
                pre_diff_idx = i;
            case 'Post'
                col = epoch_colors(3,:);
                thk = 2;
                post_diff_idx = i;
        end
                
        % Make scatter plot
        set(0, 'CurrentFigure', f1)
        subplot(length(Monkeys),length(sessions_to_use),...
            length(sessions_to_use)*(m-1) + mod(i-1, length(sessions_to_use)) + 1)
        hold on
        scatter(dog_props, cat_props, 'bo')
%         scatterDiagHist(dog_props, cat_props, -0.2:0.025:0.2, 'bo')
        if i == 1 && m == 1
            xlabel('Fraction dogs causing exc.')
            ylabel('Fraction cats causing exc.')
        end
        xlim([0 1])
        ylim([0 1])
        xticks([0 0.5 1])
        yticks([0 0.5 1])
        title(sprintf('%s', Monkeys(m).Sessions(sessn).ShortName), 'FontWeight', 'normal')
%         title(Monkeys(m).XTickLabs{i},'FontWeight', 'normal')
        formatPlot(gca, f1)
        axis square
        
        % Color in example points
%         if i == 1 && m == 1
%             scatter(dog_props(14), cat_props(14), 'go', 'filled')
%             scatter(dog_props(2), cat_props(2), 'go', 'filled')
%         end

        % Plot reference lines
        set(0, 'CurrentFigure', f1)
        plot(xvals, xvals, 'k--', 'LineWidth', 2)
        plot(xvals, intercept + xvals*slope, 'r--', 'LineWidth', 2)
        
        % Make histograms to add to diag of scatters
        set(0, 'CurrentFigure', f2)
%         subplot(length(Monkeys),length(sessions_to_use),...
%             length(sessions_to_use)*(m-1) + mod(i-1, length(sessions_to_use)) + 1)
        subplot(length(Monkeys), 1, m)
        hold on
%         histogram(dog_props - cat_props, 'BinEdges', -0.2:0.025:0.2)
        histogram(perpendicular_resids, 'BinEdges', 0:0.025:0.2)
        ax = gca;
        formatPlot(ax, f2)
        ax.YAxis.Visible = 'off';
        ax.YGrid = 'off';
        ax.XGrid = 'off';
%         ax.XLabel.FontSize = 40;
    end
    
%     % Informative sgtitle if desired
%     sgtitle(sprintf('Area %s, %d to %d, image VR propns,(alpha = %0.2f)',...
%         area_to_plot, interval_to_plot(1), interval_to_plot(2), vr_alpha), 'Interpreter', 'none')

    % Run variance tests on diffs
    [h,pval] = vartest2(diffs{pre_diff_idx}, diffs{post_diff_idx});
    fprintf('DIFFS Pre var: %0.3g, post var: %0.3g. H = %d, p = %d \n', ...
            var(diffs{pre_diff_idx}), var(diffs{post_diff_idx}), h, pval)
        
    % Ditto on residuals
    [h,pval] = vartest2(resids{pre_diff_idx}, resids{post_diff_idx});
    fprintf('RESIDUALS Pre var: %0.3g, post var: %0.3g. H = %d, p = %d \n', ...
            var(resids{pre_diff_idx}), var(resids{post_diff_idx}), h, pval)
end

% saveas(f1, fullfile(figureSavePath, sprintf('VR_scatter_%s', propn_id)), 'epsc')
% saveas(f2, fullfile(figureSavePath, sprintf('VR_hist_%s', propn_id)), 'epsc')

%% Plot sparseness over sessions
rSessionsByMonk = {[7 9] [6 7]};

% rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};
rArrayLocs = {'te'};
ranksum_alpha = 0.05;
show_shuffle = false;  % not actually that helpful

idx_in_sparseness_mat = 1;  % oops, overwrote w just 175-275, but w shuffle
% idx_in_sparseness_mat = 2;  % 1,2,3 are 75-175, 175-275, 275-375

close all
figure(1)
figure(2)

for m = 1:length(Monkeys)

rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            % Get indices of units in the pValues matrix to use in calculating
            % proprotion of units with signf GLMs
            if strcmp(loc, 'te')
                units = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, TE_LOCS);
            else
                units = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, loc);
            end
            
            sparsenesss_indexes = Monkeys(m).Sessions(sessn).Unit_sparseness(units, idx_in_sparseness_mat); 
            if show_shuffle
                sparsenesss_indexes_shuff = Monkeys(m).Sessions(sessn).Unit_sparseness_SHUFFLE(units, idx_in_sparseness_mat, :); 
            end
            if regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre.*')
               pre = sparsenesss_indexes; 
               color = mlc(1);
               lw = 1.5;
            elseif regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post.*')
                post = sparsenesss_indexes;
                color = mlc(2);
                lw = 1.5;
            else
                color = [0.4 0.4 0.4];
                lw = 0.75;
            end
            
            
            % Prep the ecdf plot
            [f, x] = ecdf(sparsenesss_indexes);
            
            % Make the plot
            figure(1)
            subplot(2, length(rArrayLocs), sub2ind([length(rArrayLocs), 2], iLoc, m))
            hold on
            
%             plot(x, f, 'Color', mlc(i), 'DisplayName', Monkeys(m).Sessions(sessn).ShortName)
            plot(x, f, 'Color', color, 'LineWidth', lw)
            title(sprintf('%s, %s', Monkeys(m).Name, loc))
            xlim([0 1])
%             xlabel('"Sparseness index" (inverse sparseness)')
            xlabel('s')
            
            % Show shuffle if requested
            if show_shuffle
                for iShuff=1:nShuffles
                    [f_shuff, x_shuff] = ecdf(sparsenesss_indexes_shuff(:,:,iShuff));
                    plot(x_shuff, f_shuff, 'Color', color, 'LineWidth', lw-1, 'LineStyle', '--')
                end
            end
            
        end
        % Stats test pre vs post: test for broadening of tuning curves (ie
        % greater sparseness index [i know, it's confusing that going up
        % means less sparseness...that's just how vogels defined it]).
        [pval, h] = ranksum(pre, post, 'alpha', ranksum_alpha, 'tail', 'left');
        fprintf('%s, %s, Nst pre vs post (ranksum): h = %d, p = %0.2f \n',...
        Monkeys(m).Name, loc, h, pval)
        formatPlot(gca, gcf)
    end
end

%% Plot sparseness within sessions, cat vs dog
rSessionsByMonk = {[7 9] [6 7]};
rArrayLocs = {'te'};
ranksum_alpha = 0.05;

idx_in_sparseness_mat = 2;  % 1,2 are 75-175, 175-275
% idx_in_sparseness_mat = 1;

close all
figure(1)

for m = 1:length(Monkeys)

rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            % Get indices of units in the pValues matrix to use in calculating
            % proprotion of units with signf GLMs
            if strcmp(loc, 'te')
                units = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, TE_LOCS);
            else
                units = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, loc);
            end
            
            sparsenesss_indexes = Monkeys(m).Sessions(sessn).Unit_sparseness(units, idx_in_sparseness_mat); 
            sparsenesss_indexes_cat = Monkeys(m).Sessions(sessn).Unit_sparseness_cat(units, idx_in_sparseness_mat); 
            sparsenesss_indexes_dog = Monkeys(m).Sessions(sessn).Unit_sparseness_dog(units, idx_in_sparseness_mat); 
            
            if regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre.*')
               pre = sparsenesss_indexes; 
               color = mlc(1);
               lw = 1.5;
               lstyle = '--';
            elseif regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post.*')
                post = sparsenesss_indexes;
                color = mlc(2);
                lw = 1.5;
                lstyle = '-';
            else
                color = [0.4 0.4 0.4];
                lw = 0.75;
            end
            
            
            % Make the plot
            figure(1)
%             subplot(2, length(rSessions), sub2ind([length(rSessions), 2], i, m))
            subplot(2,1,m)
            hold on
            
            % Plot overall popn sparseness
%             [f, x] = ecdf(sparsenesss_indexes);
%             plot(x, f, 'Color', 'k', 'LineWidth', lw, 'DisplayName', 'overall', 'LineStyle', lstyle)
            
            % Category-specific
            [f, x] = ecdf(sparsenesss_indexes_cat);
            plot(x, f, 'Color', 'red', 'LineWidth', lw, 'DisplayName', 'overall', 'LineStyle', lstyle)
            [f, x] = ecdf(sparsenesss_indexes_dog);
            plot(x, f, 'Color', 'blue', 'LineWidth', lw, 'DisplayName', 'overall', 'LineStyle', lstyle)
            
            % Test for catg differences
            [pval, h] = ranksum(sparsenesss_indexes_cat, sparsenesss_indexes_dog,...
                'alpha', ranksum_alpha, 'tail', 'left');  % test dogs less sparse --> dogs > cats --> Y > X --> left tail
            fprintf('%s, %s, %s, sparseness dog > cat (ranksum): h = %d, p = %0.2f \n',...
                Monkeys(m).Name, loc, Monkeys(m).Sessions(sessn).ShortName, h, pval)
            
            [pval, h] = ranksum(sparsenesss_indexes_cat, sparsenesss_indexes_dog,...
                'alpha', ranksum_alpha);  % test two sided
            fprintf('%s, %s, %s, sparseness dog vs cat (ranksum): h = %d, p = %0.2f \n',...
                Monkeys(m).Name, loc, Monkeys(m).Sessions(sessn).ShortName, h, pval)
            
            
            title(sprintf('%s, %s', Monkeys(m).Name, loc),...
                'Interpreter', 'none')
            xlim([0 1])
%             xlabel('"Sparseness index" (inverse sparseness)')
            xlabel('s')
            formatPlot(gca, gcf)
            
        end
        % Stats test pre vs post: test for broadening of tuning curves (ie
        % greater sparseness index [i know, it's confusing that going up
        % means less sparseness...that's just how vogels defined it]).
        [pval, h] = ranksum(pre, post, 'alpha', ranksum_alpha, 'tail', 'left');
        fprintf('%s, %s, overall sparseness pre vs post (ranksum): h = %d, p = %0.2f \n',...
            Monkeys(m).Name, loc, h, pval)
        
        fprintf("\n")
    end
end



%% Functions
function formatPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',20, 'YGrid', 'on', 'XGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end