% figs for the cat/dog manuscript, 18 Oct 2023

%% Load behavioral data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
svmRecordPath = 'XMA2/Monkey_structs/SVM_Records.mat';
pv_path = 'XMA2/Monkey_structs';
fullSVMPath = fullfile(EXT_HD, pv_path, 'SVM_results_%g.mat');

% Load behavioral data
load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_behav_and_metaNI.mat')) % behavior and neural summaries, but w/o spike times

%% Global values

TE_LOCS = {'anterior', 'middle', 'posterior'};

img_sets = {[1:20 261:280], [21:260 281:520], [1:260, 261:520]}; % {[cats_set1 dogs_set1], [cats_set2...]}
img_set_names = {'train', 'test', 'all'};

% General SVM Parameters
random_seed = 10; % for reproducibility 
ignoreVal = 20; % if neuron has less than this num spikes, do not use it.

manuallyReduceIntervals = true; % test a subset of all the intervals for faster testing
%     manualIntervals = {[75 175], [175 275]}; % NB, also need to change variable fname_base to grab file containing desired intervals
%     manualIntervals = {[175 275]}; % NB, also need to change variable fname_base to grab file containing desired intervals

% pre-computed spike count data. 
% See pv_get_interval_spike_counts if you want to add more intervals.
step = 5; % Interval parameters (for loading the correct spike counts file)
width = 100;
X_fname_base = sprintf('%%s_allNeurons_step%d_wd%d.mat', step, width); % contains full time-course of spike counts (Fig 2b)
spikeCountPath = 'XMA2/Spike_count_mats';

% Note that folds are not generated purely randomly -- we impose
% constraints about not re-using the same image in train vs test ("abstract
% category") and about balancing the number of cats and dogs in each fold
% (so that the shuffle comes out right at 50 %) (except that it still isnt' precisely 50%,
% probably due to repeating of fixation error trials).
kfold_num = 5; % k-fold cross validation. 
catg2_ind1 = 261;  % imgs 1-260 are cats, 261-520 are dogs


%% Fig 2B prep: Params

% Params
ranking_interval = {[175 275]};
svm_train_set_idx = 3; % 1 for 20/20 training set, 2 for 240/240 testing set, 3 for all
svm_test_set_idx = 3;
rArrayLocs = {'te'};
% rArrayLocs = {'anterior', 'middle', 'posterior'};
rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)
% rSessionsByMonk = {1:9, 1:7};
% rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};  % Fig 2C

% for SVM timecourses (Fig 2B)
starts = -100:10:300;
manualIntervals = arrayfun(@(x) [x x+100], starts, 'UniformOutput', false);
    
% Load data
ID = 921995; % ID: 921995 -- full timecourse (step=10 ms), shuff x5, pre/post, [10 25 50 100] unit subsets. (Fig 2B)
load(sprintf(fullSVMPath, ID), 'data');
SVM = data;
clear data

%% Fig 2B prep: cluster-based permutation statistics
% See: https://www.nature.com/articles/s41593-018-0148-7, "Cluster-based
% permutation procedure" in methods section.
% 1. run ttests at each time point
% 2. Generate clusters and cluster statistics (real and shuffled)
% 3. Compare ranked statistics from shuffle to real to find pvals (mult
% comp correction)

% params
top_n_units = 100;
cluster_alpha = 0.001;

kfl_base_name = sprintf('KFL_SingleUnitSVMRanking_sparsa_Top_%d', top_n_units);
for m = 1:length(Monkeys)
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

%% Fig 2B prep: Check when timecourse first signf above shuffle
above_chance_alpha = 0.05;
kfl_base_name = sprintf('KFL_SingleUnitSVMRanking_sparsa_Top_%d', top_n_units);

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

%% Fig 2B: plot timecourse

% Plotting params
plot_alpha = 0.4; % transparency of sem fill
mkYLims = {[0.45 0.85], [0.45 0.65]};


figure2('Position', [400 400 1000 600])
% tiledlayout(length(Monkeys), length(rArrayLocs))
for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    % Pre-allocate vectors for plotting
    accMeans = zeros(length(rSessions), length(manualIntervals), length(rArrayLocs));
    accSems = zeros(length(rSessions), length(manualIntervals), length(rArrayLocs));
    
    % Prepare figure
    subplot(2,1,m)
    
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        hold on
        
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            for iInt = 1:length(manualIntervals)
                interval = manualIntervals{iInt};
                
                % Get field names
                kflID = get_good_interval_name2(interval, loc, kfl_base_name);
                
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
            formatSVMPlot(gca, gcf)
        end
    end
end

%% Fig 2C prep: load data

% Load data
ID = 316514; % all sessions, [10 25 50 100], {[175 275]}, no shuffle.
load(sprintf(fullSVMPath, ID), 'data');
SVM = data;
clear data

%% Fig 2C prep: params
rArrayLocs = {'te'};

%% Fig 2C: plot a given subset/interval over sessions

rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};  % Fig 2C
top_n_units = 100;
interval_to_plot = [175 275];
kfl_base_name = sprintf('KFL_SingleUnitSVMRanking_sparsa_Top_%d', top_n_units);  % top 10, 25, 50, or 100


figure2('Position', [400 400 1000 600])
mkYLims = {[0.675 0.8], [0.57 0.63]};

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    % Pre-allocate vectors for plotting
    accMeans = zeros(length(rSessions), length(rArrayLocs));
    accSems = zeros(length(rSessions), length(rArrayLocs));
    
    % Prepare figure
    subplot(2,1,m)
    hold on
    
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};

        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            % Get field names
            kflID = get_good_interval_name2(interval_to_plot, loc, kfl_base_name);

            % Get data
            kfls = SVM(m).Sessions(sessn).(kflID);
            kfls = kfls(:); % reshape into 1 x numel
            accMeans(i, iLoc) = mean(1 - kfls);
            accSems(i, iLoc) = std(1 - kfls) / numel(kfls);
        end
        
        bar(accMeans(:, iLoc))
        errorbar(accMeans(:, iLoc), accSems(:,iLoc), 'o',...
            'LineWidth', 2, 'Color', 'k',...
            'MarkerSize', 0.1)
        
        % Make graphs have the same y-axes within each monkey
        ylim(mkYLims{m})

        % Make the plot look nice
        formatSVMPlot(gca, gcf, 16)
    end
end

%% Fig 2D prep: load data

% ID: 198063 -- num_best_units_to_use = [1:1:50 52:2:200]; {[175 275]}, pre/post.
% Also has analysis where we add/remove the best units one at a time. (Fig 2D)
% Still weird sga / sparsa switch, for now seemingly nothing to be done
% about it -- just use top 100 neurons in general.
ID = 198063;
load(sprintf(fullSVMPath, ID), 'data');
SVM = data;
clear data

%% Fig 2D prep: params
kfl_base_name = 'KFL_SingleUnitSVMRanking_sparsa_Top_%d';  % 2D
rArrayLocs = {'te'};
num_best_units_to_use = [1:1:50 52:2:200];

manualIntervals = {[175 275]};
rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)

%% Fig 2D prep: plot raw performance wrt adding units

for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        array = rArrayLocs{iLoc};

        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};

            % pre alloc vector to hold mean + std for accuracies
            means = nan(length(rSessions), length(num_best_units_to_use));
            stds = nan(length(rSessions), length(num_best_units_to_use));
            
            
            for i = 1:length(rSessions)
                sessn = rSessions(i);
                disp(sessn)
                
                for iUnitSubset = 1:length(num_best_units_to_use)
                    unit_subset = 1:(num_best_units_to_use(iUnitSubset));
                    kflID = get_good_interval_name2(interval, array, sprintf(kfl_base_name, num_best_units_to_use(iUnitSubset)));
                    if ~ismember(kflID, fields(SVM(m).Sessions(sessn)))
                        warning('No KFL field found -- did you use wrong field name template?')
                        continue
                    end
                    kfls = SVM(m).Sessions(sessn).(kflID); 
                    disp(kfls)

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

%% Fig 2D: curve fitting for SVM while adding units

sigmoid = fittype('a1/(1 + exp(-b1*(x-c1))) + d',...
    'dependent', 'y', 'independent', 'x',...
    'coefficients', {'a1', 'b1', 'c1', 'd'});

lower_n_units_val_by_monk = [3 2];  % to allow a decent sigmoid fit
ylims_by_monk = {[0.63, 0.8], [0.54, 0.63]};

% lower_n_units_val_by_monk = [2 1]; % to see what a bad fit looks like...

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
                    kflID = get_good_interval_name2(interval, array, sprintf(kfl_base_name, num_best_units_to_use(iUnitSubset)));
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
            
            % Get the few data points at begining that
            % don't fit nicely into the sigmoid curve.
            idx_excluded_lower = find((num_best_units_to_use < lower_n_units_val_by_monk(m)));
            excluded_xvals_lower = log10(num_best_units_to_use(idx_excluded_lower));
            excluded_means = all_means(:, idx_excluded_lower);
            
            figure
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
                scatter(xvals, means(i,:), 20, mlc(i), 'filled', ...
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
            set(gca, 'XScale', 'log')
            xlabel('Log10 num best units used')
            ylabel('Accuracy')
            title(sprintf('%s, interval %d - %d', Monkeys(m).Name, interval(1), interval(2)), 'Interpreter', 'none')
%             legend
            ylim(ylims_by_monk{m})
            formatSVMPlot(gca, gcf, 16)
            
            % Compare pre/post coeff conf ints
            n_coeffs = 5;
            coeff_names = {'Amp', 'Shape', 'Midpoint', 'Intercept', 'Midpt Slope'};
            figure
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
            end
            sgtitle(sprintf('%s, interval %d - %d', Monkeys(m).Name, interval(1), interval(2)), 'Interpreter', 'none')
            
        end
    end
end

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