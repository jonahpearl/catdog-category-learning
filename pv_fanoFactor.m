% goal of this analysis is to look at the fano factor of neural responses
% before and after learning, as a proxy for neural reliability. Perhaps
% reliability increases after learning as another mechanism of improving
% category coding.

%% Load behavioral data
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

% General SVM Parameters
random_seed = 10; % for reproducibility 
rArrayLocs = {'te'};
rSessionsByMonk = {[7 9], [6 7]};
% rSessionsByMonk = {1:9, 1:7};
ignoreVal = 20; % if neuron has less than this num spikes, do not use it.
manuallyReduceIntervals = true; % test a subset of all the intervals for faster testing
%     manualIntervals = {[75 175], [175 275]}; % NB, also need to change variable fname_base to grab file containing desired intervals
    manualIntervals = {[175 275]}; % NB, also need to change variable fname_base to grab file containing desired intervals

   
% Interval parameters (for loading the correct spike counts file)
step = 5;
width = 100;


% Other Parameters
spikeCountPath = 'XMA2/Spike_count_mats';
TE_LOCS = {'anterior', 'middle', 'posterior'};
catg2_ind1 = 261;

% path to spike counts. See pv_get_interval_spike_counts if you want to add
% more intervals.
% fname_base = sprintf('%%s_allNeurons_step%d_wd%d.mat', step, width); % double %% escapes the first %s
X_fname_base = sprintf('%%s_allNeurons_variableBin_1.mat'); % contains 75-175, 175-225, 175-275, and 175-350.

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


%% Calculate fano factor for each unit (fast)

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
        Y = Monkeys(m).Sessions(sessn).Session_Y_catg; % list of categories for each trial (image) in X (either 1 or 2)
        imgs_shown = Monkeys(m).Sessions(sessn).Session_Y_imageID;  % list of imgs for each trial
        
        % For finding intervals in X_full
        rIntervals_original = rIntervals;
        
        for iInt = 1:length(manualIntervals)
            interval = manualIntervals{iInt};

            % Find index in X's 3rd dim for requested interval
            int_idx = find(cellfun(@(a) all(a == interval), rIntervals_original));

            for iUnit = 1:length(Monkeys(m).Sessions(sessn).UnitInfo)
                units_to_use = [iUnit];

                % Subset data to just this unit / interval
                X_subset = X_full(:, units_to_use, int_idx);
                Y_subset = Y;

                ffs = zeros(length(Monkeys(m).Sessions(sessn).CueInfo), 1);
                % Loop through images
                for iImg = 1:length(Monkeys(m).Sessions(sessn).CueInfo)
                    img = Monkeys(m).Sessions(sessn).CueInfo(iImg).CueID;
                    trials = Monkeys(m).Sessions(sessn).CueInfo(img).TrialNum;
                    spike_counts = X_subset(trials, :);
                    ffs(iImg) = var(spike_counts) / mean(spike_counts);  % this simple FF metric tends to be biased for low spike rates, let's see how it looks here.
                end

                % Also calculate an overall ff for this neuron
                ff_overall = var(X_subset) / mean(X_subset);

                % Store data.
                ffID = get_good_interval_name2(interval, '', 'FF_standard_by_img');
                Monkeys(m).Sessions(sessn).UnitInfo(iUnit).(ffID) = ffs;
                
                ffID = get_good_interval_name2(interval, '', 'FF_overall');
                Monkeys(m).Sessions(sessn).UnitInfo(iUnit).(ffID) = ff_overall;
                

                % Report progress.
                fprintf('Done with %s interval # %d, unit %d, session %d \n', Monkeys(m).Name, iInt, iUnit, sessn)
            end
            
        end
    end
end

%% Plot the results

rArrayLocs = {'te'};
use_only_vr_units = false;
vr_alpha = 0.05; 

figure
for m = 1:length(Monkeys)
    rng(random_seed)
    rSessions = rSessionsByMonk{m};
    
    subplot(2,1,m)
    hold on
        
    for iInt = 1:length(manualIntervals)
        interval = manualIntervals{iInt};
        for iLoc = 1:length(rArrayLocs)
            array = rArrayLocs{iLoc};
            
            ff_mats = {};
                
            for i = 1:length(rSessions)
                sessn = rSessions(i);

                % Get data from appropriate storage place
                % X_full is (trials) x (units) x (intervals)
                MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
                fileName = sprintf(X_fname_base, MonkID);
                [X_full, rIntervals] = load_interval_data(fullfile(EXT_HD, spikeCountPath, fileName)); % X is spike counts, rIntervals is list of each interval
                Y = Monkeys(m).Sessions(sessn).Session_Y_catg; % list of categories for each trial (image) in X (either 1 or 2)
                imgs_shown = Monkeys(m).Sessions(sessn).Session_Y_imageID;  % list of imgs for each trial

                % For finding intervals in X_full
                rIntervals_original = rIntervals;
                
                % get various boolean vectors describing units to use
                if use_only_vr_units
                    test_int = test_intervals{p};
                    baseline_int = baseline_intervals{p};
                    tid = get_good_interval_name2(test_int, 'full', 'VisResp_all_TTEST'); % with a t-test
                    bid = get_good_interval_name2(baseline_int, '', '');
                    vr_id = strcat(tid,bid);
                    vr_bool = Monkeys(m).Sessions(sessn).(vr_id) < vr_alpha;
                else
                    vr_bool = ones(1, length(Monkeys(m).Sessions(sessn).UnitInfo));
                end
                
                if strcmp(rArrayLocs{iLoc}, 'te')
                    units_to_use_bool = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, TE_LOCS);
                else
                    units_to_use_bool = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, array);
                end

                % apply the boolean vectors
                units_to_use_list = find(vr_bool & units_to_use_bool);
                
                % retrieve fano factor data
%                 ffID = get_good_interval_name2(interval, '', 'FF_standard_by_img');
%                 ff_mat = [Monkeys(m).Sessions(sessn).UnitInfo(units_to_use_list).(ffID)];
                ffID = get_good_interval_name2(interval, '', 'FF_overall');
                ff_mat = [Monkeys(m).Sessions(sessn).UnitInfo.(ffID)];
                
                ff_mats = [ff_mats ff_mat];
                
                histogram(ff_mat(:), [0:0.1:4], 'FaceColor', mlc(i), 'DisplayName', Monkeys(m).Sessions(sessn).ShortName, 'Normalization', 'probability')
%                 title(sprintf('%s, pooled across units/imgs', Monkeys(m).Name), 'Interpreter', 'none')
                xlabel("Fano factor")
                
                % Draw mean / std on plots
                mu = nanmean(ff_mat(:));
                sigma = nanstd(ff_mat(:));
                yl = ylim;
                marker_yval = yl(2) * 1.02 + 0.005 * (i-1);
                scatter(mu, marker_yval , 100, mlc(i), 'filled', 'v',  'MarkerEdgeColor', 'k')
                plot([mu-sigma mu+sigma], repelem(marker_yval, 2), '-', 'Color', mlc(i), 'LineWidth', 1.5)
                
                xticks([0 1 2 3])
%                 xlim([0 3])
%                 set(gca, 'yscale', 'log')

%                 legend
                pause(0.5)
                formatGLMPlot(gca, gcf)
            end
            [p, h] = ranksum(ff_mats{1}(:), ff_mats{2}(:));
            disp(p)
            nanmedian(ff_mats{1}(:)), nanmedian(ff_mats{2}(:))
        end
    end
end


%% Functions
function formatGLMPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',28, 'YGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end


