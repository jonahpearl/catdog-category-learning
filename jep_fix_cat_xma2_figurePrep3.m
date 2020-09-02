% script for compiling many analyses into figures


%% ******* Neural data summaries *******

%% Load data
clearvars
close all
path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = 'MaxMarta_xma2_temporalActvnSuppn.mat';
Data = load(fullfile(path1, fname));
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Get session Y, spike counts, add 'area' label to unitinfo

% get session Y
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

% spike counts
catgs = {1:260, 261:520};
intervals_to_test = {[-175 0], [175 350], [75 175]};
for m = 1:length(Monkeys)
    sessions_to_use = 1:length(Monkeys(m).Sessions);
    
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        units_to_use = 1:length(Monkeys(m).Sessions(sessn).UnitInfo);
        
        for p = 1:length(intervals_to_test)
            interval_to_test = intervals_to_test{p};
            
            catg1 = Monkeys(m).Sessions(sessn).CueInfo(catgs{1});
            catg2 = Monkeys(m).Sessions(sessn).CueInfo(catgs{2});
            catg1_timeson = vertcat(catg1.Times_on);
            catg2_timeson = vertcat(catg2.Times_on);
            
            % Pre allocate X
            X = zeros(length(catg1_timeson)+length(catg2_timeson), length(units_to_use));
            
            % collect data in this loop
            for j = 1:length(units_to_use)
                
                % collect all spike times
                unit = units_to_use(j);
                allspike_times = sort(Monkeys(m).Sessions(sessn).UnitInfo(unit).Spike_times);
                
                % collect spike counts for catg1
                catg1_sc = zeros(length(catg1_timeson), 1);
                for k = 1:length(catg1_timeson)
                    [lower_idx, upper_idx] = binarySearch_window(allspike_times, catg1_timeson(k) + interval_to_test(1), catg1_timeson(k) + interval_to_test(2)); % allspike_times must be sorted!
                    if lower_idx == -1 || upper_idx == -1 || upper_idx - lower_idx == -1 % ie no spikes in window
                        catg1_sc(k) = 0;
                        continue
                    else
                        catg1_sc(k) = upper_idx - lower_idx + 1;
                    end
                end
                
                % collect spike counts for catg2
                catg2_sc = zeros(length(catg2_timeson), 1);
                for k = 1:length(catg2_timeson)
                    [lower_idx, upper_idx] = binarySearch_window(allspike_times, catg2_timeson(k) + interval_to_test(1), catg2_timeson(k) + interval_to_test(2)); % allspike_times must be sorted!
                    if lower_idx == -1 || upper_idx == -1 || upper_idx - lower_idx == -1 % ie no spikes in window
                        catg2_sc(k) = 0;
                        continue
                    else
                        catg2_sc(k) = upper_idx - lower_idx + 1;
                    end
                end
                % concat into X
                X(1:length(catg1_sc),j) = catg1_sc;
                X(length(catg1_sc)+1:end,j) = catg2_sc;
            end
            % store
            xid = get_good_interval_name2(interval_to_test, 'full', 'X');
            Monkeys(m).Sessions(sessn).(xid) = X;
        end
    end
end

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

%% Params list
% baseline_windows = {[-200 -180] [-100 -80]};
% window_sz = 20;
% window_slide = 10;
% slide_bounds = [-50 450]; % pls ensure slide fits neatly into bounds
% alpha = 0.001;
% window_starts = slide_bounds(1):window_slide:(slide_bounds(2)-window_sz);

%% Plot average spikes
interval_to_use = [175 350];
limit_to_window_results_summary = 1; % -1 to use all, or pick an outcome:
%{
Window results summary key:
0: suppressed
1: nothing
2: exicted
3: suppressed then excited
4: excited then suppressed
%}

figure
hold on

for m = 1:length(Monkeys)
    % select sessions to plot
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
        sessions_to_use = [1 2 7 9]; % skip base04 (outlier) and post01 (low trial count)
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
        sessions_to_use = [1 2 6 7]; % This is base 01, base02, pre and post. skip base03,04,05 and sub03 (low trial count)
    end
    
    % exclude Max TEO units
    subset = 'te';
    
    % pre-allocate
    ecdf_cell = cell(length(sessions_to_use), 2);
    mean_spikes = zeros(1, length(sessions_to_use));
    std_spikes = zeros(1, length(sessions_to_use));
    
    for i = 1:length(sessions_to_use)
        % find specified units, to exclude Max TEO
        sessn = sessions_to_use(i);
        if strcmp(subset, 'te')
            area_bool = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Area}, {'te'});
        else
            error('Unexpected keyword for subset')
        end
        if limit_to_window_results_summary >= 0
            window_results_bool = ([Monkeys(m).Sessions(sessn).UnitInfo.WindowResults_summary] == limit_to_window_results_summary);
        else
            window_results_bool = ones(1, length(Monkeys(m).Sessions(sessn).UnitInfo));
        end
        units_to_use = find(area_bool & window_results_bool);
        
        % get data
        xid = get_good_interval_name2(interval_to_test, 'full', 'X');
        X = Monkeys(m).Sessions(sessn).(xid);
        spikes = reshape(X(:, units_to_use), 1, size(X,1)*length(units_to_use));
        fprintf('%s, session %d, used %d of %d units\n', Monkeys(m).Name, sessn, length(units_to_use), length(Monkeys(m).Sessions(sessn).UnitInfo));
        % get mean / std
        mean_spikes(i) = mean(spikes);
        std_spikes(i) = std(spikes);
        
        % get ecdf and plot
        [f, x] = ecdf(spikes);
        ecdf_cell(i,:) = {x, f};
        subplot(2,2,2*m)
        hold on
        plot(x, f)
    end
    
    % plot mean +/- std
    subplot(2,2,2*m-1)
    errorbar(mean_spikes, std_spikes)
    
    % format ecdf plots
    subplot(2,2,2*m)
    legend
end

%% Set sessions to use and xtick names and colors
xtickcolors = cbrewer('qual', 'Dark2', 3);
xtickcolors = xtickcolors([3 2 1], :);
for m = 1:length(Monkeys)
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
        Monkeys(m).Sessions_to_use = [1 2 3 5 6 7 9];
        Monkeys(m).XTickLabs= {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
        Monkeys(m).Code = 'R';
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
%         Monkeys(m).Sessions_to_use = [1 2 6 7];
        Monkeys(m).Sessions_to_use = 1:7;
        Monkeys(m).XTickLabs = {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
%         Monkeys(m).XTickLabs = {'Base 1', 'Base 2', 'Pre', 'Post'};
        Monkeys(m).Code = 'X';
    end
end

%% Bar plot of unit array / amount across days

figure
hold on

for m = 1:length(Monkeys)
    
    % set up figure
    subplot(2,1,m)
    xlabel('Session')
    ylabel('Num units')
    title(sprintf('Distribution of units by location, %s', Monkeys(m).Name), 'Interpreter', 'none')
    
    % get mk specific info
    sessions_to_use = Monkeys(m).Sessions_to_use;
    xticklabs = Monkeys(m).XTickLabs;
%     all_arrays = {'anterior', 'middle', 'posterior'};
    all_arrays = {'anterior'};
    
    % pre allocate
    data = zeros(length(sessions_to_use), length(all_arrays));
    
    % get data
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        arrays = unique({Monkeys(m).Sessions(sessn).UnitInfo.Location});
        for j = 1:length(arrays)
            idx = find(strcmp(all_arrays, arrays{j}));
            data(i,idx) = sum(strcmp(arrays{j}, {Monkeys(m).Sessions(sessn).UnitInfo.Location}));
        end
    end
    
    % plot
    bar(data, 'stacked')
    xticks(1:length(data))
    xticklabels({})
%     xticklabels(get_xtick_labs_colored(xticklabs, short_names, xtickcolors))
%     xtickangle(45)
%     legend(all_arrays)
    ylim([0 75])
    yticks(0:25:75)
%     ylabel('Num. units')
    set(gca, 'FontSize', 30)
end

%% ******* GLMs *******

%% Load data
clearvars
close all
% path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
path1 = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = 'MaxMarta_GLMs_VR_contrast.mat';
Data = load(fullfile(path1, fname));
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Set sessions to use and xtick names and colors
xtickcolors = cbrewer('qual', 'Dark2', 3);
xtickcolors = xtickcolors([3 2 1], :);
for m = 1:length(Monkeys)
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
%         Monkeys(m).Sessions_to_use = [1 2 3 5 6 7 9];
%         Monkeys(m).XTickLabs= {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
%         Monkeys(m).Sessions_to_use = [1 7 9];
%         Monkeys(m).XTickLabs = {'Base 1', 'Pre', 'Post'};
        Monkeys(m).Sessions_to_use = [7 9];
        Monkeys(m).XTickLabs = {'Pre', 'Post'};
        Monkeys(m).Code = 'R';
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
%         Monkeys(m).Sessions_to_use = [1 2 3 4 5 6 7];
%         Monkeys(m).XTickLabs = {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
%         Monkeys(m).Sessions_to_use = [1 6 7];
%         Monkeys(m).XTickLabs = {'Base 1', 'Pre', 'Post'};
        Monkeys(m).Sessions_to_use = [6 7];
        Monkeys(m).XTickLabs = {'Pre', 'Post'};
        Monkeys(m).Code = 'X';
    end
end

%% Compile GLM and VR data into bools
% intervals_to_use = {[75 175], [175 350]};
% baseline_intervals = {[-150 -50], [-175 0]};
intervals_to_use = {[175 350]};
baseline_intervals = {[-175 0]};


glm_alpha = 0.05;
VR_alpha = 0.05;

for m = 1:length(Monkeys)
    for p = 1:length(intervals_to_use)
        interval = intervals_to_use{p};
        baseline_int = baseline_intervals{p};
        glmid = get_good_interval_name2(interval, 'full', 'GLM_coeffs');
        catg_signf_id = get_good_interval_name2(interval, 'full', 'catg_signf');
        which_catg_id = get_good_interval_name2(interval, 'full', 'which_catg'); % does neuron fire more for cats or dogs? Negative coeff for catg_2 means cats, pos means dogs
        tid = get_good_interval_name2(interval, 'full', 'VisResp_all_test');
        bid = get_good_interval_name2(baseline_int, '', '');
        full_id = strcat(tid,bid);
        
        % get data
%         for i = 1:length(Monkeys(m).Sessions)
        for i = Monkeys(m).Sessions_to_use
            for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
                coef_table = Monkeys(m).Sessions(i).UnitInfo(j).(glmid);
                
                % mark if signf, and if so, whether stronger for cats or
                % dogs
                if isempty(coef_table)
                    Monkeys(m).Sessions(i).UnitInfo(j).(catg_signf_id) = 0;
                    Monkeys(m).Sessions(i).UnitInfo(j).(which_catg_id) = 0;
                elseif coef_table{'catg_2', 'pValue'} < glm_alpha
                    Monkeys(m).Sessions(i).UnitInfo(j).(catg_signf_id) = 1;
                    if coef_table{'catg_2', 'Estimate'} <= 0
                        Monkeys(m).Sessions(i).UnitInfo(j).(which_catg_id) = 1; % cats
                    else
                        Monkeys(m).Sessions(i).UnitInfo(j).(which_catg_id) = 2; % dogs
                    end
                else
                    Monkeys(m).Sessions(i).UnitInfo(j).(catg_signf_id) = 0;
                    Monkeys(m).Sessions(i).UnitInfo(j).(which_catg_id) = 0;
                end
            end
            
            % concat data across units
            signf_vals = [Monkeys(m).Sessions(i).UnitInfo.(catg_signf_id)];
            which_catg_vals = [Monkeys(m).Sessions(i).UnitInfo.(which_catg_id)];
            vr_vals = ([Monkeys(m).Sessions(i).UnitInfo.(full_id)] < VR_alpha);
            vr_vals(isnan(vr_vals)) = 0;
            
            % store data
            sigid = get_good_interval_name2(interval, 'full', 'catg_signf_MAT');
            catid = get_good_interval_name2(interval, 'full', 'which_catg_MAT');
            vrid = get_good_interval_name2(interval, 'full', 'VR_all_ttest2_MAT');
            Monkeys(m).Sessions(i).(sigid) = signf_vals;
            Monkeys(m).Sessions(i).(catid) = which_catg_vals;
            Monkeys(m).Sessions(i).(vrid) = vr_vals;
        end
    end
end

%% (defunct) Ditto but for trial-balanced controls
% 
% % intervals_to_use = {[75 175], [175 350]};
% % baseline_intervals = {[-150 -50], [-175 0]};
% intervals_to_use = {[175 350]};
% baseline_intervals = {[-175 0]};
% 
% 
% glm_alpha = 0.05;
% VR_alpha = 0.05;
% 
% for m = 1:length(Monkeys)
%     for p = 1:length(intervals_to_use)
%         interval = intervals_to_use{p};
%         baseline_int = baseline_intervals{p};
%         glmid_tb = get_good_interval_name2(interval, 'full', 'GLM_coeffs_TRIALBALANCED');
%         
%         for i = Monkeys(m).Sessions_to_use
%             for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
%                 coef_table_tb = Monkeys(m).Sessions(i).UnitInfo(j).(glmid_tb);
%                 
%                 % mark if signf, and if so, whether stronger for cats or
%                 % dogs
%                 if isempty(coef_table_tb)
%                     Monkeys(m).Sessions(i).UnitInfo(j).(catg_signf_id_tb) = 0;
%                     Monkeys(m).Sessions(i).UnitInfo(j).(which_catg_id_tb) = 0;
%                 elseif coef_table_tb{'catg_2', 'pValue'} < glm_alpha
%                     Monkeys(m).Sessions(i).UnitInfo(j).(catg_signf_id_tb) = 1;
%                     if coef_table_tb{'catg_2', 'Estimate'} <= 0
%                         Monkeys(m).Sessions(i).UnitInfo(j).(which_catg_id_tb) = 1; % cats
%                     else
%                         Monkeys(m).Sessions(i).UnitInfo(j).(which_catg_id_tb) = 2; % dogs
%                     end
%                 else
%                     Monkeys(m).Sessions(i).UnitInfo(j).(catg_signf_id_tb) = 0;
%                     Monkeys(m).Sessions(i).UnitInfo(j).(which_catg_id_tb) = 0;
%                 end
%             end
%  
%             % concat data across units and store
%             signf_vals = [Monkeys(m).Sessions(i).UnitInfo.(catg_signf_id_tb)];
%             which_catg_vals = [Monkeys(m).Sessions(i).UnitInfo.(which_catg_id_tb)];
%             sigid_tb = get_good_interval_name2(interval, 'full', 'catg_signf_MAT_TB');
%             catid_tb = get_good_interval_name2(interval, 'full', 'which_catg_MAT_TB');
%             Monkeys(m).Sessions(i).(sigid_tb) = signf_vals;
%             Monkeys(m).Sessions(i).(catid_tb) = which_catg_vals;
%         end
%     end
% end

%% Proportion VR over days
interval = [175 350];
subset = 'te';
vrid = get_good_interval_name2(interval, 'full', 'VR_all_ttest2_MAT');
figure('Position', [200 200 800 600])
for m = 1:length(Monkeys)
    subplot(2,1,m)
    sessions_to_use = Monkeys(m).Sessions_to_use;
    xticklabs = Monkeys(m).XTickLabs;
    short_names = {Monkeys(m).Sessions(sessions_to_use).ShortName};
    
    % pre allocate
    barmat = zeros(1, length(sessions_to_use));
    
    % get data
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        if strcmp(subset, 'te')
            units_to_use = find(ismember({Monkeys(m).Sessions(sessn).UnitInfo.Area}, {'te'}));
        else
            error('Unexpected keyword for subset')
        end
        vr_bools = Monkeys(m).Sessions(sessn).(vrid);
        barmat(i) = sum(vr_bools(units_to_use), 'omitnan') / numel(units_to_use);
    end
    
    % plot
    bar(barmat)
    xticks(1:length(sessions_to_use))
    xticklabels(get_xtick_labs_colored(xticklabs, short_names, xtickcolors))
    xtickangle(45)
    ylabel('Fraction VR / all')
    yticks(0:0.25:0.75)
    ylim([0 0.8])
%     title(sprintf('Monkey %s', Monkeys(m).Code))
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',36, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
end

set(gcf, 'Color', 'white')

%% Plot prop signf (out of all units) for TE for both monks
colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3    ;
    0          0        0    ];

% params
interval_to_plot = [175 350];
% interval_to_plot = [75 175];
overall_chi_sq_baseline_alpha = 0.05;
cat_dog_chi_sq_baseline_alpha = 0.05;
num_comp = 1;
overall_chi_sq_alpha = overall_chi_sq_baseline_alpha / num_comp;
cat_dog_chi_sq_alpha = cat_dog_chi_sq_baseline_alpha / num_comp;

% set field ids
sigid = get_good_interval_name2(interval_to_plot, 'full', 'catg_signf_MAT');
catid = get_good_interval_name2(interval_to_plot, 'full', 'which_catg_MAT');
% sigid = get_good_interval_name2(interval_to_plot, 'full', 'catg_signf_MAT_TB');
% catid = get_good_interval_name2(interval_to_plot, 'full', 'which_catg_MAT_TB');

% prep fig
figure('Position', [200 200 500 600])
hold on

for m = 1:length(Monkeys)
    subplot(2,1,m)
    hold on
    
    sessions_to_use = Monkeys(m).Sessions_to_use;
    xticklabs = Monkeys(m).XTickLabs;
    short_names = {Monkeys(m).Sessions(sessions_to_use).ShortName};
    
    % exclude Max TEO units
    subset = 'te';
    
    % pre allocate
    catg_mat_X = zeros(length(sessions_to_use), 1);
    catg_mat_N = zeros(length(sessions_to_use), 1);
    which_catg_mat_X = zeros(length(sessions_to_use), 2);
    which_catg_prop_tests = zeros(1, length(Monkeys(m).Sessions));
    
    % get data
    baseline_i_vals = [];
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        signf_vals = Monkeys(m).Sessions(sessn).(sigid); % 1 x num units
        which_catg_vals = Monkeys(m).Sessions(sessn).(catid); % 1 x num units
        
        % recall which day is post later on, for stats
        if ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post', 'once'))
            post_i_val = i;
        elseif ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre', 'once'))
            pre_i_val = i;
        elseif ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Base', 'once'))
            baseline_i_vals = [baseline_i_vals i];
        end
        
        % find specified units, to exclude Max TEO
        if strcmp(subset, 'te')
            units_to_use = find(ismember({Monkeys(m).Sessions(sessn).UnitInfo.Area}, {'te'}));
        else
            error('Unexpected keyword for subset')
        end
        signf_vals = signf_vals(units_to_use);
        which_catg_vals = which_catg_vals(units_to_use);
        
        % Store num signf (X) out of total eligible units (N).
        % This type of storage, rather than fraction (X/N) makes chi sq
        % tests easier later on.
        catg_mat_X(i,1) = sum(signf_vals == 1);
        catg_mat_N(i,1) = length(signf_vals);
        which_catg_mat_X(i,1) = sum(which_catg_vals == 1);
        which_catg_mat_X(i,2) = sum(which_catg_vals == 2);
        which_catg_prop_tests(i) = prop_test([sum(which_catg_vals == 1) sum(which_catg_vals == 2)], repelem(length(which_catg_vals),2), false, cat_dog_chi_sq_alpha);
        if length(which_catg_vals) ~= length(signf_vals)
            error('Length of vectors don''t match')
        end
    end
    
    data = catg_mat_X ./ catg_mat_N;
    
    % plot baseline data
    plot(baseline_i_vals, data(baseline_i_vals), '-o', ...
        'MarkerSize', 8, 'MarkerFaceColor', colors(1,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s, %s', Monkeys(m).Name, subset))
    
    % plot pre/post data
    plot([pre_i_val post_i_val], data([pre_i_val post_i_val]), '-o', ...
        'MarkerSize', 8, 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s, %s', Monkeys(m).Name, subset))

    % mark pre/post
%     plot([3.5 3.5], [0 0.7], 'k--', 'LineWidth', 1.5)
    
    % == Cat/dog bar plots ==
%     bar(which_catg_mat_X ./ catg_mat_N, 'grouped')
%     signf_catg_diffs = find(which_catg_prop_tests);
%     yvals = max(which_catg_mat_X(signf_catg_diffs, 1:2) ./ catg_mat_N(signf_catg_diffs, 1),[],2);
%     scatter(signf_catg_diffs, 0.025 + yvals, 'k*')
%     line([signf_catg_diffs-0.25; signf_catg_diffs+0.25], repmat((0.02+yvals)', 2,1), 'Color', 'k')  % mark signf cat/dog differences
    
    % overall stats testing
%     for i = 1:length(sessions_to_use)
%         if i == post_i_val
%             continue
%         elseif prop_test(catg_mat_X([i post_i_val]), catg_mat_N([i post_i_val]), false, overall_chi_sq_alpha)
%             scatter(i*1.02, data(i) * 1.1, 60, '*k', 'HandleVisibility', 'off')
%         end
%     end


    % pre/post stats testing
    if prop_test(catg_mat_X([pre_i_val post_i_val]), catg_mat_N([pre_i_val post_i_val]), false, overall_chi_sq_alpha)
        scatter(mean([pre_i_val post_i_val]), max(data([pre_i_val post_i_val])) * 1.2, 60, '*k', 'HandleVisibility', 'off')
        plot([pre_i_val post_i_val], 1.1*repelem(max(data([pre_i_val post_i_val])),1,2), 'k-', 'HandleVisibility', 'off')
        
        [h, pval, chisq, df] = prop_test(catg_mat_X([pre_i_val post_i_val]), catg_mat_N([pre_i_val post_i_val]), false, overall_chi_sq_alpha);
        fprintf('%s, pre vs post GLM-signf. proportions are signf. different (p = %d, chisq %0.5f, df %d) \n',...
                Monkeys(m).Name, pval, chisq, df)
        
    end
    % format plot
    ylabel(sprintf('Fraction of units significant \n (p < %0.3f)', glm_alpha))
    yticks(0:0.25:0.75)
    ylim([0 0.8])
    xlim([0.5 length(sessions_to_use)+0.5])
    xticks(1:length(sessions_to_use))
    xticklabels(xticklabs)
%     xticklabels(get_xtick_labs_colored(xticklabs, short_names, xtickcolors))
%     xtickangle(45)
%     title(sprintf('Monkey %s', Monkeys(m).Code))
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',18, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
    %     legend({'Combined', 'Cats', 'Dogs'})
end
set(gcf,'color','w');

%% Same as above but trial-balanced
colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3    ;
    0          0        0    ];

% params
interval_to_plot = [175 350];
% interval_to_plot = [75 175];
overall_chi_sq_baseline_alpha = 0.05;
cat_dog_chi_sq_baseline_alpha = 0.05;
num_comp = 1;
overall_chi_sq_alpha = overall_chi_sq_baseline_alpha / num_comp;
cat_dog_chi_sq_alpha = cat_dog_chi_sq_baseline_alpha / num_comp;

% set field ids
boolid = get_good_interval_name2(interval_to_plot, 'full','GLM_TRIALBALANCED_signfBool');

% prep fig
figure('Position', [200 200 500 600])
hold on

for m = 1:length(Monkeys)
    subplot(2,1,m)
    hold on
    
    sessions_to_use = Monkeys(m).Sessions_to_use;
    xticklabs = Monkeys(m).XTickLabs;
    short_names = {Monkeys(m).Sessions(sessions_to_use).ShortName};
    
    % exclude Max TEO units
    subset = 'te';
    
    % pre allocate
    catg_mat_X = zeros(length(sessions_to_use), 1);
    catg_mat_N = zeros(length(sessions_to_use), 1);
    
    % get data
    baseline_i_vals = [];
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        
        % recall which day is post later on, for stats
        if ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post', 'once'))
            post_i_val = i;
        elseif ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre', 'once'))
            pre_i_val = i;
        elseif ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Base', 'once'))
            baseline_i_vals = [baseline_i_vals i];
        end
        
        % find specified units, to exclude Max TEO
        if strcmp(subset, 'te')
            units_to_use = find(ismember({Monkeys(m).Sessions(sessn).UnitInfo.Area}, {'te'}));
        else
            error('Unexpected keyword for subset')
        end
        
        bools = [Monkeys(m).Sessions(sessn).UnitInfo(units_to_use).(boolid)];
        catg_mat_X(i,1) = sum(bools== 1);
        catg_mat_N(i,1) = length(bools);
    end
    
    data = catg_mat_X ./ catg_mat_N;
    
    % plot baseline data
    plot(baseline_i_vals, data(baseline_i_vals), '-o', ...
        'MarkerSize', 8, 'MarkerFaceColor', colors(1,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s, %s', Monkeys(m).Name, subset))
    
    % plot pre/post data
    plot([pre_i_val post_i_val], data([pre_i_val post_i_val]), '-o', ...
        'MarkerSize', 8, 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s, %s', Monkeys(m).Name, subset))
    
    % pre/post stats testing
    if prop_test(catg_mat_X([pre_i_val post_i_val]), catg_mat_N([pre_i_val post_i_val]), false, overall_chi_sq_alpha)
        scatter(mean([pre_i_val post_i_val]), max(data([pre_i_val post_i_val])) * 1.2, 60, '*k', 'HandleVisibility', 'off')
        plot([pre_i_val post_i_val], 1.1*repelem(max(data([pre_i_val post_i_val])),1,2), 'k-', 'HandleVisibility', 'off')
        
        [h, pval, chisq, df] = prop_test(catg_mat_X([pre_i_val post_i_val]), catg_mat_N([pre_i_val post_i_val]), false, overall_chi_sq_alpha);
        fprintf('%s, pre vs post GLM-signf. proportions are signf. different (p = %d, chisq %0.5f, df %d) \n',...
                Monkeys(m).Name, pval, chisq, df)
    end
    % format plot
    ylabel(sprintf('Fraction of units significant \n (p < %0.3f)', glm_alpha))
    yticks(0:0.25:0.75)
    ylim([0 0.8])
    xlim([0.5 length(sessions_to_use)+0.5])
    xticks(1:length(sessions_to_use))
    xticklabels(get_xtick_labs_colored(xticklabs, short_names, xtickcolors))
    xtickangle(45)
%     title(sprintf('Monkey %s', Monkeys(m).Code))
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',18, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
    %     legend({'Combined', 'Cats', 'Dogs'})
end
set(gcf,'color','w');

%% (old) Plot number significant units over days, array breakdown, separate by mk

% Params
interval_to_plot = [175 350];
glm_alpha = 0.05;
chi_sq_alpha = 0.05;
subsets_to_use = {'te', 'posterior', 'middle', 'anterior'}; % put TE first, or scaling won't work.
all_greens = cbrewer('seq', 'Greens', 10);
greens_to_use = all_greens([4 7 10], :);
colors_to_use = [1 1 1; greens_to_use];
scale_by = 'array'; % area or array
overlay_te = false;

% make figure
figure
hold on

% set field ids
sigid = get_good_interval_name2(interval_to_plot, 'full', 'catg_signf_MAT');
catid = get_good_interval_name2(interval_to_plot, 'full', 'which_catg_MAT');


for m = 1:length(Monkeys)
    subplot(2,1,m)
    hold on
    
    % select sessions to plot
%     if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_use = [1 2 7 9]; % skip base04 (outlier) and post01 (low trial count)
%     elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_use = [1 2 6 7]; % This is base 01, base02, pre and post. skip base03,04,05 and sub03 (low trial count)
%     end
    sessions_to_use = Monkeys(m).Sessions_to_use;    
    
    for sub = 1:length(subsets_to_use)
        subset = subsets_to_use{sub};
        
        % pre allocate
        catg_mat_X = zeros(length(sessions_to_use), 1);
        catg_mat_N = zeros(length(sessions_to_use), 1);
        which_catg_mat_X = zeros(length(sessions_to_use), 2);
        
        % get data
        for i = 1:length(sessions_to_use)
            sessn = sessions_to_use(i);
            signf_vals = Monkeys(m).Sessions(sessn).(sigid); % 1 x num units
            which_catg_vals = Monkeys(m).Sessions(sessn).(catid); % 1 x num units
            
            % recall which day is post later on, for stats
            if ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post', 'once'))
                post_i_val = i;
            end
            
            % find specified units
            if strcmp(subset, 'te')
                units_to_use = find(ismember({Monkeys(m).Sessions(sessn).UnitInfo.Area}, {'te'}));
            elseif strcmp(subset, 'teo')
                units_to_use = find(ismember({Monkeys(m).Sessions(sessn).UnitInfo.Area}, {'teo'}));
                if isempty(units_to_use)
                    fprintf('Skipping session %d monkey %s bc no %s units \n', i, Monkeys(m).Name, units_to_count)
                    continue
                end
            elseif ismember(subset, {'anterior', 'middle', 'posterior'})
                units_to_use = find(ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, subset));
            else
                error('Unexpected keyword for subset')
            end
            
            % filter data
            signf_vals = signf_vals(units_to_use);
            which_catg_vals = which_catg_vals(units_to_use);
            
            % Store num signf (X) out of total eligible units (N).
            % This type of storage, rather than percent (X/N) makes chi sq
            % tests easier later on.
            catg_mat_X(i,1) = sum(signf_vals == 1);
            catg_mat_N(i,1) = length(signf_vals);
            which_catg_mat_X(i,1) = sum(which_catg_vals == 1);
            which_catg_mat_X(i,2) = sum(which_catg_vals == 2);
            if length(which_catg_vals) ~= length(signf_vals)
                error('Length of vectors don''t match')
            end
        end
        
        % format plot
        %         if strcmp(subset, 'te')
        %             wd = 2;
        %             te_N = catg_mat_N;
        %         else
        %             wd = 1;
        %         end
        wd = 2;
        
        % plot
        if ~overlay_te && strcmp(subset, 'te')
            continue
        elseif strcmp(scale_by, 'array')
            plot(catg_mat_X ./ catg_mat_N, '-o', 'LineWidth', wd, ...
                'Color', colors_to_use(sub, :), 'DisplayName', subset)
        elseif strcmp(scale_by, 'area')
            plot(catg_mat_X ./ te_N, '-o', 'LineWidth', wd,...
                'Color', colors_to_use(sub, :), 'DisplayName', subset)
        else
            error('Unexpected keyword for scale_by')
        end
        
        % stats testing
        for i = 1:length(sessions_to_use)
            if i == post_i_val
                continue
            elseif prop_test(catg_mat_X([i post_i_val]), catg_mat_N([i post_i_val]), false, chi_sq_alpha)
                if strcmp(scale_by, 'array')
                    scatter(i*1.02, catg_mat_X(i) / catg_mat_N(i) * 1.04, 'vk', 'LineWidth', wd, 'HandleVisibility', 'off')
                elseif strcmp(scale_by, 'area')
                    scatter(i*1.02, catg_mat_X(i) / te_N(i) * 1.04, 'vk', 'LineWidth', wd, 'HandleVisibility', 'off')
                else
                    error('Unexpected keyword for scale_by')
                end
            end
        end
    end
    
    % format subplot
    %     title(sprintf('%s', Monkeys(m).Name), 'Interpreter', 'none')
    ylim([0 0.9])
    xlabel('Session')
    xticks(1:length(catg_mat_X))
    xticklabels({Monkeys(m).Sessions(sessions_to_use).ShortName})
    
    %     ylabel('Fraction units significant in GLM')
    yticks(0:0.2:0.8)
    
    xtickangle(45)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], ...
        'XMinorTick', 'off', 'YMinorTick', 'off', 'XGrid', 'on', 'YGrid', 'on',...
        'LineWidth', 1,  'fontsize',15, ...
        'fontname', 'Helvetica','FontWeight','Bold', 'LineWidth', 2,...
        'XColor', 'black', 'YColor', 'black')
    legend('Location', 'eastoutside')
end
% format figure
set(gcf,'color','w');

%% (old) Break down by cat/dog, TE or TEO

% Params
interval_to_plot = [175 350];
subsets_to_plot = {'te'};

% make figure
figure
hold on

% set field ids
sigid = get_good_interval_name2(interval_to_plot, 'full', 'catg_signf_MAT');
catid = get_good_interval_name2(interval_to_plot, 'full', 'which_catg_MAT');

for m = 1:length(Monkeys)
    subplot(1,2,m)
    hold on
    
    % select sessions to plot
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
        sessions_to_use = [1 2 7 9]; % skip base04 (outlier) and post01 (low trial count)
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
        sessions_to_use = [1 2 6 7]; % This is base 01, base02, pre and post. skip base03,04,05 and sub03 (low trial count)
    end
    
    for sub = 1:length(subsets_to_plot)
        subset = subsets_to_plot{sub};
        
        % pre allocate
        which_catg_mat_X = zeros(length(sessions_to_use), 2);
        which_catg_mat_N = zeros(length(sessions_to_use), 1);
        
        % get data
        for i = 1:length(sessions_to_use)
            sessn = sessions_to_use(i);
            which_catg_vals = Monkeys(m).Sessions(sessn).(catid); % 1 x num units
            
            % find specified units
            if strcmp(subset, 'te')
                units_to_use = find(ismember({Monkeys(m).Sessions(sessn).UnitInfo.Area}, {'te'}));
            elseif strcmp(subset, 'teo')
                units_to_use = find(ismember({Monkeys(m).Sessions(sessn).UnitInfo.Area}, {'teo'}));
                if isempty(units_to_use)
                    fprintf('Skipping session %d monkey %s bc no %s units \n', i, Monkeys(m).Name, units_to_count)
                    continue
                end
            end
            
            % filter data
            which_catg_vals = which_catg_vals(units_to_use);
            which_catg_mat_N(i,1) = length(which_catg_vals);
            which_catg_mat_X(i,1) = sum(which_catg_vals == 1);
            which_catg_mat_X(i,2) = sum(which_catg_vals == 2);
        end
        
        % plot
        bar(which_catg_mat_X ./ which_catg_mat_N, 'grouped')
        legend({'Cats', 'Dogs'})
    end
    
    % format subplot
    title(sprintf('%s', Monkeys(m).Name), 'Interpreter', 'none')
    ylim([0 0.9])
    xlabel('Session')
    ylabel('Proportion units signf. in GLM')
    xticks(1:length(catg_mat))
    xticklabels({Monkeys(m).Sessions(sessions_to_use).ShortName})
    xtickangle(45)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], ...
        'XMinorTick', 'off', 'YMinorTick', 'on', 'XGrid', 'on', 'YGrid', 'on',...
        'LineWidth', 1,  'fontsize',15, ...
        'fontname', 'Helvetica','FontWeight','Bold', 'LineWidth', 2,...
        'XColor', 'black', 'YColor', 'black')
    %     legend
end
% format figure
set(gcf,'color','w');

%% (old) Break down cat/dogs even further, by array

% Params
interval_to_plot = [175 350];
subsets_to_plot = {'te', 'posterior', 'middle', 'anterior'}; % put TE first, or scaling won't work.
scale_by = 'area'; % area (normalize to num units in TE for all subsets) or array (normalize each subset to itself)

% make figure
figure('Position', [200 200 650 1200])
hold on

% set field ids
sigid = get_good_interval_name2(interval_to_plot, 'full', 'catg_signf_MAT');
catid = get_good_interval_name2(interval_to_plot, 'full', 'which_catg_MAT');

for sub = 1:length(subsets_to_plot)
    subset = subsets_to_plot{sub};
    
    for m = 1:length(Monkeys)
        % select sessions to plot
        if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
            sessions_to_use = [1 2 7 9]; % skip base04 (outlier) and post01 (low trial count)
        elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
            sessions_to_use = [1 2 6 7]; % This is base 01, base02, pre and post. skip base03,04,05 and sub03 (low trial count)
        end
        
        subplot(length(subsets_to_plot), length(Monkeys), ...
            length(Monkeys)*(sub-1) + mod(m-1, length(Monkeys)) + 1)
        hold on
        
        % pre allocate
        which_catg_mat_X = zeros(length(sessions_to_use), 2);
        which_catg_mat_N = zeros(length(sessions_to_use), 1);
        
        % get data
        for i = 1:length(sessions_to_use)
            sessn = sessions_to_use(i);
            which_catg_vals = Monkeys(m).Sessions(sessn).(catid); % 1 x num units
            
            % find specified units
            if strcmp(subset, 'te')
                units_to_use = find(ismember({Monkeys(m).Sessions(sessn).UnitInfo.Area}, {'te'}));
            elseif strcmp(subset, 'teo')
                units_to_use = find(ismember({Monkeys(m).Sessions(sessn).UnitInfo.Area}, {'teo'}));
                if isempty(units_to_use)
                    fprintf('Skipping session %d monkey %s bc no %s units \n', i, Monkeys(m).Name, units_to_count)
                    continue
                end
            elseif ismember(subset, {'anterior', 'middle', 'posterior'})
                units_to_use = find(ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, subset));
            end
            
            % filter data
            which_catg_vals = which_catg_vals(units_to_use);
            which_catg_mat_N(i,1) = length(which_catg_vals);
            which_catg_mat_X(i,1) = sum(which_catg_vals == 1);
            which_catg_mat_X(i,2) = sum(which_catg_vals == 2);
        end
        
        % store TE data, in case we need to normalize by it
        if strcmp(subset, 'te')
            te_N = which_catg_mat_N;
        end
        
        % plot
        if strcmp(scale_by, 'array')
            bar(which_catg_mat_X ./ which_catg_mat_N, 'grouped')
        elseif strcmp(scale_by, 'area')
            bar(which_catg_mat_X ./ te_N, 'grouped')
        end
        
        % format subplot
        title(sprintf('%s', subset))
        ylim([0 0.65])
        if sub == length(subsets_to_plot) && m == 1
            ylabel(sprintf('Proportion units \n signf. in GLM \n (normalized by %s)',scale_by))
            xlabel(sprintf('Session, %s', Monkeys(m).Name), 'Interpreter', 'none')
            xticks(1:length(catg_mat))
            xticklabels({Monkeys(m).Sessions(sessions_to_use).ShortName})
            xtickangle(45)
        elseif sub == length(subsets_to_plot)
            xlabel(sprintf('Session, %s', Monkeys(m).Name), 'Interpreter', 'none')
            xticks(1:length(catg_mat))
            xticklabels({Monkeys(m).Sessions(sessions_to_use).ShortName})
            xtickangle(45)
        else
            xticklabels({})
        end
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], ...
            'XMinorTick', 'off', 'YMinorTick', 'on', 'XGrid', 'on', 'YGrid', 'on',...
            'LineWidth', 1,  'fontsize',15, ...
            'fontname', 'Helvetica','FontWeight','Bold', 'LineWidth', 2,...
            'XColor', 'black', 'YColor', 'black')
        %         legend
    end
end
% format figure
set(gcf,'color','w');

%% ******* SVMs *******

%% Load and plot: Supp Fig 3, trial-fixed SVMs (load)
% load data
clearvars
close all
fname = 'MaxMarta_fix_cat_xma2_SVM4_full_TrialSubset.mat'; % 75-175 and 175-350, arrays, fixed trial nums across sessions
path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
Data = load(fullfile(path1, fname));
[status, SVM] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Load and plot: Supp Fig 3, trial-fixed SVMs (plot)

% overall params
xtickcolors = cbrewer('qual', 'Dark2', 3);
xtickcolors = xtickcolors([3 2 1], :);
for m = 1:length(SVM)
    if strcmp(SVM(m).Name, 'Marta_fix_cat_xma2')
        SVM(m).Sessions_to_use = [1 2 3 5 6 7 9];
        SVM(m).XTickLabs= {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
        SVM(m).Code = 'R';
    elseif strcmp(SVM(m).Name, 'Max_fix_cat_xma2')
%         SVM(m).Sessions_to_use = [1 2 3 4 5 6 7];
%         SVM(m).XTickLabs = {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
        SVM(m).Sessions_to_use = [1 2 6 7];
        SVM(m).XTickLabs = {'Base 1', 'Base 2', 'Pre', 'Post'};
        SVM(m).Code = 'X';
    end
end

% plotting params
marta_subsets5 = {'2500TrialsSubset_SC_filtered', '2500TrialsSubset_posterior_SC_filtered',...
    '2500TrialsSubset_middle_SC_filtered', '2500TrialsSubset_anterior_SC_filtered'};
max_subsets5 = {'2500TrialsSubset_No_TEO_SC_filtered', '2500TrialsSubset_posterior_SC_filtered',...
    '2500TrialsSubset_middle_SC_filtered', '2500TrialsSubset_anterior_SC_filtered'};
colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3;
    0          0        0    ];
interval_to_test = [175 350];
include_shuffle = false;
overlay_te = true;
baseline_ttest_alpha = 0.05;
num_comps = 18;
ttest_alpha = baseline_ttest_alpha / num_comps;
neural_popn_subsets_monk = {marta_subsets5, max_subsets5};
neural_popn_subset_names_monk = {marta_subsets5, max_subsets5};



% prep fig
figure('Position', [200 200 500 600]) % small for talks
% figure('Position', [200 200 1000 1200]) % large for fig
hold on

for m = 1:length(SVM)
    subplot(2,1,m)
    hold on
    
    sessions_to_use = SVM(m).Sessions_to_use;
    xticklabs = SVM(m).XTickLabs;
    short_names = {SVM(m).Sessions(sessions_to_use).ShortName};
    
    % select neural subsets to plot
    neural_popn_subsets = neural_popn_subsets_monk{m};
    neural_popn_subset_names = neural_popn_subset_names_monk{m};
    %     neural_popn_subset_colors = neural_popn_subset_colors_monk{m};
    
     
    % get baseline/pre/post session inds
    baseline_i_vals = [];
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        if ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Post', 'once'))
            post_i_val = i;
        elseif ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Pre', 'once'))
            pre_i_val = i;
        elseif ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Base', 'once'))
            baseline_i_vals = [baseline_i_vals i];
        end
    end
        
    % plot each subset, and run stats testing
    for sub = 1:length(neural_popn_subsets)
        neuron_subset = neural_popn_subsets{sub};
        neuron_subset_name = neural_popn_subset_names{sub};
        
        if ~overlay_te && strcmp(neuron_subset_name, 'All TE')
            continue
        end
        
        % get data and plot, all sessions to plot at once
        kfl_id = get_good_interval_name2(interval_to_test, neuron_subset, 'KFL_abcat');
        kfl_shuffled_id = get_good_interval_name2(interval_to_test, neuron_subset, 'KFL_abcat_SHUFFLED');
        acc_mat = 1 - [SVM(m).Sessions(sessions_to_use).(kfl_id)]; % each kfl_vec is a column vector, can horiz concat across days. plot accuracy instead of loss
        if include_shuffle
            shuffled_acc_mat = 1 - [SVM(m).Sessions(sessions_to_use).(kfl_shuffled_id)];
        end

        % plot the data
        wd = 1;
        col = colors(sub, :);
        eb1 = errorbar(baseline_i_vals, mean(acc_mat(:,baseline_i_vals)),...
            std(acc_mat(:,baseline_i_vals)), ...
            'Color', col, 'LineWidth', 2,...
            'DisplayName', neural_popn_subset_names{sub}); % takes mean across iterations
        eb2 = errorbar([pre_i_val post_i_val], mean(acc_mat(:,[pre_i_val post_i_val])),...
            std(acc_mat(:,[pre_i_val post_i_val])),...
            'Color', col, 'LineWidth', 2,...
            'DisplayName', neural_popn_subset_names{sub}); 
        
        % add shuffled data if desired
        if include_shuffle
            ebs = errorbar(mean(shuffled_acc_mat), std(shuffled_acc_mat), '--', 'LineWidth', 0.5,...
                'DisplayName', neural_popn_subset_names{sub}, 'Color', col); % takes mean across iterations
            ebs.Bar.LineStyle = 'dotted';
            uistack(ebs, 'bottom')
        end
        
        % bring the TE line to the top, so it can be seen clearly
        if overlay_te && strcmp(neuron_subset_name, 'All TE')
            uistack(eb1, 'top')
            uistack(eb2, 'top')
        elseif overlay_te
            uistack(eb1, 'bottom')
            uistack(eb2, 'bottom')
        end
        
        % add pre/post stats testing
        pre = acc_mat(:, pre_i_val);
        post = acc_mat(:, post_i_val);
        if ttest2(pre, post, 'tail', 'left', 'Alpha', ttest_alpha)
            max_y = max([mean(pre) mean(post)]);
            scatter(mean([pre_i_val post_i_val]), max_y*1.06, '*', ...
                'MarkerEdgeColor', col, 'MarkerFaceColor', col,  'HandleVisibility', 'off')
            plot([pre_i_val post_i_val], 1.04*repelem(max_y,1,2), 'k-', 'HandleVisibility', 'off')
        end
    end
    
    % format subplot
    ylabel('Accuracy')
    xlim([0.5 length(sessions_to_use)+0.5])
    xticks(1:length(sessions_to_use))
    xticklabels(get_xtick_labs_colored(xticklabs, short_names, xtickcolors))
    xtickangle(45)
%     title(sprintf('Monkey %s', SVM(m).Code), 'Interpreter', 'none')
%     make_custom_patch_legend(colors(1:4,:), {'All', 'Posterior', 'Middle', 'Anterior'}, 'FontSize', 36)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',28, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
end
set(gcf,'color','w');

%% Load and plot: Fig 2 (TE) and Supp Fig 3C (arrays) (load)
clearvars
close all
fname = 'MaxMarta_fix_cat_xma2_SVM4_full.mat'; % 75-175 and 175-350.
path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
Data = load(fullfile(path1, fname));
[status, SVM] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Load and plot: Fig 2 (TE) and Supp Fig 3C (plot)
% overall params
xtickcolors = cbrewer('qual', 'Dark2', 3);
xtickcolors = xtickcolors([3 2 1], :);
for m = 1:length(SVM)
    if strcmp(SVM(m).Name, 'Marta_fix_cat_xma2')
        SVM(m).Sessions_to_use = [1 2 3 5 6 7 9];
        SVM(m).XTickLabs= {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
        SVM(m).Code = 'R';
    elseif strcmp(SVM(m).Name, 'Max_fix_cat_xma2')
        SVM(m).Sessions_to_use = [1 2 3 4 5 6 7];
        SVM(m).XTickLabs = {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
%         SVM(m).Sessions_to_use = [1 2 6 7];
%         SVM(m).XTickLabs = {'Base 1', 'Base 2', 'Pre', 'Post'};
        SVM(m).Code = 'X';
    end
end

% plotting params
marta_te_data = {'SC_filtered'};
marta_te_name = {'All TE'};
marta_subsets1 = {'SC_filtered', 'posterior_SC_filtered', 'middle_SC_filtered', 'anterior_SC_filtered'};
marta_subsets1_names = {'All TE', 'Posterior', 'Middle', 'Anterior'};
max_te_data = {'No_TEO_SC_filtered'};
max_te_name = {'All TE'};
max_subsets1 = {'No_TEO_SC_filtered', 'posterior_SC_filtered', 'middle_SC_filtered', 'anterior_SC_filtered'}; % arrays
max_subsets1_names = {'All TE', 'Posterior', 'Middle', 'Anterior'};
colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3;
    0          0        0    ];
interval_to_test = [175 350];
include_shuffle = false;
overlay_te = false;
baseline_ttest_alpha = 0.05;
num_comps = 18;
ttest_alpha = baseline_ttest_alpha / num_comps;
neural_popn_subsets_monk = {marta_subsets1, max_subsets1};
neural_popn_subset_names_monk = {marta_subsets1_names, max_subsets1_names};
% neural_popn_subsets_monk = {marta_te_data, max_te_data};
% neural_popn_subset_names_monk = {marta_te_name, max_te_name};


% prep fig
figure('Position', [200 200 500 600])
hold on

for m = 1:length(SVM)
    subplot(2,1,m)
    hold on
    sessions_to_use = SVM(m).Sessions_to_use;
    xticklabs = SVM(m).XTickLabs;
    short_names = {SVM(m).Sessions(sessions_to_use).ShortName};
    
    % select neural subsets to plot
    neural_popn_subsets = neural_popn_subsets_monk{m};
    neural_popn_subset_names = neural_popn_subset_names_monk{m};
    %     neural_popn_subset_colors = neural_popn_subset_colors_monk{m};
    
     
    % get baseline/pre/post session inds
    baseline_i_vals = [];
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        if ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Post', 'once'))
            post_i_val = i;
        elseif ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Pre', 'once'))
            pre_i_val = i;
        elseif ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Base', 'once'))
            baseline_i_vals = [baseline_i_vals i];
        end
    end
        
    % plot each subset, and run stats testing
    for sub = 1:length(neural_popn_subsets)
        neuron_subset = neural_popn_subsets{sub};
        neuron_subset_name = neural_popn_subset_names{sub};
        
        if ~overlay_te && strcmp(neuron_subset_name, 'All TE')
            continue
        end
        
        % get data and plot, all sessions to plot at once
        kfl_id = get_good_interval_name2(interval_to_test, neuron_subset, 'KFL_abcat');
        kfl_shuffled_id = get_good_interval_name2(interval_to_test, neuron_subset, 'KFL_abcat_SHUFFLED');
        acc_mat = 1 - [SVM(m).Sessions(sessions_to_use).(kfl_id)]; % each kfl_vec is a column vector, can horiz concat across days. plot accuracy instead of loss
        if include_shuffle
            shuffled_acc_mat = 1 - [SVM(m).Sessions(sessions_to_use).(kfl_shuffled_id)];
        end

        % plot the data
        wd = 1;
        col = colors(sub, :);
%         col = neural_popn_subset_colors(sub,:);
%         eb = errorbar(mean(acc_mat), std(acc_mat), 'LineWidth', wd,...
%             'DisplayName', neural_popn_subset_names{sub}, 'Color', col); % takes mean across iterations
        eb1 = errorbar(baseline_i_vals, mean(acc_mat(:,baseline_i_vals)),...
            std(acc_mat(:,baseline_i_vals)), ...
            'Color', col, 'LineWidth', 2,...
            'DisplayName', neural_popn_subset_names{sub}); % takes mean across iterations
        eb2 = errorbar([pre_i_val post_i_val], mean(acc_mat(:,[pre_i_val post_i_val])),...
            std(acc_mat(:,[pre_i_val post_i_val])),...
            'Color', col, 'LineWidth', 2,...
            'DisplayName', neural_popn_subset_names{sub}); 
        
        % add shuffled data if desired
        if include_shuffle
            ebs = errorbar(mean(shuffled_acc_mat), std(shuffled_acc_mat), '--', 'LineWidth', 0.5,...
                'DisplayName', neural_popn_subset_names{sub}, 'Color', col); % takes mean across iterations
            ebs.Bar.LineStyle = 'dotted';
            uistack(ebs, 'bottom')
        end
        
        % bring the TE line to the top, so it can be seen clearly
        if overlay_te && strcmp(neuron_subset_name, 'All TE')
            uistack(eb1, 'top')
            uistack(eb2, 'top')
        elseif overlay_te
            uistack(eb1, 'bottom')
            uistack(eb2, 'bottom')
        end
        
        % add pre/post stats testing
        pre = acc_mat(:, pre_i_val);
        post = acc_mat(:, post_i_val);
        if ttest2(pre, post, 'tail', 'left', 'Alpha', ttest_alpha)
            max_y = max([mean(pre) mean(post)]);
            scatter(mean([pre_i_val post_i_val]), max_y*1.06, '*', ...
                'MarkerEdgeColor', col, 'MarkerFaceColor', col,  'HandleVisibility', 'off')
            plot([pre_i_val post_i_val], 1.04*repelem(max_y,1,2), 'k-', 'HandleVisibility', 'off')
        end
    end
    
    % format subplot
    ylabel('Accuracy')
    ylim([0.5 0.85])
    ylabel('SVM Accuracy')
    xlim([0.5 length(sessions_to_use)+0.5])
    xticks(1:length(sessions_to_use))
    xticklabels(get_xtick_labs_colored(xticklabs, short_names, xtickcolors))
    xtickangle(45)
%     title(sprintf('Monkey %s', SVM(m).Code), 'Interpreter', 'none')
%     make_custom_patch_legend(colors(2:4,:), {'Posterior', 'Middle', 'Anterior'})
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',28, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
end
set(gcf,'color','w');

%% Load and plot: supp fig 3, intervals (barplot) (load)
% load data
clearvars
close all
fname = 'MaxMarta_fix_cat_xma2_SVM4_balancedInts.mat'; % 75-175, 175-275, 275-375.
path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
Data = load(fullfile(path1, fname));
[status, SVM] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Load and plot: supp fig 3, intervals (barplot) (plot)


% params
intervals_to_test = {[-150 -50], [75 175], [175 275], [275 375]}; % time window from cue on
marta_te_data = {'SC_filtered'};
marta_te_name = {'All TE'};
max_te_data = {'No_TEO_SC_filtered'};
max_te_name = {'All TE'};
neural_popn_subsets_monk = {marta_te_data, max_te_data};
neuron_subset_names_monk = {marta_te_name, max_te_name};

figure('Position', [300 300 630 500])
hold on

for m = 1:length(SVM)
    neural_popn_subsets = neural_popn_subsets_monk{m};
    neuron_subset_name = neuron_subset_names_monk{m};
    subplot(2,1,m)
    hold on
    
    if strcmp(SVM(m).Name, 'Marta_fix_cat_xma2')
        sessions_to_use = [7 9];
    elseif strcmp(SVM(m).Name, 'Max_fix_cat_xma2')
        sessions_to_use = [6 7];
    end
    
    for sub = 1:length(neural_popn_subsets)
        neuron_subset = neural_popn_subsets{sub};
        neuron_subset_name = neuron_subset_name{sub};
        
        bar_mat_means = zeros(length(intervals_to_test), 2);
        bar_mat_stds = zeros(length(intervals_to_test), 2);
        all_data = zeros(50*length(intervals_to_test), 2);
        xlabs = cell(1, length(intervals_to_test));
        ii = 1;
        for p = 1:length(intervals_to_test)
            interval_to_test = intervals_to_test{p};
            
            % get data
            kfl_id = get_good_interval_name2(interval_to_test, neuron_subset, 'KFL_abcat');
            acc_mat = 1 - [SVM(m).Sessions(sessions_to_use).(kfl_id)];
            
            % store
            bar_mat_means(p,:) = mean(acc_mat);
            bar_mat_stds(p,:) = std(acc_mat);
            all_data(ii:(ii+49),:) = acc_mat;
            xlabs{p} = sprintf('%d to %d', interval_to_test(1), interval_to_test(2));
            ii = ii + 50;
        end
        
        b = bar(bar_mat_means, 'grouped');
        pause(0.5);
        xoff = vertcat(b.XOffset);
        xvals = repmat(1:length(intervals_to_test),2,1) + xoff;
        errorbar(xvals', bar_mat_means, bar_mat_stds, 'ko');
        xticks(1:length(intervals_to_test))
        xticklabels(xlabs)
        ylabel('SVM accuracy')
        ylim([0.48 0.8])
        pause(0.5)
        legend({'Pre', 'Post'})
        
    end
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',26, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
end

%% ANOVA2 for Supp Fig 3
% params
intervals_to_test = {[-150 -50], [75 175], [175 275], [275 375]}; % time window from cue on
marta_te_data = {'SC_filtered'};
marta_te_name = {'All TE'};
max_te_data = {'No_TEO_SC_filtered'};
max_te_name = {'All TE'};
neural_popn_subsets_monk = {marta_te_data, max_te_data};
neuron_subset_names_monk = {marta_te_name, max_te_name};

for m = 1:length(SVM)
    neural_popn_subsets = neural_popn_subsets_monk{m};
    neuron_subset_name = neuron_subset_names_monk{m};
    if strcmp(SVM(m).Name, 'Marta_fix_cat_xma2')
        sessions_to_use = [7 9];
    elseif strcmp(SVM(m).Name, 'Max_fix_cat_xma2')
        sessions_to_use = [6 7];
    end
    
    for sub = 1:length(neural_popn_subsets)
        neuron_subset = neural_popn_subsets{sub};
        neuron_subset_name = neuron_subset_name{sub};
        
        bar_mat_means = zeros(length(intervals_to_test), 2);
        bar_mat_stds = zeros(length(intervals_to_test), 2);
        all_data = zeros(50*length(intervals_to_test), 2);
        xlabs = cell(1, length(intervals_to_test));
        ii = 1;
        for p = 1:length(intervals_to_test)
            interval_to_test = intervals_to_test{p};
            
            % get data
            kfl_id = get_good_interval_name2(interval_to_test, neuron_subset, 'KFL_abcat');
            acc_mat = 1 - [SVM(m).Sessions(sessions_to_use).(kfl_id)];
            
            % store
            bar_mat_means(p,:) = mean(acc_mat);
            bar_mat_stds(p,:) = std(acc_mat);
            all_data(ii:(ii+49),:) = acc_mat;
            xlabs{p} = sprintf('%d to %d', interval_to_test(1), interval_to_test(2));
            ii = ii + 50;
        end 
    end
    [pval, tbl, stats] = anova2(all_data, 50);
    multcompare(stats, 'estimate', 'row')
end

%% Load data (results of various analyses are available as .MAT files)
clearvars
close all
% fname = 'MaxMarta_fix_cat_xma2_SVM4_full.mat'; % 75-175 and 175-350.
% fname = 'MaxMarta_fix_cat_xma2_SVM4_full_TrialSubset.mat'; % 75-175 and 175-350, arrays, fixed trial nums across sessions
% fname = 'MaxMarta_fix_cat_xma2_SVM4_balancedInts.mat'; % 75-175, 175-275, 275-375.
% fname = 'MaxMarta_fix_cat_xma2_SVM4_only_GLM_Alpha05.mat';
% fname = 'MaxMarta_fix_cat_xma2_SVM4_BestGLMSubsets_TE.mat';
% fname = 'MaxMarta_fix_cat_xma2_SVM4_VisExc_TE.mat'; % using only units that were visually excited in non-parametric fine grained analysis
path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
Data = load(fullfile(path1, fname));
[status, SVM] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Set sessions to use and xtick names and colors
xtickcolors = cbrewer('qual', 'Dark2', 3);
xtickcolors = xtickcolors([3 2 1], :);
for m = 1:length(SVM)
    if strcmp(SVM(m).Name, 'Marta_fix_cat_xma2')
        SVM(m).Sessions_to_use = [1 2 3 5 6 7 9];
        SVM(m).XTickLabs= {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
        SVM(m).Code = 'R';
    elseif strcmp(SVM(m).Name, 'Max_fix_cat_xma2')
%         SVM(m).Sessions_to_use = [1 2 3 4 5 6 7];
%         SVM(m).XTickLabs = {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
        SVM(m).Sessions_to_use = [1 2 6 7];
        SVM(m).XTickLabs = {'Base 1', 'Base 2', 'Pre', 'Post'};
        SVM(m).Code = 'X';
    end
end

%% Define data params
all_greens = cbrewer('seq', 'Greens', 10);
greens_to_use = all_greens([4 7 10], :);


% MARTA
% arrays
marta_te_data = {'SC_filtered'};
marta_te_name = {'All TE'};
marta_subsets1 = {'SC_filtered', 'posterior_SC_filtered', 'middle_SC_filtered', 'anterior_SC_filtered'};
marta_subsets1_names = {'All TE', 'Posterior', 'Middle', 'Anterior'};
marta_subsets1_colors = [1 1 1; greens_to_use];
% GLM subsets
marta_subsets2 = {'SC_filtered', 'GLM_signf_1en3_all', 'GLM_NONsignf_1en3_all'};
marta_subsets2_names = {'All TE', 'TE, units siginifcant in GLM', 'TE, units NOT siginifcant in GLM'};
% image subsets
marta_subsets3 = { '240cd_all_NumImgMatched', '20cd_all_NumImgMatched', '240cd_all_NONNumImgMatched', '20cd_all_NONNumImgMatched'};
marta_subsets4 = {'VisExc_2only_all', 'VisExc_234_all'};
marta_subsets5 = {'2500TrialsSubset_SC_filtered', '2500TrialsSubset_posterior_SC_filtered',...
    '2500TrialsSubset_middle_SC_filtered', '2500TrialsSubset_anterior_SC_filtered'};

best_glm_subs = [1 2 4 8 16 32];
bestN_subsets = arrayfun(@(l) sprintf('BestGLMSubset_%d', l), best_glm_subs, 'UniformOutput', false);
marta_bestN_subsets = [bestN_subsets 'SC_filtered'];
marta_bestN_subsets_names = {'1', '2', '4', '8', '16', '32', 'TE'};

% MAX
% arrays
max_te_data = {'No_TEO_SC_filtered'};
max_te_name = {'All TE'};
max_subsets1 = {'No_TEO_SC_filtered', 'posterior_SC_filtered', 'middle_SC_filtered', 'anterior_SC_filtered'}; % arrays
max_subsets1_names = {'All TE', 'Posterior', 'Middle', 'Anterior'};
max_subsets1_colors = [1 1 1; greens_to_use];
max_subsets1a = {'No_TEO_SC_filtered', 'teo_SC_filtered'};
max_subsets1a_names = {'All TE', 'TEO'};
% GLM subsets
max_subsets2 = {'No_TEO_SC_filtered', 'GLM_signf_1en3_No_TEO', 'GLM_NONsignf_1en3_No_TEO'};
max_subsets2_names = {'All TE', 'TE, units siginifcant in GLM', 'TE, units NOT siginifcant in GLM'};
max_subsets2a = {'GLM_signf_1en3_teo', 'GLM_signf_1en3_No_TEO', 'GLM_NONsignf_1en3_teo', 'GLM_NONsignf_1en3_No_TEO'};
% image set subsets. These are spike-count filtered tho they don't say so.
max_subsets3 = { '240cd_teo_NumImgMatched', '240cd_No_TEO_NumImgMatched', '20cd_teo_NumImgMatched', '20cd_No_TEO_NumImgMatched'};
max_subsets4 = {'VisExc_2only_No_TEO', 'VisExc_234_No_TEO'};
max_subsets5 = {'2500TrialsSubset_No_TEO_SC_filtered', '2500TrialsSubset_posterior_SC_filtered',...
    '2500TrialsSubset_middle_SC_filtered', '2500TrialsSubset_anterior_SC_filtered'};
max_bestN_subsets = [bestN_subsets 'No_TEO_SC_filtered'];
max_bestN_subsets_names = {'1', '2', '4', '8', '16', '32', 'TE'};

%% (Supplemental fig) Plot SVM data across days, all intervals (diff subsets on diff subplots)

% Plotting Params
neural_popn_subsets_monk = {marta_subsets4, max_subsets4};
neural_popn_subset_names_monk = {marta_subsets4, max_subsets4};
% intervals_to_test = {[-175 0], [75 175], [175 350]}; % time window from cue on
intervals_to_test = {[75 175], [175 350]}; % time window from cue on

for m = 1:length(SVM)
    neural_popn_subsets = neural_popn_subsets_monk{m};
    neural_popn_subset_names = neural_popn_subset_names_monk{m};
    
    % select sessions to plot
    if strcmp(SVM(m).Name, 'Marta_fix_cat_xma2')
        sessions_to_use = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
    elseif strcmp(SVM(m).Name, 'Max_fix_cat_xma2')
        sessions_to_use = [1 2 3 5 6 7 8 9]; % skip base04 and sub03 (low trial count)
    end
    
    % set up figure and legend labels
%     figure('Position', [500 500 600 800])
    figure
    hold on
    labs = {};
    
    for sub = 1:length(neural_popn_subsets)
        neuron_subset = neural_popn_subsets{sub};
        neuron_subset_name = neural_popn_subset_names{sub};
        
        % configure subplots
        subplot(length(neural_popn_subsets), 1, sub)
        hold on
        %         if sub == length(neural_popn_subsets)
        %             subplot(length(neural_popn_subsets)+1, 1, [sub sub+1])
        %             hold on
        %         else
        %             subplot(length(neural_popn_subsets)+1, 1, sub)
        %             hold on
        %         end
        
        for p = 1:length(intervals_to_test)
            interval_to_test = intervals_to_test{p};
            
            % skip precue interval for GLM signf subsets
            if all(interval_to_test == [-175 0]) && (strcmp(neuron_subset, 'GLM_signf_1en3_teo') || strcmp(neuron_subset, 'GLM_signf_1en3_No_TEO'))
                plot(1:length(SVM(m).Sessions), repelem(0.5, length(SVM(m).Sessions)), '--')
                labs = [labs {sprintf('%d to %d ms',interval_to_test(1), interval_to_test(2))}];
                continue
            end
            
            % get data
            kfl_id = get_good_interval_name2(interval_to_test, neuron_subset, 'KFL_abcat');
            acc_mat = 1 - [SVM(m).Sessions(sessions_to_use).(kfl_id)]; % each kfl_vec is a column vector, can horiz concat across days. plot accuracy instead of loss
            
            % plot
            errorbar(mean(acc_mat), std(acc_mat), 'LineWidth', 2); % takes mean across iterations
            
            % format plot
            labs = [labs {sprintf('%d to %d ms',interval_to_test(1), interval_to_test(2))}];
            title(sprintf('%s',neuron_subset_name), 'Interpreter', 'none')
            ylim([0.45 0.80])
            xticklabels({})
            set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], ...
                'XMinorTick', 'off', 'YMinorTick', 'on', 'XGrid', 'on', 'YGrid', 'on',...
                'LineWidth', 1,  'fontsize',15, ...
                'fontname', 'Helvetica','FontWeight','Bold', 'LineWidth', 2,...
                'XColor', 'black', 'YColor', 'black')
            if sub == length(neural_popn_subsets)
                xlabel(sprintf('Session, %s', SVM(m).Name), 'Interpreter', 'none')
                ylabel('Accuracy')
                xticks(1:size(acc_mat,2))
                xticklabels({SVM(m).Sessions(sessions_to_use).ShortName})
                xtickangle(45)
                legend(labs, 'Location', 'eastoutside')
            end
        end
%         sgtitle(sprintf('%s KFL over sessions', SVM(m).Name), 'Interpreter', 'none')
        set(gcf,'color','w');
    end
end

%% (Supplemental fig) Plot SVM 175-350, across days/arrays
colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3;
    0          0        0    ];

% for best N 
% colors = [jet(6); [0 0 0]];
% colors(5,:) = colors(5,:)-[0.2 0.2 0.2]; % darken the yellow

% params
interval_to_test = [175 350];
include_shuffle = false;
overlay_te = true;
baseline_ttest_alpha = 0.05;
num_comps = 18;
ttest_alpha = baseline_ttest_alpha / num_comps;


% neural_popn_subsets_monk = {marta_bestN_subsets, max_bestN_subsets};
% neural_popn_subset_names_monk = {marta_bestN_subsets_names, max_bestN_subsets_names};

% neural_popn_subsets_monk = {marta_subsets4, max_subsets4};
% neural_popn_subset_names_monk = {{'Visually excited only', 'Visually excited any time'}, {'Visually excited only', 'Visually excited any time'}};

neural_popn_subsets_monk = {marta_subsets1, max_subsets1};
neural_popn_subset_names_monk = {marta_subsets1_names, max_subsets1_names};
% neural_popn_subset_colors_monk = {marta_subsets1_colors, max_subsets1_colors};

% neural_popn_subsets_monk = {marta_subsets5, max_subsets5};
% neural_popn_subset_names_monk = {marta_subsets5, max_subsets5};


% neural_popn_subsets_monk = {marta_subsets2, max_subsets2};
% neural_popn_subset_names_monk = {marta_subsets2_names, max_subsets2_names};

% prep fig
figure('Position', [200 200 500 600])
hold on

for m = 1:length(SVM)
    subplot(2,1,m)
    hold on
    
    
    % select sessions to plot
%     if strcmp(SVM(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_use = [1 2 7 9]; % skip base04 (outlier) and post01 (low trial count)
%     elseif strcmp(SVM(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_use = [1 2 6 7]; % This is base 01, base02, pre and post. skip base03,04,05 and sub03 (low trial count)
%     end
    sessions_to_use = SVM(m).Sessions_to_use;
    xticklabs = SVM(m).XTickLabs;
    short_names = {SVM(m).Sessions(sessions_to_use).ShortName};
    
    % select neural subsets to plot
    neural_popn_subsets = neural_popn_subsets_monk{m};
    neural_popn_subset_names = neural_popn_subset_names_monk{m};
    %     neural_popn_subset_colors = neural_popn_subset_colors_monk{m};
    
     
    % get baseline/pre/post session inds
    baseline_i_vals = [];
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        if ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Post', 'once'))
            post_i_val = i;
        elseif ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Pre', 'once'))
            pre_i_val = i;
        elseif ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Base', 'once'))
            baseline_i_vals = [baseline_i_vals i];
        end
    end
        
    % plot each subset, and run stats testing
    for sub = 1:length(neural_popn_subsets)
        neuron_subset = neural_popn_subsets{sub};
        neuron_subset_name = neural_popn_subset_names{sub};
        
        if ~overlay_te && strcmp(neuron_subset_name, 'All TE')
            continue
        end
        
        % get data and plot, all sessions to plot at once
        kfl_id = get_good_interval_name2(interval_to_test, neuron_subset, 'KFL_abcat');
        kfl_shuffled_id = get_good_interval_name2(interval_to_test, neuron_subset, 'KFL_abcat_SHUFFLED');
        acc_mat = 1 - [SVM(m).Sessions(sessions_to_use).(kfl_id)]; % each kfl_vec is a column vector, can horiz concat across days. plot accuracy instead of loss
        if include_shuffle
            shuffled_acc_mat = 1 - [SVM(m).Sessions(sessions_to_use).(kfl_shuffled_id)];
        end

        % plot the data
        wd = 1;
        col = colors(sub, :);
%         col = neural_popn_subset_colors(sub,:);
%         eb = errorbar(mean(acc_mat), std(acc_mat), 'LineWidth', wd,...
%             'DisplayName', neural_popn_subset_names{sub}, 'Color', col); % takes mean across iterations
        eb1 = errorbar(baseline_i_vals, mean(acc_mat(:,baseline_i_vals)),...
            std(acc_mat(:,baseline_i_vals)), ...
            'Color', col, 'LineWidth', 2,...
            'DisplayName', neural_popn_subset_names{sub}); % takes mean across iterations
        eb2 = errorbar([pre_i_val post_i_val], mean(acc_mat(:,[pre_i_val post_i_val])),...
            std(acc_mat(:,[pre_i_val post_i_val])),...
            'Color', col, 'LineWidth', 2,...
            'DisplayName', neural_popn_subset_names{sub}); 
        
        % add shuffled data if desired
        if include_shuffle
            ebs = errorbar(mean(shuffled_acc_mat), std(shuffled_acc_mat), '--', 'LineWidth', 0.5,...
                'DisplayName', neural_popn_subset_names{sub}, 'Color', col); % takes mean across iterations
            ebs.Bar.LineStyle = 'dotted';
            uistack(ebs, 'bottom')
        end
        
        % bring the TE line to the top, so it can be seen clearly
        if overlay_te && strcmp(neuron_subset_name, 'All TE')
            uistack(eb1, 'top')
            uistack(eb2, 'top')
        elseif overlay_te
            uistack(eb1, 'bottom')
            uistack(eb2, 'bottom')
        end
        
        % add pre/post stats testing
        pre = acc_mat(:, pre_i_val);
        post = acc_mat(:, post_i_val);
        if ttest2(pre, post, 'tail', 'left', 'Alpha', ttest_alpha)
            max_y = max([mean(pre) mean(post)]);
            scatter(mean([pre_i_val post_i_val]), max_y*1.06, '*', ...
                'MarkerEdgeColor', col, 'MarkerFaceColor', col,  'HandleVisibility', 'off')
            plot([pre_i_val post_i_val], 1.04*repelem(max_y,1,2), 'k-', 'HandleVisibility', 'off')
        end
        
        % now go session by session and do stats
        % first get post data
%         for i = 1:length(sessions_to_use)
%             sessn = sessions_to_use(i);
%             if ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Post', 'once'))
%                 post_data = 1 - SVM(m).Sessions(sessn).(kfl_id);
%             end
%         end
%         
%         % run stats testing
%         for i = 1:length(sessions_to_use)
%             sessn = sessions_to_use(i);
%             switch regexp(SVM(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
%                 case 'Base'
%                     data_to_compare = 1 - SVM(m).Sessions(sessn).(kfl_id);
%                     if ttest2(data_to_compare, post_data, 'tail', 'left', 'Alpha', ttest_alpha)
%                         scatter(i*1.02, mean(data_to_compare)*1.02, 40, col, 'filled', 'v', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off')
%                     end
%                 case 'Pre'
%                     data_to_compare = 1 - SVM(m).Sessions(sessn).(kfl_id);
%                     if ttest2(data_to_compare, post_data, 'tail', 'left', 'Alpha', ttest_alpha)
%                         scatter(i*1.02, mean(data_to_compare)*1.02, 40, col, 'filled', 'v', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off')
%                     end
%                 case 'Post'
%                     % dont test against itself
%                     continue
%             end
%         end
%         
    end
    
    % format subplot
    ylabel('Accuracy')
%     legend('Location', 'northeastoutside', 'Interpreter', 'none')
    %     yl = get(gca, 'YLim');
    %     yticks(0:.1:round(yl(2),1))
    xlim([0.5 length(sessions_to_use)+0.5])
    xticks(1:length(sessions_to_use))
    xticklabels(get_xtick_labs_colored(xticklabs, short_names, xtickcolors))
    xtickangle(45)
    title(sprintf('Monkey %s', SVM(m).Code), 'Interpreter', 'none')
%     make_custom_patch_legend(colors(2:4,:), {'Posterior', 'Middle', 'Anterior'})
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',18, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
end

% format figure
% sgtitle(sprintf('Interval %d to %d', interval_to_test(1), interval_to_test(2)))
set(gcf,'color','w');

%% Plot SVM 175-350, just TE
colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3;
    0          0        0    ];

% just full data
neural_popn_subsets_monk = {marta_te_data, max_te_data};
neural_popn_subset_names_monk = {marta_te_name, max_te_name};

% neural_popn_subsets_monk = {marta_bestN_subsets, max_bestN_subsets};
% neural_popn_subset_names_monk = {marta_bestN_subsets_names, max_bestN_subsets_names};

% neural_popn_subsets_monk = {marta_subsets1, max_subsets1};
% neural_popn_subset_names_monk = {marta_subsets1_names, max_subsets1_names};

include_shuffle = false;
plot_trial_balanced = false;
num_trials_to_use = 2500;
interval_to_test = [175 350];
% interval_to_test = [75 175];
% intervals_to_test = {[-150 -50], [75 175], [175 275], [275 375]}; % time window from cue on
baseline_ttest_alpha = 0.05;

num_comps = 9;
ttest_alpha = baseline_ttest_alpha / num_comps;

figure('Position', [200 200 500 600])
hold on

for m = 1:length(SVM)
    sessions_to_use = SVM(m).Sessions_to_use;
    xticklabs = SVM(m).XTickLabs;
    short_names = {SVM(m).Sessions(sessions_to_use).ShortName};
    
    subplot(2,1,m)
    hold on
    
    % select neural subsets to plot
    neural_popn_subsets = neural_popn_subsets_monk{m};
    neural_popn_subset_names = neural_popn_subset_names_monk{m};
    baseline_i_vals = [];
    
    % get baseline/pre/post session inds
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        if ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Post', 'once'))
            post_i_val = i;
        elseif ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Pre', 'once'))
            pre_i_val = i;
        elseif ~isempty(regexp(SVM(m).Sessions(sessn).ShortName, 'Base', 'once'))
            baseline_i_vals = [baseline_i_vals i];
        end
    end
        
    % plot each subset, and run stats testing
    for sub = 1:length(neural_popn_subsets)
        neuron_subset = neural_popn_subsets{sub};
        neuron_subset_name = neural_popn_subset_names{sub};
        
        % get data and plot, all sessions to plot at once
        if plot_trial_balanced
            kfl_id = get_good_interval_name2(interval_to_test, neuron_subset, sprintf('KFL_abcat_%dTrialsSubset',num_trials_to_use));
        else
            kfl_id = get_good_interval_name2(interval_to_test, neuron_subset, 'KFL_abcat');
        end
        acc_mat = 1 - [SVM(m).Sessions(sessions_to_use).(kfl_id)]; % each kfl_vec is a column vector, can horiz concat across days. plot accuracy instead of loss
        errorbar(baseline_i_vals, mean(acc_mat(:,baseline_i_vals)), std(acc_mat(:,baseline_i_vals)), 'Color', colors(sub,:), 'LineWidth', 2, 'DisplayName', neural_popn_subset_names{sub}); % takes mean across iterations
        errorbar([pre_i_val post_i_val], mean(acc_mat(:,[pre_i_val post_i_val])), std(acc_mat(:,[pre_i_val post_i_val])), 'Color', colors(sub,:), 'LineWidth', 2, 'DisplayName', neural_popn_subset_names{sub}); % takes mean across iterations
        
        if include_shuffle
            kfl_shuffled_id = get_good_interval_name2(interval_to_test, neuron_subset, 'KFL_abcat_SHUFFLED');
            shuffled_acc_mat = 1 - [SVM(m).Sessions(sessions_to_use).(kfl_shuffled_id)];
            errorbar(mean(shuffled_acc_mat), std(shuffled_acc_mat), '--', 'LineWidth', 0.5,...
                'DisplayName', neural_popn_subset_names{sub}, 'Color', colors(sub,:)); % takes mean across iterations
        end
        
        % pre post stats testing
        pre = acc_mat(:, pre_i_val);
        post = acc_mat(:, post_i_val);
        if ttest2(pre, post, 'tail', 'left', 'Alpha', ttest_alpha)
            max_y = max([mean(pre) mean(post)]);
            scatter(mean([pre_i_val post_i_val]), max_y*1.06, '*k',  'HandleVisibility', 'off')
            plot([pre_i_val post_i_val], 1.04*repelem(max_y,1,2), 'k-', 'HandleVisibility', 'off')
        end
        
        % run overall stats testing
%         for i = 1:length(sessions_to_use)
%             sessn = sessions_to_use(i);
%             switch regexp(SVM(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
%                 case 'Base'
%                     data_to_compare = 1 - SVM(m).Sessions(sessn).(kfl_id);
%                     if ttest2(data_to_compare, post_data, 'tail', 'left', 'Alpha', ttest_alpha)
%                         scatter(i*1.03, mean(data_to_compare)*1.03, '*k',  'HandleVisibility', 'off')
%                     end
%                 case 'Pre'
%                     data_to_compare = 1 - SVM(m).Sessions(sessn).(kfl_id);
%                     if ttest2(data_to_compare, post_data, 'tail', 'left', 'Alpha', ttest_alpha)
%                         scatter(i*1.03, mean(data_to_compare)*1.03, '*k', 'HandleVisibility', 'off')
%                     end
%                 case 'Post'
%                     post_data = 1 - SVM(m).Sessions(sessn).(kfl_id);
%                     continue
%             end
%         end
        
    end
    
    % format subplot
%     title(sprintf('Monkey %s', SVM(m).Code), 'Interpreter', 'none')
    ylim([0.5 0.85])
    ylabel('SVM Accuracy')
    xlim([0.5 length(sessions_to_use)+0.5])
    xticks(1:length(sessions_to_use))
    xticklabels(get_xtick_labs_colored(xticklabs, short_names, xtickcolors))
    xtickangle(45)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',18, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
    %     legend
end

% format figure
set(gcf,'color','w');

%% ******* Scatter plots of VR *******

%% Load data
clearvars
close all
path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
% fname = 'MaxMarta_GLMs_VR_contrast.mat';
fname = 'MaxMarta_VR_img_TTest.mat';
Data = load(fullfile(path1, fname));
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Set sessions to use and xtick names and colors
xtickcolors = cbrewer('qual', 'Dark2', 3);
epoch_colors = xtickcolors([3 2 1], :);
matlab_colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3;
    0          0        0    ];
for m = 1:length(Monkeys)
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
%         Monkeys(m).Sessions_to_use = [1 2 3 5 6 7 9];
%         Monkeys(m).XTickLabs= {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
        Monkeys(m).Sessions_to_use = [1 7 9];
        Monkeys(m).XTickLabs = {'Base 1', 'Pre', 'Post'};
        Monkeys(m).Code = 'R';
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
%         Monkeys(m).Sessions_to_use = [1 2 3 4 5 6 7];
%         Monkeys(m).XTickLabs = {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
        Monkeys(m).Sessions_to_use = [1 6 7];
        Monkeys(m).XTickLabs = {'Base 1', 'Pre', 'Post'};
        Monkeys(m).Code = 'X';
    end
end

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

marta_xls = './RecordingMarta.xlsx';
max_xls = './RecordingMax.xlsx';
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

%% Load waveform info if desired

WF = load('XMA2/Monkey_structs/Marta_fix_cat_xma2_ni_wf_summaryOnly');
WF = WF.WF;
class_bounds = [0 0.3 0.4 1];
for m = 1
    for i = 1:9
        for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
            wid = WF(m).Sessions(i).UnitInfo(j).Waveform_width_ms;
            if wid < class_bounds(1)
                Monkeys(m).Sessions(i).UnitInfo(j).Waveform_wid_class = 1;
                
            elseif wid >= class_bounds(1) && wid < class_bounds(2)
                Monkeys(m).Sessions(i).UnitInfo(j).Waveform_wid_class  = 2;
                
            elseif wid >= class_bounds(2) && wid < class_bounds(3)
                Monkeys(m).Sessions(i).UnitInfo(j).Waveform_wid_class  = 3;
                
            elseif wid >= class_bounds(3) && wid < class_bounds(4)
                Monkeys(m).Sessions(i).UnitInfo(j).Waveform_wid_class  = 4;
                
            elseif wid >= class_bounds(4)
                Monkeys(m).Sessions(i).UnitInfo(j).Waveform_wid_class  = 5;
            end
        end
        
    end
end

%% (raw data) Scatter avg cat/dog spike rates
interval_to_plot = [175 350];
area_to_plot = 'te';
% array_cols =  [0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250;
%     0.4940    0.1840    0.5560];
colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3    ;
    0          0        0    ];
color_by_array = true;
figure('Position', [400 400 1300 1000])
hold on

for m = 1:length(Monkeys)
    
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
        %         sessions_to_use = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_use = [1 2 7 9];
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
        %         sessions_to_use = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
        sessions_to_use = [1 2 6 7];
    end
    
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        
        % only take te data
        arrays = {Monkeys(m).Sessions(sessn).UnitInfo.Location};
        areas = {Monkeys(m).Sessions(sessn).UnitInfo.Area};
        area_bool = strcmp(areas, area_to_plot);
        
        % get data
        xid = get_good_interval_name2(interval_to_plot, 'full', 'X');
        X = Monkeys(m).Sessions(sessn).(xid);
        Y = Monkeys(m).Sessions(sessn).Session_Y_catg;
        cat_means = mean(X(Y==1, area_bool));
        dog_means = mean(X(Y==2, area_bool));
        xvals = 0:10;
        
        % get array data
        arrays_subset = arrays(area_bool);
        uq = unique(arrays_subset);
        cols = zeros(length(arrays_subset),3);
        for j = 1:length(uq)
            inds = find(strcmp(arrays_subset, uq{j}));
            cols(inds, :) = repmat(colors(j,:), numel(inds), 1);
        end
        
        % plot data
        subplot(length(Monkeys),length(sessions_to_use),...
            length(sessions_to_use)*(m-1) + mod(i-1, length(sessions_to_use)) + 1)
        hold on
        
        if color_by_array
            scatter(dog_means, cat_means, 30, cols, 'o', 'filled')
            %             make_custom_patch_legend(colors(1:length(uq), :), uq)
        else
            scatter(dog_means, cat_means, 'bo')
        end
        plot(xvals, xvals, 'k--', 'HandleVisibility', 'off')
        %         plot(xvals, intercept + xvals*slope, 'r--')
        
        
        % format plots
        if i == 1 && m == 1
            xlabel('Avg dog spike count')
            ylabel('Avg cat spike count')
        end
        if m == 1
            xlim([0 6])
            ylim([0 6])
        elseif m == 2
            xlim([0 10])
            ylim([0 10])
        end
        xticks([0 5 10])
        yticks([0 5 10])
        title(sprintf('%s',Monkeys(m).Sessions(sessn).ShortName), 'FontWeight', 'normal')
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',18, ...
            'fontname', 'Helvetica', 'fontweight', 'normal', ...
            'XColor', 'black', 'YColor', 'black')
        make_custom_patch_legend(colors(1:length(uq), :), uq);
    end
    sgtitle(sprintf('Area %s, %d to %d, image VR propns,(alpha = %0.2f)',...
        area_to_plot, interval_to_plot(1), interval_to_plot(2), vr_alpha), 'Interpreter', 'none')
    set(gcf, 'Color', 'w')
end

%% (OLD) Run linear regression and plot the data
alpha = 0.05;
test_intervals = {[75 175]};
baseline_intervals = {[-150 -50]};
for m = 1:length(Monkeys)
    sessions_to_use = 1:length(Monkeys(m).Sessions);
    for p = 1:length(test_intervals)
        test_int = test_intervals{p};
        baseline_int = baseline_intervals{p};
        tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img');
        bid = get_good_interval_name2(baseline_int, '', '');
        vr_id = strcat(tid,bid);
        
        areas_to_plot = unique({Monkeys(m).Sessions(sessions_to_use(1)).UnitInfo.Area});
        for a = 1:length(areas_to_plot)
            
            f1 = figure('Position', [400 400 1300 1000]);
            hold on
            
            fqq = figure('Position', [400 400 1300 1000]);
            hold on
            
            for i = 1:length(sessions_to_use)
                sessn = sessions_to_use(i);
                
                % get data
                data_mat = Monkeys(m).Sessions(sessn).(vr_id); % stimuli x units
                area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, areas_to_plot{a});
                data_mat = data_mat(:, area_bool);
                cat_props = sum(data_mat(1:260,:) < alpha) / 260;
                dog_props = sum(data_mat(261:520,:) < alpha) / 260;
                
                % run linear model and plot results
                figure(f1)
                subplot(3,4,i)
                hold on
                scatter(dog_props, cat_props, 'bo')
                plot(0:0.01:1, -0:0.01:1, 'k--')
                xlabel('Proportion dogs VR')
                ylabel('Proportion cats VR')
                
                % run model
                mdl = fitlm(dog_props, cat_props);
                xvals = 0:0.01:1;
                plot(xvals, xvals*mdl.Coefficients{'x1', 'Estimate'} + mdl.Coefficients{'(Intercept)', 'Estimate'}, 'r--')
                title(sprintf('Session %s, R2 = %0.2f', Monkeys(m).Sessions(sessn).ShortName, mdl.Rsquared.Ordinary))
                xlim([0 1])
                ylim([0 1])
                
                % qq plot
                figure(fqq)
                subplot(3,4,i)
                hold on
                qqplot(mdl.Residuals.Raw)
                title(sprintf('Session %s, R2 = %0.2f', Monkeys(m).Sessions(sessn).ShortName, mdl.Rsquared.Ordinary))
            end
            figure(f1)
            sgtitle(sprintf('%s, area %s, %d to %d ms, image VR propns,(alpha = %0.2f)', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), alpha), 'Interpreter', 'none')
            figure(fqq)
            sgtitle(sprintf('%s, area %s, %d to %d ms, qq plot (alpha = %0.2f)', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), alpha), 'Interpreter', 'none')
        end
    end
end

%% (OLD) Calculate cat/dog selectivity indexes (cat prop - dog prop)

alpha = 0.05;
alpha_string_scale_factor = 100;
test_intervals = {[75 175], [175 350]};
baseline_intervals = {[-150 -50], [-175 0]};
for m = 1:length(Monkeys)
    sessions_to_use = 1:length(Monkeys(m).Sessions);
    for p = 1:length(test_intervals)
        test_int = test_intervals{p};
        baseline_int = baseline_intervals{p};
        tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img');
        bid = get_good_interval_name2(baseline_int, '', '');
        vr_id = strcat(tid,bid);
        
        areas_to_plot = unique({Monkeys(m).Sessions(sessions_to_use(1)).UnitInfo.Area});
        for a = 1:length(areas_to_plot)
            for i = 1:length(sessions_to_use)
                sessn = sessions_to_use(i);
                
                % get data
                data_mat = Monkeys(m).Sessions(sessn).(vr_id); % stimuli x units
                area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, areas_to_plot{a});
                data_mat = data_mat(:, area_bool);
                cat_props = sum(data_mat(1:260,:) < alpha) / 260;
                dog_props = sum(data_mat(261:520,:) < alpha) / 260;
                
                % store differences as measure of spread around y=x
                diff_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropDiffs_%s_Alpha%d', areas_to_plot{a}, alpha*alpha_string_scale_factor));
                Monkeys(m).Sessions(sessn).(diff_id) = cat_props - dog_props;
            end
        end
    end
end

%% (Prop VR) Use major-axis regression to assess slopes, store residuals
% MAR alpha is not dynamic for figures -- if you want to use a different
% MAR alpha (for slope CI's), come back here and re-run the analysis first.

vr_alpha = 0.05;
mar_alpha = 0.05;
alpha_string_scale_factor = 100;
test_intervals = {[75 175], [175 350]};
baseline_intervals = {[-150 -50], [-175 0]};
% test_intervals = {[175 350]};
% baseline_intervals = {[-175 0]};
make_plots = false;


for m = 1:length(Monkeys)
    % for m = 1
    sessions_to_use = 1:length(Monkeys(m).Sessions);
    %     if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
    % %         sessions_to_use = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
    %         sessions_to_use = [1 2 7 9];
    %     elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
    % %         sessions_to_use = [1 2 3 4 5 6 7 9]; % skip base04 and sub03 (low trial count)
    %         sessions_to_use = [1 2 6 7];
    %     end
    
    for p = 1:length(test_intervals)
        test_int = test_intervals{p};
        baseline_int = baseline_intervals{p};
        %         tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img');
        tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img_TTEST');
        bid = get_good_interval_name2(baseline_int, '', '');
        vr_id = strcat(tid,bid);
        
        areas_to_plot = unique({Monkeys(m).Sessions(sessions_to_use(1)).UnitInfo.Area});
        for a = 1:length(areas_to_plot)
            if make_plots
                figure('Position', [400 400 1300 1000])
                hold on
            end
            for i = 1:length(sessions_to_use)
                sessn = sessions_to_use(i);
                
                % get data
                data_mat = Monkeys(m).Sessions(sessn).(vr_id); % stimuli x units
                area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, areas_to_plot{a});
                data_mat = data_mat(:, area_bool);
                cat_props = sum(data_mat(1:260,:) < vr_alpha) / 260;
                dog_props = sum(data_mat(261:520,:) < vr_alpha) / 260;
                
                % calculate correlation
                x = corrcoef(cat_props, dog_props);
                fprintf('%s, %s, %s, interval %d to %d, correlation %0.2f \n', ...
                    Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName, areas_to_plot{a}, test_int(1), test_int(2), x(1,2))
                
                
                % get angle of best-fit line with major axis regression
                % method. I checked this function against a reference, it's
                % correct.
                [coeffs,coeff_ints,~,~] = maregress(dog_props, cat_props, mar_alpha);
                theta = rad2deg(atan(coeffs(2)));
                theta_bounds = rad2deg(atan(coeff_ints(2,:)));
                intercept = coeffs(1);
                intercept_bounds = coeff_ints(1,:);
                % calculate distance from points to regression line:
                % D = abs(ax0+by0+c) / sqrt(a^2+b^2), where point is
                % (x0,y0), and line is ax+by+c = 0.
                % We have y=mx+c form, so a=coeffs(2), b = -1, and c =
                % coeffs(1).
                perpendicular_residuals = abs(coeffs(2)*dog_props + -1*cat_props + coeffs(1)) / sqrt(coeffs(2)^2 + 1);
                
                % store data
                propn_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_CatDogPropns_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                slope_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Slope_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                intercept_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Intercept_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                resid_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Resids_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                
                Monkeys(m).Sessions(sessn).(propn_id) = [cat_props; dog_props]; % 2 x units
                Monkeys(m).Sessions(sessn).(resid_id) = perpendicular_residuals; % 1 x units
                Monkeys(m).Sessions(sessn).(slope_id) = [theta_bounds(1) theta theta_bounds(2)]; % 1 x 3
                Monkeys(m).Sessions(sessn).(intercept_id) = [intercept_bounds(1) intercept intercept_bounds(2)]; % 1 x 3
                
                
                % ditto but with 95 % outliers removed. Fit the line
                % without outliers, but still calculate residuals using all
                % points.
                diffs = cat_props - dog_props;
                outliers = abs(diffs) > prctile(abs(diffs), 95);
                [coeffs,coeff_ints,~,~] = maregress(dog_props(~outliers), cat_props(~outliers), mar_alpha);
                theta = rad2deg(atan(coeffs(2)));
                theta_bounds = rad2deg(atan(coeff_ints(2,:)));
                intercept = coeffs(1);
                intercept_bounds = coeff_ints(1,:);
                perpendicular_residuals = abs(coeffs(2)*dog_props + -1*cat_props + coeffs(1)) / sqrt(coeffs(2)^2 + 1);
                slope_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Slope_%s_Alpha%d_OutliersRm', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                intercept_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Intercept_%s_Alpha%d_OutliersRm', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                resid_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Resids_%s_Alpha%d_OutliersRm', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                Monkeys(m).Sessions(sessn).(resid_id) = perpendicular_residuals; % 1 x units
                Monkeys(m).Sessions(sessn).(slope_id) = [theta_bounds(1) theta theta_bounds(2)]; % 1 x 3
                Monkeys(m).Sessions(sessn).(intercept_id) = [intercept_bounds(1) intercept intercept_bounds(2)]; % 1 x 3
                
                
                % plot
                if make_plots
                    subplot(3,4,i)
                    hold on
                    scatter(dog_props, cat_props, 'bo')
                    xvals = 0:0.01:1;
                    plot(xvals, xvals, 'k--')
                    plot(xvals, coeffs(1) + xvals*coeffs(2), 'r--')
                    xlabel('Proportion dogs VR')
                    ylabel('Proportion cats VR')
                    title(sprintf('Session %s, theta = %0.1f (%0.1f - %0.1f)', Monkeys(m).Sessions(sessn).ShortName, theta, theta_bounds(1), theta_bounds(2)))
                end
                % with interaction method
                % make data table.
                %                 xvals = [dog_props cat_props];
                %                 group_vals = [repelem(1, numel(dog_props)) repelem(2, numel(dog_props))];
                %                 yvals = [cat_props cat_props];
                %                 tbl = table(xvals', group_vals', yvals', 'VariableNames', {'Dog_prop', 'Group', 'Cat_prop'});
                %                 % run linear model
                %                 mdl = fitlm(tbl, 'Cat_prop ~ Dog_prop * Group');
            end
            if make_plots
                sgtitle(sprintf('%s, area %s, %d to %d, image VR propns,(alpha = %0.2f)', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), vr_alpha), 'Interpreter', 'none')
            end
        end
    end
end

%% Notes: VR units to investigate img selectivity of
%{
-Max, post(7): index in TE: 145, 66, 129, 103, 122
-Marta, post (9): 20, 76, 18, 74 (all dogs), 43, 25, 81; 88, 91

%}

%% Transform max TE inds into unitinfo inds

inds_in_te = [145, 66, 129, 103, 122];
m = 2;
sessn = 7;
unit_locs = {Monkeys(m).Sessions(sessn).UnitInfo.Area};
te_in_full = find(strcmp(unit_locs, 'te'));
inds_in_full = te_in_full(inds_in_te);

%% Fig 3 / Supp Fig 4: plot selected MAR results with scatterplots

interval_to_plot = [175 350];
vr_alpha = 0.05;
alpha_string_scale_factor = 100;
area_to_plot = 'te';
color_by_waveform_class = false;
colors = cbrewer('qual', 'Set1', 5);
figure('Position', [400 400 800 500])
% figure('Position', [400 400 1600 500])
hold on

for m = 1:length(Monkeys)
% for m = 1
    
%     if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
%         %         sessions_to_use = [1 2 3 5 6 7 9]; % skip base04 (outlier) and post 01 (low trial count)
%         %         sessions_to_use = [1 2 7 9];
%         sessions_to_use = [1 7 9];
%     elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
%         %         sessions_to_use = [1 2 3 5 6 7]; % skip base04 due to all zeros on scatter plot (low trial count)
%         %         sessions_to_use = [1 2 6 7];
%         sessions_to_use = [1 6 7];
%     end
    sessions_to_use = Monkeys(m).Sessions_to_use;
    
%     sessions_to_use = Monkeys(m).Sessions_to_use;
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        
        % get data
        propn_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_CatDogPropns_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        slope_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_Slope_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        intercept_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_Intercept_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
%         slope_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Slope_%s_Alpha%d_OutliersRm', area_to_plot, vr_alpha*alpha_string_scale_factor));
%         intercept_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Intercept_%s_Alpha%d_OutliersRm', area_to_plot, vr_alpha*alpha_string_scale_factor));
%         propn_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_CatDogPropns_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        t = Monkeys(m).Sessions(sessn).(propn_id);
        cat_props = t(1,:);
        dog_props = t(2,:);
        slope = tan(deg2rad(Monkeys(m).Sessions(sessn).(slope_id)(2)));
        intercept = Monkeys(m).Sessions(sessn).(intercept_id)(2);
        xvals = 0:0.01:1;
        
        % get color / plotting params
        switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
            case 'Base'
                col = epoch_colors(1,:);
                thk = 1;
            case 'Pre'
                col = epoch_colors(2,:);
                thk = 2;
                pre_data = perpendicular_residuals;
            case 'Post'
                col = epoch_colors(3,:);
                thk = 2;
                post_data = perpendicular_residuals;
        end
                
        
        % plot data
        subplot(length(Monkeys),length(sessions_to_use),...
            length(sessions_to_use)*(m-1) + mod(i-1, length(sessions_to_use)) + 1)
        %         subplot(1, length(sessions_to_use), i)
        hold on
        
        % red, blue, green, purple, orange
        if color_by_waveform_class
            area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, area_to_plot);
            cols = [Monkeys(m).Sessions(sessn).UnitInfo(area_bool).Waveform_wid_class];
            scatter(dog_props, cat_props, 20, colors(cols-1,:), 'filled', 'o')
        else
            scatter(dog_props, cat_props, 'bo')
        end
        
        
        % color in example points
%         if i == 1 && m == 1
%             scatter(dog_props(14), cat_props(14), 'go', 'filled')
%             scatter(dog_props(2), cat_props(2), 'go', 'filled')
%         end


        plot(xvals, xvals, 'k--', 'LineWidth', 2)
        plot(xvals, intercept + xvals*slope, 'r--', 'LineWidth', 2)
        
        % format plots
        if i == 1 && m == 1
            %         if i == 1
            xlabel('Proportion dogs VR')
            ylabel('Proportion cats VR')
        end
        xlim([0 1])
        ylim([0 1])
        xticks([0 0.5 1])
        yticks([0 0.5 1])
        %         title(sprintf('%s', Monkeys(m).Sessions(sessn).ShortName), 'FontWeight', 'normal')
        title(Monkeys(m).XTickLabs{i}, 'Color', col, 'FontWeight', 'normal')
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',20, 'YGrid', 'on', 'XGrid', 'on', ...
            'XColor', 'black', 'YColor', 'black')
    end
    %     sgtitle(sprintf('Area %s, %d to %d, image VR propns,(alpha = %0.2f)',...
    %         area_to_plot, interval_to_plot(1), interval_to_plot(2), vr_alpha), 'Interpreter', 'none')
    set(gcf, 'Color', 'w')
end

%% Extra: Same as above VR plot, separate by array

interval_to_plot = [175 350];
tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img_TTEST');
bid = get_good_interval_name2(baseline_int, '', '');
vr_id = strcat(tid,bid);
vr_alpha = 0.05;
areas_to_plot = {'anterior', 'middle', 'posterior'};
% array_cols =  [0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250;
%     0.4940    0.1840    0.5560];
colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3    ;
    0          0        0    ];


for m = 1:length(Monkeys)
    
    % new fig for each monk
    figure('Position', [400 400 1300 1000])
    hold on
    
    
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
        %         sessions_to_use = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_use = [1 2 7 9];
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
        %         sessions_to_use = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
        sessions_to_use = [1 2 6 7];
    end
    
    for a = 1:length(areas_to_plot)
        area = areas_to_plot{a};
        
        for i = 1:length(sessions_to_use)
            sessn = sessions_to_use(i);

            % get data
            data_mat = Monkeys(m).Sessions(sessn).(vr_id); % stimuli x units
            area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Location}, area);
            data_mat = data_mat(:, area_bool);
            cat_props = sum(data_mat(1:260,:) < vr_alpha) / 260;
            dog_props = sum(data_mat(261:520,:) < vr_alpha) / 260;
            xvals = 0:0.1:1;
            
            % plot data
            subplot(length(areas_to_plot),length(sessions_to_use),...
                length(sessions_to_use)*(a-1) + mod(i-1, length(sessions_to_use)) + 1)
            hold on
            scatter(dog_props, cat_props, 30, colors(a,:), 'o')
            plot(xvals, xvals, 'k--')
%             make_custom_patch_legend(array_cols, uq)

            % format plots
            if i == 1 && m == 1
                xlabel('Proportion dogs VR')
                ylabel('Proportion cats VR')
            end
%             if a == 1
%                 xlim([0 0.2])
%                 ylim([0 0.2])
%                 xticks([0 0.1 0.2])
%                 yticks([0 0.1 0.2])
%             else
                xlim([0 1])
                ylim([0 1])
                xticks([0 0.5 1])
                yticks([0 0.5 1])
%             end
            
            %         set(gca, 'XScale', 'log', 'YScale', 'log')
            title(sprintf('%s, %s', Monkeys(m).Sessions(sessn).ShortName, area), 'FontWeight', 'normal')
            set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
                'XMinorTick', 'off', 'YMinorTick', 'off',...
                'fontsize',18, ...
                'fontname', 'Helvetica', 'fontweight', 'normal', ...
                'XColor', 'black', 'YColor', 'black')
            %         make_custom_patch_legend(colors(1:length(uq), :), uq);
        end
    end
    sgtitle(sprintf('Monkey %s, %d to %d, image VR propns,(alpha = %0.2f)',...
        Monkeys(m).Name, interval_to_plot(1), interval_to_plot(2), vr_alpha), 'Interpreter', 'none')
    set(gcf, 'Color', 'w')
end

%% Extra: Same as above VR plot, color by array
interval_to_plot = [175 350];
vr_alpha = 0.05;
alpha_string_scale_factor = 100;
area_to_plot = 'te';
% array_cols =  [0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250;
%     0.4940    0.1840    0.5560];
colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3    ;
    0          0        0    ];
figure('Position', [400 400 1300 1000])
hold on

for m = 1:length(Monkeys)
    
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
        %         sessions_to_use = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_use = [1 2 7 9];
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
        %         sessions_to_use = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
        sessions_to_use = [1 2 6 7];
    end
    
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        
        % get data
        propn_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_CatDogPropns_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        slope_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_Slope_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        intercept_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_Intercept_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        t = Monkeys(m).Sessions(sessn).(propn_id);
        cat_props = t(1,:);
        dog_props = t(2,:);
        slope = tan(deg2rad(Monkeys(m).Sessions(sessn).(slope_id)(2)));
        intercept = Monkeys(m).Sessions(sessn).(intercept_id)(2);
        xvals = 0:0.01:1;
        
        % get array data
        arrays = {Monkeys(m).Sessions(sessn).UnitInfo.Location};
        areas = {Monkeys(m).Sessions(sessn).UnitInfo.Area};
        area_bool = strcmp(areas, area_to_plot);
        arrays_subset = arrays(area_bool);
        uq = unique(arrays_subset);
        cols = zeros(length(arrays_subset),3);
        for j = 1:length(uq)
            inds = find(strcmp(arrays_subset, uq{j}));
            cols(inds, :) = repmat(colors(j,:), numel(inds), 1);
        end
        
        % plot data
        subplot(length(Monkeys),length(sessions_to_use),...
            length(sessions_to_use)*(m-1) + mod(i-1, length(sessions_to_use)) + 1)
        hold on
        scatter(dog_props, cat_props, 30, cols, 'o', 'filled')
        plot(xvals, xvals, 'k--')
        %         plot(xvals, intercept + xvals*slope, 'r--')
        plot(xvals, intercept + xvals*slope, '--', 'Color', [0.4 0.4 0.4])
        %         make_custom_patch_legend(array_cols, uq)
        
        % format plots
        if i == 1 && m == 1
            xlabel('Proportion dogs VR')
            ylabel('Proportion cats VR')
        end
        xlim([0 1])
        ylim([0 1])
        xticks([0 0.5 1])
        yticks([0 0.5 1])
        %         set(gca, 'XScale', 'log', 'YScale', 'log')
        title(sprintf('%s',Monkeys(m).Sessions(sessn).ShortName), 'FontWeight', 'normal')
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',18, ...
            'fontname', 'Helvetica', 'fontweight', 'normal', ...
            'XColor', 'black', 'YColor', 'black')
        %         make_custom_patch_legend(colors(1:length(uq), :), uq);
    end
    sgtitle(sprintf('Area %s, %d to %d, image VR propns,(alpha = %0.2f)',...
        area_to_plot, interval_to_plot(1), interval_to_plot(2), vr_alpha), 'Interpreter', 'none')
    set(gcf, 'Color', 'w')
end

%% Fig 3: Plot CDFs of MAR residuals, and do stats

vr_alpha = 0.05;
alpha_string_scale_factor = 100;
cdf_alpha = 0.05;
ranksum_alpha = 0.025;
test_intervals = {[175 350]};

figure('Position', [400 400 260 430])
hold on

for m = 1:length(Monkeys)
    subplot(2,1,m)
    hold on
    
%     if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_use = [1 2 3 5 6 7 9]; % skip base04 (outlier) and post 01 (low trial count)
%         %         sessions_to_use = [1 2 7 9];
%     elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_use = [1 2 3 5 6 7]; % skip base04 due to all zeros on scatter plot (low trial count)
%         %         sessions_to_use = [1 2 6 7];
%     end
    
    sessions_to_use = Monkeys(m).Sessions_to_use;
    
    for p = 1:length(test_intervals)
        test_int = test_intervals{p};
        %         areas_to_plot = unique({Monkeys(m).Sessions(sessions_to_use(1)).UnitInfo.Area});
        areas_to_plot = {'te'};
        
        for a = 1:length(areas_to_plot)
            for i = 1:length(sessions_to_use)
                sessn = sessions_to_use(i);
                
                % use selectivity index (cats - dogs)
                %                 diff_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropDiffs_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                %                 diffs = abs(Monkeys(m).Sessions(sessn).(diff_id));
                %                 [f, x] = ecdf(diffs);
                %                 xlabel_str = 'Selectivity (|prop_cat - prop_dog|)';
                
                % use residuals from major axis regression
                propn_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_CatDogPropns_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                resid_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Resids_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                t = Monkeys(m).Sessions(sessn).(propn_id);
                cat_props = t(1,:);
                dog_props = t(2,:);
                perpendicular_residuals = Monkeys(m).Sessions(sessn).(resid_id);
                
                % control: restrict to only neurons with low num VR imgs
%                 units_to_use = (cat_props/2 + dog_props/2) < 0.25;
%                 perpendicular_residuals = perpendicular_residuals(units_to_use);
                
                [f, x] = ecdf(perpendicular_residuals);
                xlabel_str = 'Residuals from major axis regression';
                
                switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
                    case 'Base'
%                         col = [0.4 0.4 0.4];
                        col = epoch_colors(1,:);
                        thk = 1;
                    case 'Pre'
                        %                         col = 'r';
                        col = [0.4 0.4 0.4];
                        col = epoch_colors(2,:);
                        thk = 2;
                        pre_data = perpendicular_residuals;
                    case 'Post'
%                         col = 'g';
                        col = epoch_colors(3,:);
                        thk = 2;
                        post_data = perpendicular_residuals;
                    case 'Sub'
                        col = [0.8 0.8 0.8];
                        thk = 1;
                end
                plot(x, f, 'Color', col, 'LineWidth', thk, 'DisplayName', Monkeys(m).Sessions(sessn).ShortName)
            end
            
            [pval,h] = ranksum(pre_data, post_data, 'Alpha', ranksum_alpha, 'Tail', 'left');
            fprintf('%s, %s, %d to %d, pre vs. post ranksum (one-sided non-parametric), hyp %d, pval of %d \n', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), h, pval)
            
            [h, pval] = ttest2(pre_data, post_data, 'tail', 'left');
            fprintf('%s, %s, %d to %d, pre vs. post ttest2 (one-sided), hyp %d, pval of %d \n', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), h, pval)
            
            %             title(sprintf('%s, %s, %d to %d, VR alpha %0.2f', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), vr_alpha), 'Interpreter', 'none')
            %             legend
            xlabel('Residuals')
            xlim([0 0.2])
            ylabel('Empirical CDF')
            set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
                'XMinorTick', 'off', 'YMinorTick', 'off',...
                'fontsize',18, ...
                'fontname', 'Helvetica', 'fontweight', 'normal', ...
                'XColor', 'black', 'YColor', 'black')
        end
    end
end
set(gcf, 'Color', 'w')

%% Fig 3: Compare slopes of MAR regressions across days
vr_alpha = 0.05;
alpha_string_scale_factor = 100;
test_intervals = {[175 350]};
area_to_plot = 'te';

% make fig
figure('Position', [400 400 600 800]);
hold on
for m = 1:length(Monkeys)
    subplot(2,1,m)
    hold on
    
    % select sessions
%     if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_use = [1 2 3 5 6 7 9]; % skip base04 (outlier) and post 01 (low trial count)
%         %         sessions_to_use = [1 2 7 9];
%     elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_use = [1 2 3 4 5 6 7]; % skip base04 due to all zeros on scatter plot (low trial count)
%         %         sessions_to_use = [1 2 6 7];
%     end
    sessions_to_use = Monkeys(m).Sessions_to_use;
    baseline_i_vals = [];
    
    % get baseline/pre/post session inds
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        if ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post', 'once'))
            post_i_val = i;
        elseif ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre', 'once'))
            pre_i_val = i;
        elseif ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Base', 'once'))
            baseline_i_vals = [baseline_i_vals i];
        end
    end
    
    for p = 1:length(test_intervals)
        test_int = test_intervals{p};
        slope_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Slope_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
%         slope_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Slope_%s_Alpha%d_OutliersRm', area_to_plot, vr_alpha*alpha_string_scale_factor));
        
        % pre allocate
        theta_means = zeros(1,length(sessions_to_use));
        theta_bounds = zeros(2,length(sessions_to_use));
        
        % collect data
        for i = 1:length(sessions_to_use)
            sessn = sessions_to_use(i);
            t = Monkeys(m).Sessions(sessn).(slope_id);
            theta_means(i) = t(2);
            theta_bounds(:,i) = [t(1); t(3)];
            if ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post', 'once'))
                post_data = t;
            end
        end
        
        % plot data
        errorbar(baseline_i_vals, theta_means(baseline_i_vals), diff(theta_bounds(:,baseline_i_vals))/2, 'LineWidth', 2, ...
            'DisplayName', sprintf('%s, %d to %d', Monkeys(m).Name, test_int(1), test_int(2)),...
            'Color', matlab_colors(1,:))
        errorbar([pre_i_val post_i_val], theta_means([pre_i_val post_i_val]), diff(theta_bounds(:, [pre_i_val post_i_val]))/2, 'LineWidth', 2, ...
            'DisplayName', sprintf('%s, %d to %d', Monkeys(m).Name, test_int(1), test_int(2)),...
            'Color', matlab_colors(1,:))
        
        % stats testing
        for i = 1:length(sessions_to_use)
            sessn = sessions_to_use(i);
            switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
                case {'Pre'}
                    t = Monkeys(m).Sessions(sessn).(slope_id);
                    test_bounds = t([1 3]);
                    post_bounds = post_data([1 3]);
                    % look for overlap in the confidence intervals.
                    if isempty(intersect(round(test_bounds(1),1):0.1:round(test_bounds(2),1), round(post_bounds(1),1):0.1:round(post_bounds(2),1)))
                        scatter(mean([pre_i_val post_i_val]), t(1)*1.05, '*k', 'HandleVisibility', 'off')
                        plot([pre_i_val post_i_val], repelem(t(1)*1.075, 1, 2) , '-k', 'HandleVisibility', 'off')
                    end
                case 'Post'
                    % dont test against itself
                    continue
            end
        end
    end
    
    % format plot
    xlabel('Session')
    xticks(1:length(sessions_to_use))
    xticklabs = Monkeys(m).XTickLabs;
    xticklabels(get_xtick_labs_colored(xticklabs, short_names, epoch_colors))
    xtickangle(45)
    ylabel('Slope (degrees)')
    ylim([40 51])
    yticks(40:5:50)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
        'fontsize', 26, ...
        'fontname', 'Helvetica', 'fontweight', 'normal', ...
        'XColor', 'black', 'YColor', 'black')
    % legend({'Monkey R', 'Monkey X'}, 'Interpreter', 'none')
    set(gcf, 'Color', 'w')
end
% sgtitle(sprintf('%s, slopes of major-axis regression, VR alpha %0.2f', area_to_plot, vr_alpha), 'Interpreter', 'none')

%% Extra: Compare intercepts of MAR regressions across days
vr_alpha = 0.05;
alpha_string_scale_factor = 100;
test_intervals = {[75 175], [175 350]};
for m = 1:length(Monkeys)
    % for m = 1
    
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
        sessions_to_use = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
        sessions_to_use = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
    end
    areas_to_plot = unique({Monkeys(m).Sessions(sessions_to_use(1)).UnitInfo.Area});
    
    for a = 1:length(areas_to_plot)
        figure('Position', [400 400 600 400]);
        hold on
        
        for p = 1:length(test_intervals)
            test_int = test_intervals{p};
            
            % pre allocate
            intercept_means = zeros(1,length(sessions_to_use));
            intercept_bounds = zeros(2,length(sessions_to_use));
            
            for i = 1:length(sessions_to_use)
                sessn = sessions_to_use(i);
                
                % use residuals from major axis regression
                int_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Intercept_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                t = Monkeys(m).Sessions(sessn).(int_id);
                intercept_means(i) = t(2);
                intercept_bounds(:,i) = [t(1); t(3)];
            end
            errorbar(intercept_means, diff(intercept_bounds)/2, 'DisplayName', sprintf('%d to %d', test_int(1), test_int(2)))
        end
        xticks(1:length(intercept_means))
        xticklabels({Monkeys(m).Sessions(sessions_to_use).ShortName})
        xtickangle(90)
        xlabel('Session')
        ylabel('Intercept (propn cat VR)')
        set(gca, 'FontSize', 20)
        legend
        sgtitle(sprintf('%s, %s, intercept of major-axis regression, VR alpha %0.2f', Monkeys(m).Name, areas_to_plot{a}, vr_alpha), 'Interpreter', 'none')
    end
end

%% (Spike count) Use MAR for same analysis but with avg cat/dog spike counts

mar_alpha = 0.05;
test_intervals = {[75 175], [175 350]};
make_plots = true;

for m = 1:length(Monkeys)
    % for m = 1
    sessions_to_use = 1:length(Monkeys(m).Sessions);
    %     if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
    % %         sessions_to_use = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
    %         sessions_to_use = [1 2 7 9];
    %     elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
    % %         sessions_to_use = [1 2 3 4 5 6 7 9]; % skip base04 and sub03 (low trial count)
    %         sessions_to_use = [1 2 6 7];
    %     end
    
    for p = 1:length(test_intervals)
        interval_to_plot = test_intervals{p};
        
        if make_plots
            figure('Position', [400 400 1300 1000])
            hold on
        end
        
        for i = 1:length(sessions_to_use)
            sessn = sessions_to_use(i);
            
            % only take te data
            arrays = {Monkeys(m).Sessions(sessn).UnitInfo.Location};
            areas = {Monkeys(m).Sessions(sessn).UnitInfo.Area};
            area_bool = strcmp(areas, 'te');
            
            % get data
            xid = get_good_interval_name2(interval_to_plot, 'full', 'X');
            X = Monkeys(m).Sessions(sessn).(xid);
            Y = Monkeys(m).Sessions(sessn).Session_Y_catg;
            cat_means = mean(X(Y==1, area_bool));
            dog_means = mean(X(Y==2, area_bool));
            
            % get angle of best-fit line with major axis regression
            % method. I checked this function against a reference, it's
            % correct.
            [coeffs,coeff_ints,~,~] = maregress(dog_means, cat_means, mar_alpha);
            theta = rad2deg(atan(coeffs(2)));
            theta_bounds = rad2deg(atan(coeff_ints(2,:)));
            intercept = coeffs(1);
            intercept_bounds = coeff_ints(1,:);
            % calculate distance from points to regression line:
            % D = abs(ax0+by0+c) / sqrt(a^2+b^2), where point is
            % (x0,y0), and line is ax+by+c = 0.
            % We have y=mx+c form, so a=coeffs(2), b = -1, and c =
            % coeffs(1).
            perpendicular_residuals = abs(coeffs(2)*dog_means + -1*cat_means + coeffs(1)) / sqrt(coeffs(2)^2 + 1);
            
            % store data
            propn_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('SpkCtMAR_CatDogPropns_%s', 'te'));
            slope_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('SpkCtMAR_Slope_%s', 'te'));
            intercept_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('SpkCtMAR_Intercept_%s', 'te'));
            resid_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('SpkCtMAR_Resids_%s', 'te'));
            
            Monkeys(m).Sessions(sessn).(propn_id) = [cat_means; dog_means]; % 2 x units
            Monkeys(m).Sessions(sessn).(resid_id) = perpendicular_residuals; % 1 x units
            Monkeys(m).Sessions(sessn).(slope_id) = [theta_bounds(1) theta theta_bounds(2)]; % 1 x 3
            Monkeys(m).Sessions(sessn).(intercept_id) = [intercept_bounds(1) intercept intercept_bounds(2)]; % 1 x 3
            
            % plot
            if make_plots
                subplot(3,4,i)
                hold on
                scatter(dog_means, cat_means, 'bo')
                xvals = 0:0.01:1;
                plot(xvals, xvals, 'k--')
                plot(xvals, coeffs(1) + xvals*coeffs(2), 'r--')
                xlabel('Avg dog spike count')
                ylabel('Avg cat spike count')
                xlim([0 3])
                ylim([0 3])
                title(sprintf('Session %s, theta = %0.1f (%0.1f - %0.1f)', Monkeys(m).Sessions(sessn).ShortName, theta, theta_bounds(1), theta_bounds(2)))
            end
        end
        if make_plots
            sgtitle(sprintf('%s, area %s, %d to %d, avg spike counts', Monkeys(m).Name, 'te', interval_to_plot(1), interval_to_plot(2)), 'Interpreter', 'none')
            
        end
    end
end

%% (Spike count) Succinctly plot selected MAR results with scatterplots

% interval_to_plot = [75 175];
interval_to_plot = [175 350];
area_to_plot = 'te';
titles = {'Base', 'Base', 'Post'};
title_cols = [[0.4 0.4 0.4];
    [0.4  0.4   0.4];
    [0   1   0]];
% figure('Position', [400 400 1600 1000])
figure('Position', [400 400 760 430])
hold on

for m = 1:length(Monkeys)
    
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
        %         sessions_to_use = [1 2 3 5 6 7 9]; % skip base04 (outlier) and post 01 (low trial count)
        sessions_to_use = [1 2 5 6 7 9]; % skip base04 (outlier) and post 01 (low trial count)
        %         sessions_to_use = [1 2 7 9];
        %         sessions_to_use = [1 7 9];
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
        sessions_to_use = [1 2 3 5 6 7]; % skip base04 due to all zeros on scatter plot (low trial count)
        %         sessions_to_use = [1 2 6 7];
        %         sessions_to_use = [1 6 7];
    end
    
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        
        % get data
        propn_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('SpkCtMAR_CatDogPropns_%s', 'te'));
        slope_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('SpkCtMAR_Slope_%s', 'te'));
        intercept_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('SpkCtMAR_Intercept_%s', 'te'));
        t = Monkeys(m).Sessions(sessn).(propn_id);
        cat_means = t(1,:);
        dog_means = t(2,:);
        slope = tan(deg2rad(Monkeys(m).Sessions(sessn).(slope_id)(2)));
        intercept = Monkeys(m).Sessions(sessn).(intercept_id)(2);
        xvals = 0:0.01:1;
        
        % plot data
        subplot(length(Monkeys),length(sessions_to_use),...
            length(sessions_to_use)*(m-1) + mod(i-1, length(sessions_to_use)) + 1)
        %         subplot(1, length(sessions_to_use), i)
        hold on
        scatter(dog_means, cat_means, 'bo')
        plot(xvals, xvals, 'k--', 'LineWidth', 2)
        plot(xvals, intercept + xvals*slope, 'r--', 'LineWidth', 2)
        
        % format plots
        if i == 1 && m == 1
            %         if i == 1
            xlabel('Avg dog spike count')
            ylabel('Avg cat spike count')
        end
        xlim([0 5])
        ylim([0 5])
        xticks([0 2.5 5])
        yticks([0 2.5 5])
        title(sprintf('%s', Monkeys(m).Sessions(sessn).ShortName), 'FontWeight', 'normal')
        %         title(titles{i}, 'Color', title_cols(i,:), 'FontWeight', 'normal')
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',18, 'YGrid', 'on', 'XGrid', 'on', ...
            'XColor', 'black', 'YColor', 'black')
    end
    sgtitle(sprintf('Area %s, %d to %d, avg spike count',...
        area_to_plot, interval_to_plot(1), interval_to_plot(2)), 'Interpreter', 'none')
    set(gcf, 'Color', 'w')
end

%% (Spike count) Statisically and visually compare selectivity index or MAR residuals across days with CDFs

cdf_alpha = 0.05;
ranksum_alpha = 0.025;
% test_interval = [175 350];
interval_to_plot = [75 175];
area_to_plot = 'te';

figure
hold on

for m = 1:length(Monkeys)
    subplot(2,1,m)
    hold on
    
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
        sessions_to_use = [1 2 3 5 6 7 9]; % skip base04 (outlier) and post 01 (low trial count)
        %         sessions_to_use = [1 2 7 9];
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
        sessions_to_use = [1 2 3 5 6 7]; % skip base04 due to all zeros on scatter plot (low trial count)
        %         sessions_to_use = [1 2 6 7];
    end
    
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        
        % use selectivity index (cats - dogs)
        %                 diff_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropDiffs_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
        %                 diffs = abs(Monkeys(m).Sessions(sessn).(diff_id));
        %                 [f, x] = ecdf(diffs);
        %                 xlabel_str = 'Selectivity (|prop_cat - prop_dog|)';
        
        % use residuals from major axis regression
        resid_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('SpkCtMAR_Resids_%s', 'te'));
        perpendicular_residuals = Monkeys(m).Sessions(sessn).(resid_id);
        [f, x] = ecdf(perpendicular_residuals);
        xlabel_str = 'Residuals from major axis regression';
        
        switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
            case 'Base'
                col = [0.4 0.4 0.4];
                thk = 1;
            case 'Pre'
                col = 'r';
                %                         col = [0.4 0.4 0.4];
                thk = 2;
                pre_data = perpendicular_residuals;
            case 'Post'
                col = 'g';
                thk = 2;
                post_data = perpendicular_residuals;
            case 'Sub'
                col = [0.8 0.8 0.8];
                thk = 1;
        end
        plot(x, f, 'Color', col, 'LineWidth', thk, 'DisplayName', Monkeys(m).Sessions(sessn).ShortName)
    end
    
    [pval,h] = ranksum(pre_data, post_data, 'Alpha', ranksum_alpha, 'Tail', 'left');
    fprintf('%s, %s, %d to %d, pre vs. post ranksum (one-sided non-parametric), hyp %d, pval of %d \n', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), h, pval)
    
    [h, pval] = ttest2(pre_data, post_data, 'tail', 'left');
    fprintf('%s, %s, %d to %d, pre vs. post ttest2 (one-sided), hyp %d, pval of %d \n', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), h, pval)
    
    %             title(sprintf('%s, %s, %d to %d, VR alpha %0.2f', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), vr_alpha), 'Interpreter', 'none')
    %             legend
    xlabel('Residuals')
    xlim([0 0.2])
    ylabel('Empirical CDF')
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',18, ...
        'fontname', 'Helvetica', 'fontweight', 'normal', ...
        'XColor', 'black', 'YColor', 'black')
    
end
set(gcf, 'Color', 'w')

%% (Spike count) Compare slopes of MAR regressions across days
test_interval = [175 350];
area_to_plot = 'te';

% make fig
figure('Position', [400 400 600 400]);
hold on
for m = 1:length(Monkeys)
    
    % select sessions
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
        sessions_to_use = [1 2 3 5 6 7 9]; % skip base04 (outlier) and post 01 (low trial count)
        %         sessions_to_use = [1 2 7 9];
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
        sessions_to_use = [1 2 3 4 5 6 7]; % skip base04 due to all zeros on scatter plot (low trial count)
        %         sessions_to_use = [1 2 6 7];
    end
    
    slope_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('SpkCtMAR_Slope_%s', 'te'));

    % pre allocate
    theta_means = zeros(1,length(sessions_to_use));
    theta_bounds = zeros(2,length(sessions_to_use));

    % collect data
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        t = Monkeys(m).Sessions(sessn).(slope_id);
        theta_means(i) = t(2);
        theta_bounds(:,i) = [t(1); t(3)];
        if ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post', 'once'))
            post_data = t;
        end
    end

    % plot data
    errorbar(theta_means, diff(theta_bounds)/2, 'LineWidth', 2, ...
        'DisplayName', sprintf('%s, %d to %d', Monkeys(m).Name, test_int(1), test_int(2)))

    % stats testing
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
            case {'Base', 'Pre'}
                t = Monkeys(m).Sessions(sessn).(slope_id);
                test_bounds = t([1 3]);
                post_bounds = post_data([1 3]);
                % look for overlap in the confidence intervals.
                if isempty(intersect(round(test_bounds(1),1):0.1:round(test_bounds(2),1), round(post_bounds(1),1):0.1:round(post_bounds(2),1)))
                    scatter(i*1.02, t(2)*1.02, 'vk', 'filled', 'HandleVisibility', 'off')
                end
            case 'Post'
                % dont test against itself
                continue
        end
    end
end
% format plot
% xticks(1:length(theta_means))
% xticklabels({Monkeys(m).Sessions(sessions_to_use).ShortName})
% xtickangle(90)
xlabel('Session')
% mark pre/post
plot([6.5 6.5], [40 52], 'k--', 'LineWidth', 1.5)
% xticklabels({'Base', 'Base', 'Pre', 'Post'})
xticklabels({'Base', 'Base', 'Base', 'Base', 'Base', 'Pre', 'Post'})
ylabel('Slope (degrees)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',18, ...
    'fontname', 'Helvetica', 'fontweight', 'normal', ...
    'XColor', 'black', 'YColor', 'black')
% legend({'Monkey R', 'Monkey X'}, 'Interpreter', 'none')
% sgtitle(sprintf('%s, slopes of major-axis regression, VR alpha %0.2f', area_to_plot, vr_alpha), 'Interpreter', 'none')
set(gcf, 'Color', 'w')

%% (NOT USEFUL) (Prop VR) Try different transforms to plot data
% Units lying on y=x indicate general visual responsiveness, while units
% out on the 90 deg lines would be extremely specifically dog- or
% cat-responsive.

% Conclusion: most units lie on y=x, indicating broad VR. TEO has way more
% units that respond to large proportion of imgs (ie units in the top
% right).

alpha = 0.05;
test_intervals = {[175 350]};
baseline_intervals = {[-175 0]};
model_scale = 'logit'; % 'logit' or 'sqrt' or 'lin' or 'logit'
for m = 1:length(Monkeys)
    sessions_to_use = 1:length(Monkeys(m).Sessions);
    for p = 1:length(test_intervals)
        test_int = test_intervals{p};
        baseline_int = baseline_intervals{p};
        tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img');
        bid = get_good_interval_name2(baseline_int, '', '');
        vr_id = strcat(tid,bid);
        
        areas_to_plot = unique({Monkeys(m).Sessions(sessions_to_use(1)).UnitInfo.Area});
        for a = 1:length(areas_to_plot)
            
            f1 = figure('Position', [400 400 1300 1000]);
            hold on
            
            fqq = figure('Position', [400 400 1300 1000]);
            hold on
            
            for i = 1:length(sessions_to_use)
                sessn = sessions_to_use(i);
                
                % get data
                data_mat = Monkeys(m).Sessions(sessn).(vr_id); % stimuli x units
                area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, areas_to_plot{a});
                data_mat = data_mat(:, area_bool);
                cat_props = sum(data_mat(1:260,:) < alpha) / 260;
                dog_props = sum(data_mat(261:520,:) < alpha) / 260;
                
                % store differences as measure of spread around y=x
                diff_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropDiffs_%s', areas_to_plot{a}));
                Monkeys(m).Sessions(sessn).(diff_id) = cat_props - dog_props;
                
                
                % run linear model and plot results
                if strcmp(model_scale, 'logit')
                    dog_props(dog_props ==0) = 0.5/260;
                    cat_props(cat_props ==0) = 0.5/260;
                    dog_props = log(dog_props./(1 - dog_props)); % min val is log(1/260) = -5.6.
                    cat_props = log(cat_props./(1-cat_props));
                    % plot
                    figure(f1)
                    subplot(3,4,i)
                    hold on
                    xvals = -6.5:0.5:6.5;
                    scatter(dog_props, cat_props, 'bo')
                    plot(xvals, xvals, 'k--')
                    xlabel('Log-Proportion dogs VR')
                    ylabel('Log-Proportion cats VR')
                    % run model
                    mdl = fitlm(dog_props, cat_props);
                    plot(xvals, xvals*mdl.Coefficients{'x1', 'Estimate'} + mdl.Coefficients{'(Intercept)', 'Estimate'}, 'r--')
                    title(sprintf('Session %s, R2 = %0.2f', Monkeys(m).Sessions(sessn).ShortName, mdl.Rsquared.Ordinary))
                    xlim([-6.5 6.5])
                    ylim([-6.5 6.5])
                    % qq
                    figure(fqq)
                    subplot(3,4,i)
                    hold on
                    qqplot(mdl.Residuals.Raw)
                    title(sprintf('Session %s, R2 = %0.2f', Monkeys(m).Sessions(sessn).ShortName, mdl.Rsquared.Ordinary))
                elseif strcmp(model_scale, 'sqrt')
                    dog_props = sqrt(dog_props); % min val is log(1/260) = -5.6.
                    cat_props = sqrt(cat_props);
                    % plot
                    figure(f1)
                    subplot(3,4,i)
                    hold on
                    xvals =0:0.01:1;
                    scatter(dog_props, cat_props, 'bo')
                    plot(xvals, xvals, 'k--')
                    xlabel('Sqrt-Proportion dogs VR')
                    ylabel('Sqrt-Proportion cats VR')
                    % run model
                    mdl = fitlm(dog_props, cat_props);
                    plot(xvals, xvals*mdl.Coefficients{'x1', 'Estimate'} + mdl.Coefficients{'(Intercept)', 'Estimate'}, 'r--')
                    title(sprintf('Session %s, R2 = %0.2f', Monkeys(m).Sessions(sessn).ShortName, mdl.Rsquared.Ordinary))
                    xlim([0 1])
                    ylim([0 1])
                    % qq plot
                    figure(fqq)
                    subplot(3,4,i)
                    hold on
                    qqplot(mdl.Residuals.Raw)
                    title(sprintf('Session %s, R2 = %0.2f', Monkeys(m).Sessions(sessn).ShortName, mdl.Rsquared.Ordinary))
                else
                    figure(f1)
                    subplot(3,4,i)
                    hold on
                    scatter(dog_props, cat_props, 'bo')
                    plot(0:0.01:1, -0:0.01:1, 'k--')
                    xlabel('Proportion dogs VR')
                    ylabel('Proportion cats VR')
                    % run model
                    mdl = fitlm(dog_props, cat_props);
                    xvals = 0:0.01:1;
                    plot(xvals, xvals*mdl.Coefficients{'x1', 'Estimate'} + mdl.Coefficients{'(Intercept)', 'Estimate'}, 'r--')
                    title(sprintf('Session %s, R2 = %0.2f', Monkeys(m).Sessions(sessn).ShortName, mdl.Rsquared.Ordinary))
                    xlim([0 1])
                    ylim([0 1])
                    % qq plot
                    figure(fqq)
                    subplot(3,4,i)
                    hold on
                    qqplot(mdl.Residuals.Raw)
                    title(sprintf('Session %s, R2 = %0.2f', Monkeys(m).Sessions(sessn).ShortName, mdl.Rsquared.Ordinary))
                end
            end
            figure(f1)
            sgtitle(sprintf('%s, area %s, image VR propns, %s fit (alpha = %0.2f)', Monkeys(m).Name, areas_to_plot{a}, model_scale, alpha), 'Interpreter', 'none')
            figure(fqq)
            sgtitle(sprintf('%s, area %s, qq plot for %s fit (alpha = %0.2f)', Monkeys(m).Name, areas_to_plot{a}, model_scale, alpha), 'Interpreter', 'none')
        end
    end
end

%% ******* Cat/dog difference timecourses (bootstrapped KDEs) *******

%% Re load data, add Area labels, get short names of sessions
clearvars
close all
bw = 20;
% path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
path1 = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = sprintf('MaxMarta_fix_cat_xma2_bootstrappedKDEs_BW%d.mat', bw);
Data = load(fullfile(path1, fname));
[status, KDE] = stitch_monkeyStruct_from_parts(Data);
clear Data

% add Area labels
for m = 1:length(KDE)
    for i = 1:length(KDE(m).Sessions)
        for j = 1:length(KDE(m).Sessions(i).UnitInfo)
            switch KDE(m).Sessions(i).UnitInfo(j).Location
                case 'posterior'
                    KDE(m).Sessions(i).UnitInfo(j).Area = 'te';
                case 'middle'
                    KDE(m).Sessions(i).UnitInfo(j).Area = 'te';
                case 'anterior'
                    KDE(m).Sessions(i).UnitInfo(j).Area = 'te';
                case 'teo'
                    KDE(m).Sessions(i).UnitInfo(j).Area = 'teo';
            end
        end
    end
end

% Get short names
marta_xls = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto//RecordingMarta.xlsx';
max_xls = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto//RecordingMax.xlsx';

for m = 1:length(KDE)
    
    % get short session names
    date_strs = {KDE(m).Sessions.DateStr};
    if regexp(KDE(m).Name, 'Marta\w*')
        short_names = get_short_names(date_strs, marta_xls, '\w*cats\w*');
    else
        short_names = get_short_names(date_strs, max_xls, '\w*cats\w*');
    end
    for i = 1:length(KDE(m).Sessions)
        KDE(m).Sessions(i).ShortName = short_names{i};
    end
end

%% Set sessions to use and xtick names and colors
xtickcolors = cbrewer('qual', 'Dark2', 3);
xtickcolors = xtickcolors([3 2 1], :);
for m = 1:length(KDE)
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
%         KDE(m).Sessions_to_use = [1 2 3 5 6 7 9];
%         KDE(m).XTickLabs= {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
        KDE(m).Sessions_to_use = [1 2 7 9];
        KDE(m).XTickLabs= {'Base 1', 'Base 2', 'Pre', 'Post'};
        KDE(m).Code = 'R';
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
        KDE(m).Sessions_to_use = [1 2 6 7];
        KDE(m).XTickLabs = {'Base 1', 'Base 2', 'Pre', 'Post'};
        KDE(m).Code = 'X';
    end
end

%% Get xvals where true diffs are larger than bootstrapped popn (fixed boundary) 

pct = 99;
diffid = sprintf('TrueDiff_BW%d', bw);
bootid = sprintf('BootstrappedDiffs_BW%d', bw);

for m = 1:length(KDE)
% for m = 2
   if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_plot = [1 2 7 9]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_plot = [1:3 5:9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_plot = [1 2 6 7]; % This is base 01, base02, pre and post. skip base03,04,05 and sub03 (low trial count)
        sessions_to_plot = 1:7;
   end
    
    KDE(m).Name = KDE(m).Name;
%     sessions_to_plot = 6;
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        %         units_to_plot = 1:length(KDE(m).Sessions(sessn).UnitInfo);
        units_to_plot = find(ismember({KDE(m).Sessions(sessn).UnitInfo.Location}, {'anterior', 'middle', 'posterior'}));
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals;
            true_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid);
            boot_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid);
            boot_diffs = abs(boot_diffs); % just look at magnitude, not sign
            
            % reshape to prepare for fixed boundary calculation
            boot_diffs = reshape(boot_diffs, 1, numel(boot_diffs));
            
            percentile_diff = prctile(boot_diffs, pct);
            exceedid = 'ExceedBD_FixedBound';
            times = find(abs(true_diffs) > percentile_diff); % normalized times in ms where true diffs exceeds bootstrapped diffs
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid) = times;
            
            if numel(times) > 20 && sum(times > 200)>20
                fprintf('Monkey %s, session %d, unit %d \n', KDE(m).Code, sessn, unit)
            end
            
        end
    end
end

%% Plot fixed boundary example

m = 1;
sessn = 7;
unit = 168; % m1, s1, 168 10 120 102 
bw = 20;
pct = 99;


% m1, s7:  2 (whole time), 32 (later), 58 (earlier), 71 (earlier), 72, 128?, 


fprintf('Monkey %s, session %d, unit %d, loc %s \n', KDE(m).Code, sessn, unit, KDE(m).Sessions(sessn).UnitInfo(unit).Location)

srs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_SRs;
f1 = 1000 * srs(2) / KDE(m).Sessions(sessn).NumCatgTrials(1) * KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg1_raw;
f2 = 1000 * srs(3) / KDE(m).Sessions(sessn).NumCatgTrials(2) * KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg2_raw;

diffid = sprintf('TrueDiff_BW%d', bw);
bootid = sprintf('BootstrappedDiffs_BW%d', bw);
true_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid);
boot_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid);

boot_diffs_reshape = reshape(boot_diffs, 1, numel(boot_diffs));
percentile_diff = prctile(abs(boot_diffs_reshape), pct); % 1 x (num xvals) vector of the 95th percentile of bootstrapped data
t = find(abs(true_diffs) > percentile_diff);

figure
subplot(2,1,1)
hold on
plot(kde_x_vals, f1, 'r-', 'LineWidth', 2)
plot(kde_x_vals, f2, 'b-', 'LineWidth', 2)
% xlabel('Time from cue on (ms)')
% ylabel('Spikes / second / trial')
xticks({})
xlim([-200 500])
yticks([0 4])
ylim([0 4])
% yticks([0 10])
% ylim([0 15])
% title(sprintf('%s, session %s, unit %d', KDE(m).Name, KDE(m).Sessions(sessn).ShortName, unit), 'Interpreter', 'none')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',18, ...
    'fontname', 'Helvetica', 'fontweight', 'normal', ...
    'XColor', 'black', 'YColor', 'black')
% make_custom_patch_legend({'r', 'b'}, {'Cats', 'Dogs'}, 'FontSize', 16, 'Location', 'northwest')

subplot(2,1,2)
hold on
plot(kde_x_vals, abs(true_diffs), 'g-', 'LineWidth', 2)
plot(kde_x_vals, abs(boot_diffs)', '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1)
plot([-200 500], repelem(percentile_diff, 2,1), 'r--', 'LineWidth', 1.5)
yl = get(gca, 'YLim');
% ylim([yl(1) 1.1*yl(2)])
ylim([0 1.5])
scatter(kde_x_vals(t), repelem(yl(2)*1.05, numel(t)), 'go', 'filled')
xlabel('Time from cue on (ms)')
xlim([-200 500])
xticks(kde_x_vals(1):200:kde_x_vals(end))
ylabel('Spikes / second / trial')
yticks([0 1.5])
% title(sprintf('%s, session %s, unit %d', KDE(m).Name, KDE(m).Sessions(sessn).ShortName, unit), 'Interpreter', 'none')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',18, ...
    'fontname', 'Helvetica', 'fontweight', 'normal', ...
    'XColor', 'black', 'YColor', 'black')
% make_custom_patch_legend({'g', [0.4 0.4 0.4], 'r'}, {'Abs. Diff', 'Abs. Shuff. Diff', '99th percentile (static)'}, 'FontSize', 16, 'Location', 'northwest')

sgtitle(sprintf('Monkey %s, %s, unit %d (%s)', ...
    KDE(m).Code, KDE(m).Sessions(sessn).ShortName, unit, ...
    KDE(m).Sessions(sessn).UnitInfo(unit).Location),...
    'FontSize', 18)
set(gcf, 'Color', 'white')

%% Look at many example units
m = 2;
sessn = 7;
units = [5 10 12 13 14 15 20 25]; % Max sessn 6: 187, 300, 317, 
% Max sessn 7: 10
% Good examples for martha sessn 7: 2, 6, 24 (later), 58, 90 (later), 112
bw = 20; 
pct = 99;
figure
for j = 1:length(units)
    unit = units(j);
    srs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_SRs;
    f1 = 1000 * srs(2) / KDE(m).Sessions(sessn).NumCatgTrials(1) * KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg1_raw;
    f2 = 1000 * srs(3) / KDE(m).Sessions(sessn).NumCatgTrials(2) * KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg2_raw;

    diffid = sprintf('TrueDiff_BW%d', bw);
    bootid = sprintf('BootstrappedDiffs_BW%d', bw);
    true_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid);
    boot_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid);

    boot_diffs_reshape = reshape(boot_diffs, 1, numel(boot_diffs));
    percentile_diff = prctile(abs(boot_diffs_reshape), pct); % 1 x (num xvals) vector of the 95th percentile of bootstrapped data
    t = find(abs(true_diffs) > percentile_diff);

    subplot(length(units),1,j)
    hold on
    plot(kde_x_vals, f1, 'r-', 'LineWidth', 2)
    plot(kde_x_vals, f2, 'b-', 'LineWidth', 2)
    % xlabel('Time from cue on (ms)')
    % ylabel('Spikes / second / trial')
    xticks({})
    xlim([-200 500])
    title(sprintf('%s, session %s, unit %d (%s)', KDE(m).Code, KDE(m).Sessions(sessn).ShortName, unit, KDE(m).Sessions(sessn).UnitInfo(unit).Location), 'Interpreter', 'none')
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',18, ...
        'fontname', 'Helvetica', 'fontweight', 'normal', ...
        'XColor', 'black', 'YColor', 'black')
end

%% Plot results as heatmap

area_to_plot = 'te';
boundary = 'static'; % static or dynamic
nrows_insert = 3; % size of red line separating arrays
sort_by = 'after_zero_latency'; % after_zero_latency or duration
sort_empty_val = 1000; % controls if units with no signf diff go above or below others

for m = 1:length(KDE)
    
%     if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
% %         sessions_to_plot = [1 2 7 9]; % skip base04 (outlier) and post01 (low trial count)
%         sessions_to_plot = [1:3 5:9];
%     elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
% %         sessions_to_plot = [1 2 6 7]; % This is base 01, base02, pre and post. skip base03,04,05 and sub03 (low trial count)
%         sessions_to_plot = 1:7;
%     end

    sessions_to_plot = KDE(m).Sessions_to_use;
    xticklabs = KDE(m).XTickLabs;
    
    
    figure('Position', [400 400 1500 680])
    hold on
    
    % get max num units
    num_units = zeros(1, length(sessions_to_plot));
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        num_units(i) = sum(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area_to_plot));
    end
    max_num_units = max(num_units);
    
    for i = 1:length(sessions_to_plot)
        subplot(1, length(sessions_to_plot), i)
        hold on
        sessn = sessions_to_plot(i);
        kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(1).CueOnAllCues.KDEXVals;
        
        % plotting color and label
        switch regexp(KDE(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
            case 'Base'
                col = xtickcolors(1,:);
            case 'Pre'
                col = xtickcolors(2,:);
            case 'Post'
                col = xtickcolors(3,:);
        end
        
        % take only te units
        area_bool = strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area_to_plot);
        units_to_plot = find(area_bool);
        
        % sort units by array
        array_list = {KDE(m).Sessions(sessn).UnitInfo(area_bool).Location};
        [array_list, idx] = sort(array_list);
        units_to_plot = units_to_plot(idx);
        
        % get data (already sorted by array)
        heatmap_mat = zeros(length(units_to_plot), length(kde_x_vals));
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            switch boundary
                case 'static'
                    exceedid = 'ExceedBD_FixedBound';
                    heatmap_mat(j, KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid)) = 1;
                case 'dynamic'
                    exceedid = 'ExceedBD_DynamicBound';
                    heatmap_mat(j, KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid)) = 1;
                otherwise
                    error('Unexpected value for boundary')
            end
        end
        
        % sort data in each array, put lines between arrays
        uq = unique(array_list);
        sorts = cell(3,1);
        for j = 1:length(uq)
            
            % find units in array
            unit_inds = find(strcmp(array_list, uq{j}));
            
            % sort
            full_data = heatmap_mat(unit_inds, :);
            to_sort = heatmap_mat(unit_inds, 201:end);
            switch sort_by
                case 'duration'
                    durs = zeros(1, size(to_sort,1));
                    for k = 1:size(to_sort,1)
                        s = sprintf('%d', to_sort(k,:));
                        s = textscan(s, '%s', 'delimiter', '0', 'multipleDelimsAsOne',1); % all runs of 1s
                        s = s{:}; % unpack
                        if isempty(s)
                            durs(k) = sort_empty_val;
                        else
                            durs(k) = max(cellfun('length', s)); % length of each run
                        end
                    end
                    [~,idx] = sort(durs);
                    
                case 'after_zero_latency'
                    first_one = zeros(1, size(to_sort,1));
                    for k = 1:size(to_sort,1)
                        %                         first_one_ind = find(to_sort(k,:)==1,1, 'last');
                        first_one_ind = find(to_sort(k,:)==1,1, 'first');
                        if isempty(first_one_ind)
                            first_one(k) = sort_empty_val;
                        else
                            first_one(k) = first_one_ind;
                            %                             if first_one_ind < 100
                            %                                 full_data(k,full_data(k,:)==1) = 0.66;
                            %                             end
                        end
                    end
                    [~,idx] = sort(first_one);
            end
            idx = flipud(idx);
            sorts{j} = idx;
            heatmap_mat(unit_inds, :) = full_data(idx,:);
            
            % find spot for the line
            array_end_ind = max(unit_inds);
            
            % add data for the lines
            heatmap_mat = insertrows(heatmap_mat, repmat(zeros(1, length(kde_x_vals))+0.5, nrows_insert,1), array_end_ind);
            array_list = insertrows(array_list', repelem({'NA'}, nrows_insert, 1), array_end_ind)';
        end
        
        % store array sorts
        sortid = sprintf('HeatmapSortInds_%s', sort_by);
        KDE(m).Sessions(sessn).(sortid) = sorts;
        
        % add yticklabels for arrays
        ytick_spots = zeros(1, length(uq)+1);
        ytick_labs = cell(1, length(uq)+1);
        for j = 1:length(uq)
            unit_inds = find(strcmp(array_list, uq{j}));
            array_end_ind = max(unit_inds);
            ytick_spots(j) = round(median(unit_inds));
            ytick_labs{j} = sprintf('%s. (%d)', strcat(upper(uq{j}(1)), uq{j}(2:3)), numel(unit_inds));
        end
        
        % set colormap
        % 0,    0.5,    1
        % none  spacer  diff
        cm = parula(2);
        cm = [cm(1,:); [1 0 0]; cm(2,:)];
        colormap(cm);
        
        % plot data
        imagesc(heatmap_mat)
        
        % add ytick at top
        ytick_spots(end) = max_num_units;
        ytick_labs{end} = num2str(max_num_units);
        
        % sort to be ascending so MATLAB is happy
        [ytick_spots, idx] = sort(ytick_spots);
        ytick_labs = ytick_labs(idx);
        
        % format plot
        title(sprintf('%s', KDE(m).Sessions(sessn).ShortName), 'Color', col)
        if i == 1
            xlabel('Time from cue on')
            ylabel('Array (units)')
        end
        xticks(0:200:700)
        xticklabels(-200:200:500)
        ylim([0 max_num_units])
        yticks(ytick_spots)
        yticklabels(ytick_labs)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',26, ...
            'fontname', 'Helvetica', 'fontweight', 'normal', ...
            'XColor', 'black', 'YColor', 'black')
    end
    set(gcf, 'Renderer', 'Painters');
    set(gcf, 'Color', 'white');
    saveas(gcf, sprintf('catdog_f5 heatmaps mk%s.png', KDE(m).Code))
end

%% Plot avg timecourses (cat/dog diff) with days overlaid (learning)

% params
areas_to_plot = {'anterior', 'middle', 'posterior'};
% areas_to_plot = {'te'};
boundary = 'static'; % static or dynamic
% Not convinced correction is needed in this case.
% chi_sq_base_alpha = 0.05;
% chi_sq_alpha = chi_sq_base_alpha / numel(chi_sq_inds);
chi_sq_alpha = 0.05;
colors = cbrewer('qual', 'Dark2', 3);
colors = colors([3 2 1], :);
range_normalize = false;


% pre allocate
pre_data = cell(length(KDE), length(areas_to_plot));
post_data = cell(length(KDE), length(areas_to_plot));
chi_sq_inds = 200:50:700; % inds in kde_x_vals to test...corresponds to -200:50:500 in real time.



for m = 1:length(KDE)
    
    % get sessions to use
    sessions_to_plot = KDE(m).Sessions_to_use;
%     if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_plot = [7 9];
%     elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_plot = [6 7];
%     end
    
    
    % prep figure
    figure('Position', [200 200 1200 330]) % for arrays
%     figure('Position', [400 400 500 300]) % just te
    hold on
    
    for a = 1:length(areas_to_plot)
        area = areas_to_plot{a};
        
                % Plot all in one figure:
%                 switch area
%                     case 'te'
%                         subplot(2,3,1:3)
%                     otherwise
%                         subplot(2,3,2+a)
%                 end
%                 hold on
                % TE alone: don't do any subplots
                
                % just arrays
                subplot(1, 3, a)
                hold on
        
        % get data and plot
        for i = 1:length(sessions_to_plot)
            sessn = sessions_to_plot(i);
            
            % get correct units
            switch area
                case 'te'
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area));
                otherwise
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, area));
            end
            
            % get color to plot with
            if ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Base', 'once'))
                col_ind = 1;
            elseif ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Pre', 'once'))
                col_ind = 2;
            elseif ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Post', 'once'))
                col_ind = 3;
            end
            
            % get real x-axis values
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(1).CueOnAllCues.KDEXVals;
            
            % get data
            heatmap_mat = zeros(length(units_to_plot), length(kde_x_vals));
            for j = 1:length(units_to_plot)
                unit = units_to_plot(j);
                switch boundary
                    case 'static'
                        exceedid = 'ExceedBD_FixedBound';
                        heatmap_mat(j, KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid)) = 1;
                    case 'dynamic'
                        exceedid = 'ExceedBD_DynamicBound';
                        heatmap_mat(j, KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid)) = 1;
                    otherwise
                        error('Unexpected value for boundary')
                end
            end
            
            % get mean and plot
            heatmap_mean = mean(heatmap_mat); % 1 x num xvals
            if range_normalize
                heatmap_mean = normalize(heatmap_mean, 'range');
            end
            plot(heatmap_mean,  'LineWidth', 1.5, ...
                'DisplayName', KDE(m).Sessions(sessn).ShortName,...
                'Color', colors(col_ind, :))
            
            % format plot
%             title(sprintf('Monkey %s, %s', KDE(m).Code, area), 'Interpreter', 'none')
            title(sprintf('%s', area), 'Interpreter', 'none')
            xticks(0:200:600)
            xticklabels(kde_x_vals(1:200:end))
            if strcmp(area, 'te')
                xlabel('Time from cue on')
                ylabel('Fraction units')
            end
            
            % store pre / post data for stats testing
            if ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Post', 'once'))
                post_data{m,a} = heatmap_mat;
            elseif ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Pre', 'once'))
                pre_data{m,a} = heatmap_mat;
            end
        end
        
        % stats testing
        pre = pre_data{m,a};
        post = post_data{m,a};
        yl = get(gca, 'YLim');
        for ii = chi_sq_inds
            [h,p] = prop_test([sum(pre(:,ii)) sum(post(:,ii))], [size(pre,1) size(post,1)], false, chi_sq_alpha);
            if h
                scatter(ii, yl(2)*0.97, 'ko', 'filled')
            end
        end
        
        if a == 1
            ylabel('Proportion signf. units')
            xlabel('Time from cue on')
        end
        
        % format plot
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',26, ...
            'fontname', 'Helvetica', ...
            'XColor', 'black', 'YColor', 'black')
        
    end
%     sgtitle(sprintf('Monkey %s', KDE(m).Code))
end

%% Load EXCIT / SUPPRN latency data, add Area labels to it

path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = 'MaxMarta_xma2_temporalActvnSuppn.mat';
Data = load(fullfile(path1, fname));
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

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

%% Plot avg timecourse (EXCIT) with days overlaid (learning control)

areas_to_plot = {'anterior', 'middle', 'posterior'};
boundary = 'static'; % static or dynamic
% array_colors = cbrewer('qual', 'Set1', 3);
colors = cbrewer('qual', 'Dark2', 3);
colors = colors([3 2 1], :);
excit_supprn_timecourse_xvals = -50:10:430;

for m = 1:length(KDE)
    
    figure('Position', [200 200 1200 330]) % for arrays
    hold on
    
    % get sessions to use
    sessions_to_plot = KDE(m).Sessions_to_use;
    for a = 1:length(areas_to_plot)
        area = areas_to_plot{a};
        subplot(1, 3, a)
        hold on 
        
        for i = 1:length(sessions_to_plot)
            sessn = sessions_to_plot(i);
            % get correct units
            switch area
                case 'te'
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area));
                otherwise
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, area));
            end
            
            % get color to plot with
            if ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Base', 'once'))
                col_ind = 1;
            elseif ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Pre', 'once'))
                col_ind = 2;
            elseif ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Post', 'once'))
                col_ind = 3;
            end
            
            % get data
            heatmap_mat = zeros(length(units_to_plot), length(excit_supprn_timecourse_xvals));
            for j = 1:length(units_to_plot)
                unit = units_to_plot(j);
                heatmap_mat(j, [Monkeys(m).Sessions(sessn).UnitInfo(unit).WindowResults]==2) = 1;
            end
            
            %             plot(normalize(mean(heatmap_mat), 'range'), 'Color', colors(a,:), ...
            %                 'DisplayName', areas_to_plot{a}, 'LineWidth', 1.5)
            plot(excit_supprn_timecourse_xvals, mean(heatmap_mat), 'Color', colors(col_ind,:), ...
                'DisplayName', areas_to_plot{a}, 'LineWidth', 1.5)
        end
        title(sprintf('%s', area), 'Interpreter', 'none')
        xlabel('Time from cue on')
        ylabel('Fraction units')
%         ylim([0 0.7])
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',18, ...
            'fontname', 'Helvetica',  ...
            'XColor', 'black', 'YColor', 'black')
    end
    sgtitle(sprintf('Monkey %s', KDE(m).Code))
end

%% Plot cat/dog diff normalized by EXCIT timecourses, days overlaid (learning)

% params
areas_to_plot = {'anterior', 'middle', 'posterior'};
% areas_to_plot = {'te'};
boundary = 'static'; % static or dynamic
% Not convinced correction is needed in this case.
% chi_sq_base_alpha = 0.05;
% chi_sq_alpha = chi_sq_base_alpha / numel(chi_sq_inds);
chi_sq_alpha = 0.05;
colors = cbrewer('qual', 'Dark2', 3);
colors = colors([3 2 1], :);
excit_supprn_timecourse_xvals = -50:10:430;

% pre allocate
pre_data = cell(length(KDE), length(areas_to_plot));
post_data = cell(length(KDE), length(areas_to_plot));


for m = 1:length(KDE)
    
    % get sessions to use
    sessions_to_plot = KDE(m).Sessions_to_use;
%     if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_plot = [7 9];
%     elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_plot = [6 7];
%     end
    
    
    % prep figure
    figure('Position', [200 200 1200 330]) % for arrays
%     figure('Position', [400 400 500 300]) % just te
    hold on
    
    for a = 1:length(areas_to_plot)
        area = areas_to_plot{a};
        
                % Plot all in one figure:
%                 switch area
%                     case 'te'
%                         subplot(2,3,1:3)
%                     otherwise
%                         subplot(2,3,2+a)
%                 end
%                 hold on
                % TE alone: don't do any subplots
                
                % just arrays
                subplot(1, 3, a)
                hold on
        
        % get data and plot
        for i = 1:length(sessions_to_plot)
            sessn = sessions_to_plot(i);
            
            % get correct units
            switch area
                case 'te'
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area));
                otherwise
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, area));
            end
            
            % get color to plot with
            if ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Base', 'once'))
                col_ind = 1;
            elseif ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Pre', 'once'))
                col_ind = 2;
            elseif ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Post', 'once'))
                col_ind = 3;
            end
            
            % get real x-axis values
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(1).CueOnAllCues.KDEXVals;
            
            % get cat/dog diff data
            heatmap_mat_catdogdiff = zeros(length(units_to_plot), length(kde_x_vals));
            for j = 1:length(units_to_plot)
                unit = units_to_plot(j);
                switch boundary
                    case 'static'
                        exceedid = 'ExceedBD_FixedBound';
                        heatmap_mat_catdogdiff(j, KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid)) = 1;
                    case 'dynamic'
                        exceedid = 'ExceedBD_DynamicBound';
                        heatmap_mat_catdogdiff(j, KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid)) = 1;
                    otherwise
                        error('Unexpected value for boundary')
                end
            end
            heatmap_mat_catdogdiff = mean(heatmap_mat_catdogdiff);
            
            % get EXIT data
            heatmap_mat_excit = zeros(length(units_to_plot), length(excit_supprn_timecourse_xvals));
            for j = 1:length(units_to_plot)
                unit = units_to_plot(j);
                heatmap_mat_excit(j, [Monkeys(m).Sessions(sessn).UnitInfo(unit).WindowResults]==2) = 1;
            end
            heatmap_mat_excit = mean(heatmap_mat_excit);
            
            % align vectors and normalize
            inds = ismember(kde_x_vals, excit_supprn_timecourse_xvals);
            heatmap_mat_normd = heatmap_mat_catdogdiff(inds) ./ heatmap_mat_excit;
            
            % get mean and plot
            plot(excit_supprn_timecourse_xvals, heatmap_mat_normd,  'LineWidth', 1.5, ...
                'DisplayName', KDE(m).Sessions(sessn).ShortName,...
                'Color', colors(col_ind, :))
            
            % format plot
%             title(sprintf('Monkey %s, %s', KDE(m).Code, area), 'Interpreter', 'none')
            title(sprintf('%s', area), 'Interpreter', 'none')
            xlabel('Time from cue on')
            ylabel({'Fraction catdog diff', 'out of Vis Excit. units'})
            
            
            % store pre / post data for stats testing
            if ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Post', 'once'))
                post_data{m,a} = heatmap_mat;
            elseif ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Pre', 'once'))
                pre_data{m,a} = heatmap_mat;
            end
        end
        
        % format plot
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',18, ...
            'fontname', 'Helvetica', ...
            'XColor', 'black', 'YColor', 'black')
        
    end
    sgtitle(sprintf('Monkey %s', KDE(m).Code))
end
