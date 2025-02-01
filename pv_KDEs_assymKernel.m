% looks for significant differences in cat/dog responses by bootstrapping.
% uses assym kernel so we can accurately assess latencies.


%% Load data
clearvars
close all

EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs'; % pv --> "passive viewing"


fname = "MaxMarta_xma2_ni.mat";
Data = load(fullfile(EXT_HD, pv_path, fname));
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% (run this every time) Parameters

%{
-- We'll go from 100 ms before cue on to 400 ms after cue on.
-- We'll use an assymetric gaussian kernel based on Connor lab (the og idea was actually
modeled after an EPSP but they seemt to have abandonded that, probably for
computational efficiency).

%}

sessions_by_monk = {[7 9], [6 7]};

bandwidth = 15;  % msec — and then this is applied in a 1:3 ratio by the kernel, see below
buffer = 2; % how many bw of buffer to give the KDE

% ms to use before and after cue
plot_window_back = 100;
plot_window_front = 400;

% downsample a bit to speed things up
kde_x_vals = -(plot_window_back):2:plot_window_front; % exclude bw buffers now

categories = {1:260, 261:520};

%% (run this every time) Add Area + get short session names

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

%% Transfer data about total num trials
KDE = struct();
for m = 1:length(Monkeys)
    sessions_to_use = sessions_by_monk{m};
    KDE(m).Name = Monkeys(m).Name;
    
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        KDE(m).Sessions(sessn).SessionName = Monkeys(m).Sessions(sessn).SessionName;
        KDE(m).Sessions(sessn).DateStr = Monkeys(m).Sessions(sessn).DateStr;
        num_total_cat_trials = sum([Monkeys(m).Sessions(sessn).CueInfo(categories{1}).NumApp]);
        num_total_dog_trials = sum([Monkeys(m).Sessions(sessn).CueInfo(categories{2}).NumApp]);
        KDE(m).Sessions(sessn).NumCatgTrials = [num_total_cat_trials num_total_dog_trials];
    end
end

%% Calculate mean KDEs for each neuron (takes a few min)

spike_bounds_back = plot_window_back + bandwidth*buffer;
spike_bounds_front = plot_window_front + bandwidth*buffer;


for m = 1:length(Monkeys)
    sessions_to_use = sessions_by_monk{m};
    
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        units_to_use = 1:length(Monkeys(m).Sessions(sessn).UnitInfo);
        
        for j = 1:length(units_to_use) % for each unit
            unit = units_to_use(j);
            if ~ismember(Monkeys(m).Sessions(sessn).UnitInfo(j).Location, ["anterior", "middle", "posterior"])
                KDE(m).Sessions(sessn).UnitInfo(j).Location = 'not_te';  % just store a place holder string
                continue
            end
            
            max_trials = max([Monkeys(m).Sessions(sessn).CueInfo.NumApp]);
            allspike_times = sort(Monkeys(m).Sessions(sessn).UnitInfo(unit).Spike_times); % sort once to save time with many binary searches later
            
            KDE(m).Sessions(sessn).UnitInfo(unit).Electrode = Monkeys(m).Sessions(sessn).UnitInfo(unit).Electrode;
            KDE(m).Sessions(sessn).UnitInfo(unit).Electrode_unit = Monkeys(m).Sessions(sessn).UnitInfo(unit).Electrode_unit;
            KDE(m).Sessions(sessn).UnitInfo(unit).Location = Monkeys(m).Sessions(sessn).UnitInfo(unit).Location;
            
            % preallocate rough size of arrays for plotting
            % say avg rate is 25 sp/s = 0.025 sp / ms
            avg_rate = 0.25*(spike_bounds_back + spike_bounds_front);
            x_vals = zeros(max_trials*260*avg_rate, 1); % (trials x cues) = row, row x num_spikes / row = num_spikes
            y_vals = zeros(max_trials*260*avg_rate, 1);
            catg_vals = zeros(max_trials*260*avg_rate, 1);
            xx = 1; % indexer
            
            % collect data
            for k = 1:length(Monkeys(m).Sessions(sessn).CueInfo) % for each cue
                
                % NB, this times_on field ONLY contains timestamps from
                % within completed trials (see behav_analysis code).
                times_on = Monkeys(m).Sessions(sessn).CueInfo(k).Times_on;
                for q = 1:2
                    if ismember(Monkeys(m).Sessions(sessn).CueInfo(k).CueID, categories{q})
                        catg = q;
                        break
                    end
                end
                    
                bounds = [(times_on - spike_bounds_back) (times_on + spike_bounds_front)];
                for l = 1:length(bounds) % for each time that cue appeared
                    [lower_idx, upper_idx] = binarySearch_window(allspike_times, bounds(l, 1), bounds(l, 2)); % allspike_times must be sorted!
                    
                    % check edge cases
                    if lower_idx == -1 || upper_idx == -1 || upper_idx - lower_idx == -1
                        % no spikes in window
                        % lower_idx is -1 --> bounds are higher than
                        % all spike times
                        % uppder idx is -1 --> bounds are lower than spike
                        % times
                        % diff is -1 --> bounds are in between two spike
                        % times
                        continue
                    end
                    trial_spike_times = allspike_times(lower_idx:upper_idx);
                    norm_times = trial_spike_times - times_on(l);
                    x_vals(xx:xx + length(norm_times) - 1) = norm_times;
                    y_vals(xx:xx + length(norm_times) - 1) = repelem(k + 0.01*l, length(norm_times));
                    catg_vals(xx:xx + length(norm_times) - 1) = repelem(catg, length(norm_times));
                    xx = xx + length(norm_times);
                end
            end
            
            % remove extra zeros
            x_vals(xx:end) = [];
            y_vals(xx:end) = [];
            catg_vals(xx:end) = [];
            
            if isempty(x_vals)
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterXVals = x_vals;
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterYVals = y_vals;
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals = kde_x_vals;
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals = NaN;
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.BW = bandwidth;
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.Buffer = buffer;
                continue
            end
            
            % run KDE
%             f = ksdensity(x_vals, kde_x_vals, 'Bandwidth', bw, 'Kernel', @asymmetric_gaussian_kernel);
            
            % run catg specific KDEs. use all points to create KDE
            if sum(catg_vals==1) == 0
                f_catg1 = zeros(1, numel(kde_x_vals));
                f_catg2 = ksdensity(x_vals(catg_vals==2), kde_x_vals, 'Kernel', @asymmetric_gaussian_kernel, 'Bandwidth', bandwidth);
            elseif sum(catg_vals==2) == 0
                f_catg1 = ksdensity(x_vals(catg_vals==1), kde_x_vals, 'Kernel', @asymmetric_gaussian_kernel, 'Bandwidth', bandwidth);
                f_catg2 = zeros(1, numel(kde_x_vals));
            else
                f_catg1 = ksdensity(x_vals(catg_vals==1), kde_x_vals,  'Kernel', @asymmetric_gaussian_kernel, 'Bandwidth', bandwidth);
                f_catg2 = ksdensity(x_vals(catg_vals==2), kde_x_vals,  'Kernel', @asymmetric_gaussian_kernel, 'Bandwidth', bandwidth);
            end


            % store
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterXVals = x_vals;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterYVals = y_vals;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.CatgVals = catg_vals;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals = kde_x_vals;
%             KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_raw = f;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg1_raw = f_catg1;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg2_raw = f_catg2;
%             KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_SRs = [xx sum(catg_vals==1) sum(catg_vals==2)]; % [total, cats, dogs]
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_SRs = [sum(catg_vals==1) sum(catg_vals==2)]; % [cats, dogs]
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.BW = bandwidth;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.Buffer = buffer;
            
            fprintf("Done with Monkey %s, session %d, unit %d \n", Monkeys(m).Name, sessn, unit)
            
        end
    end
end

%% Calculate shuffled KDEs (takes hours)

n_shuff = 2;  % incr this eventually on an HPCC
rng(10)
categories = {1:260, 261:520};

bootid = sprintf('BootstrappedDiffs_BW%d', bandwidth);

for m = 1:length(KDE)
    sessions_to_use = sessions_by_monk{m};
    
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        units_to_use = 1:length(Monkeys(m).Sessions(sessn).UnitInfo);
        
        for j = 1:length(units_to_use)
            unit = units_to_use(j);
            
            if ~ismember(Monkeys(m).Sessions(sessn).UnitInfo(j).Location, ["anterior", "middle", "posterior"])
                continue
            end
            
            
            if isempty(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterXVals)
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid) = NaN;
                continue
            end
            
            % get data
            raster_x_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterXVals; % normalized spike times
            raster_y_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterYVals; % useful for keeping trials together in the shuffle
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals; % normalized time
            catg_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.CatgVals; % cat or dog for each spike
            num_total_cat_trials = KDE(m).Sessions(sessn).NumCatgTrials(1);
            num_total_dog_trials = KDE(m).Sessions(sessn).NumCatgTrials(2);
            
            
            % catch units with very few spikes, or none in one category
            if sum(catg_vals==1) == 0 || sum(catg_vals==2) == 0
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid) = NaN;
                continue
            elseif numel(catg_vals) < 20
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid) = NaN;
                continue
            end
            
            % get trial inds, to keep trials intact
            uq = unique(raster_y_vals); % 1 x num trials, with unique trial IDs
            trial_inds = cell(1, length(uq)); % 1 x num trials, {[spike inds in tr1], [spike inds in tr2], etc}.
            for t = 1:length(uq)
                trial_inds{t} = find(raster_y_vals == uq(t));
            end
            
            % pre-allocate / pre-calculate some params
            bootstrapped_diffs = zeros(n_shuff, length(kde_x_vals));
            shuffled_catgs = catg_vals; % to be shuffled
            n_uq = length(uq);
            
            % bootstrap
            for n = 1:n_shuff
                
                % Shuffle categories while keeping trials intact AND
                % balance cats/dogs randomly about 50/50.
                num = numel(shuffled_catgs);
                while true
                    img_shuffle = (rand(1, n_uq) > 0.5)+1;
                    for t = 1:n_uq
                        shuffled_catgs(trial_inds{t}) = img_shuffle(t);
                    end
                    
                    % if there's a low-ish number of spikes, ensure the
                    % shuffle is close to 50/50.
                    % binocdf(200,500,0.5) = 5e-6, implying that if every
                    % unit had 500 spikes only, ~all of distributions would
                    % be within 0.1 of perfect 50/50 cat/dog. So we'll say
                    % less than 500 is low-ish, and above that, I expect
                    % never to have a problem with 50*250*10*2 = 2.5e6
                    % iterations.
                    if num > 5e2
                        break
                    else
                        n1 = sum(shuffled_catgs==1);
                        if n1/num > 0.4 && n1/num < 0.6
                            break
                        end
                    end
                end
                
                % fit KDEs
                shuffled_f_catg1 = ksdensity(raster_x_vals(shuffled_catgs==1), kde_x_vals, 'Bandwidth', bandwidth);
                shuffled_f_catg2 = ksdensity(raster_x_vals(shuffled_catgs==2), kde_x_vals, 'Bandwidth', bandwidth);
                
                % NEW WAY
                % scale KDEs like so: (density / ms) x (1000 ms / s) x (num spikes) / (num trials)  = spikes / second / trial.
                shuffled_f_catg1 = shuffled_f_catg1 * 1000 * sum(shuffled_catgs==1) / num_total_cat_trials;
                shuffled_f_catg2 = shuffled_f_catg2 * 1000 * sum(shuffled_catgs==2) / num_total_dog_trials;
                
                % get diff
                bootstrapped_diffs(n,:) = shuffled_f_catg1 - shuffled_f_catg2;
            end
            
            % store
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid) = bootstrapped_diffs;
            fprintf('Done with %s, session %d, unit %d \n', KDE(m).Name, sessn, unit)
            
            % Sanity check
%             figure
%             hold on
%             plot(bootstrapped_diffs', 'Color', [.7 .7 .7])
%             f1 = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg1_raw;
%             f2 = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg2_raw;
%             f1 = f1 * 1000 * sum(catg_vals==1) / num_total_cat_trials;
%             f2 = f2 * 1000 * sum(catg_vals==2) / num_total_dog_trials;
%             plot(f1 - f2, 'Color', 'r', 'LineWidth', 3)
%             disp('hi')
            
        end
    end 
end

%% Save data
fname = sprintf('MaxMarta_fix_cat_xma2_bootstrappedKDEs_ASYMM_BW%d_lowShuffcount.mat', bandwidth);
[status, vnames] = split_monkeyStruct_in_parts(KDE);
save(fullfile(EXT_HD, pv_path, fname), vnames{:});

%% (restart here after saving) Re-load data

EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs'; % pv --> "passive viewing"

bandwidth = 15;
fname = sprintf('MaxMarta_fix_cat_xma2_bootstrappedKDEs_ASYMM_BW%d_Marta_fix_cat_xma2', bandwidth);
% fname = sprintf('MaxMarta_fix_cat_xma2_bootstrappedKDEs_ASYMM_BW%d_Max_fix_cat_xma2', bandwidth);
Data = load(fullfile(EXT_HD, pv_path, fname));
[status, KDE] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Tmp hot fix
% just working with marta's data rn, delete Max
KDE = KDE(1);
Monkeys = Monkeys(1);
rSessionsByMonk = {[7 9]};

% Max data only
% KDE = KDE(2);
% Monkeys = Monkeys(2);
% rSessionsByMonk = {[6 7]};

%% Add area labels

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

%% Get xvals where true diffs are larger than bootstrapped popn (dynamic boundary) 

pct = 99;
diffid = sprintf('TrueDiff_BW%d', bandwidth);
bootid = sprintf('BootstrappedDiffs_BW%d', bandwidth);
% rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)
rArrayLocs = {'te'};

for m = 1:length(KDE)
    sessions_to_plot = rSessionsByMonk{m};
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
%         units_to_plot = 1:length(KDE(m).Sessions(sessn).UnitInfo);
        units_to_plot = find(ismember({KDE(m).Sessions(sessn).UnitInfo.Location}, {'anterior', 'middle', 'posterior'}));
        
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            
            if isfield(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues, 'KDEYVals') && ...
                isnan(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals)
                continue
            end
            
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals;
            catg_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.CatgVals;
            num_total_cat_trials = KDE(m).Sessions(sessn).NumCatgTrials(1);
            num_total_dog_trials = KDE(m).Sessions(sessn).NumCatgTrials(2);
            
            catg1_scaled = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg1_raw * 1000 * sum(catg_vals==1) / num_total_cat_trials;
            catg2_scaled = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg2_raw * 1000 * sum(catg_vals==2) / num_total_dog_trials;
            true_diffs = catg1_scaled - catg2_scaled;
            
%             true_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid);
            boot_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid);
            boot_diffs = abs(boot_diffs); % just look at magnitude, not sign
            percentile_diffs = prctile(boot_diffs, pct, 1); % 1 x (num xvals) vector of the 95th percentile of bootstrapped data
            exceedid = 'ExceedBD_DynamicBound';
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid) = find(abs(true_diffs) > percentile_diffs); % normalized times in ms where true diffs exceeds bootstrapped diffs
            
        end
    end
end

%% Plot results as heatmap (array starts aligned) (Fig 4B)

area_to_plot = 'te';
boundary = 'dynamic'; % static or dynamic
nrows_insert = 3; % size of black line separating arrays
sort_by = 'after_zero_latency'; % after_zero_latency or duration
sort_empty_val = 1000; % controls if units with no signf diff go above or below others
divider_color = [0.5 0.5 0.5 ];
save_dir = '/Users/jonahpearl/Documents/BJR group/Catdog_paper/Feb_2023_addtl_figs/kde_diff_heatmaps';

rArrayLocs = {'anterior', 'middle', 'posterior'};

for m = 1:length(KDE)

    sessions_to_plot = rSessionsByMonk{m};
    
%     figure2('Position', [400 400 1500 680]) % four cols
    figure2('Position', [400 400 600 680]) % four cols
    axes = [];
    hold on
    
    % Get max num units, so that all plots are the same height.
    num_units_by_array = zeros(length(sessions_to_plot), length(rArrayLocs));
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        for iLoc = 1:length(rArrayLocs)
            loc = rArrayLocs{iLoc};
            num_units_by_array(i, iLoc) = sum(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, loc));
        end
    end
    
    max_num_units_per_array = reshape(max(num_units_by_array, [], 1), 1, []);
    num_rows_inserted = 0;
    
    for i = 1:length(sessions_to_plot)
        ax = subplot(1, length(sessions_to_plot), i);
        axes = [axes ax];
        hold on
        sessn = sessions_to_plot(i);
        kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(1).CueOnAllCues.KDEXVals;
        
        % Plotting color and label
        switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
            case 'Base'
%                 col = xtickcolors(1,:);
                col = [0.5 0.5 0.5];
            case 'Pre'
%                 col = xtickcolors(2,:);
                col = mlc(1);
            case 'Post'
%                 col = xtickcolors(3,:);
                col = mlc(2);
        end
        
        % Take only te units
        area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, area_to_plot);
        units_to_plot = find(area_bool);
        
        % Sort units by array
        array_list = {Monkeys(m).Sessions(sessn).UnitInfo(area_bool).Location};
        [array_list, idx] = sort(array_list);
        units_to_plot = units_to_plot(idx);
        
        % Get data (already sorted by array)
        heatmap_mat = zeros(length(units_to_plot), length(kde_x_vals));
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            
            if isfield(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues, 'KDEYVals') && ...
                isnan(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals)
                continue
            end
            
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
        
        % Sort each array's data by requested stat, and put lines between
        % arrays.
        uq = unique(array_list); 
        sorts = cell(3,1);
        for j = 1:length(uq)
            
            % Find units in array
            unit_inds = find(strcmp(array_list, uq{j}));
            
            % Sort
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
            
            % Find spot for the graphical line between arrays
            array_end_ind = max(unit_inds);
            
            
            % Add data to show the lines
            if j ~= length(uq)  % ie skip top array in the fig
                nrows_insert_real = nrows_insert + (max_num_units_per_array(j) - length(unit_inds));
                num_rows_inserted = num_rows_inserted + nrows_insert_real;
                heatmap_mat = insertrows(heatmap_mat, repmat(zeros(1, length(kde_x_vals))+0.5, nrows_insert_real, 1), array_end_ind);
                array_list = insertrows(array_list', repelem({'NA'}, nrows_insert_real, 1), array_end_ind)';
            end
        end
        
        % Store array sorts
        sortid = sprintf('HeatmapSortInds_%s', sort_by);
        KDE(m).Sessions(sessn).(sortid) = sorts;
        
        % Add yticklabels for arrays
        ytick_spots = zeros(1, length(uq)+1);
        ytick_labs = cell(1, length(uq)+1);
        for j = 1:length(uq)
            unit_inds = find(strcmp(array_list, uq{j}));
            array_end_ind = max(unit_inds);
            ytick_spots(j) = round(median(unit_inds));
            ytick_labs{j} = sprintf('%s. (%d)', strcat(upper(uq{j}(1)), uq{j}(2:3)), numel(unit_inds));
        end
        
        % Set colormap:
        % 0,    0.5,    1
        % none  spacer  diff
        cm = parula(2);
        cm = [cm(1,:); divider_color; cm(2,:)];
        colormap(cm);
        
        % Plot data
        imagesc(heatmap_mat)
        
%         % Add ytick at top
%         max_yval = sum(max_num_units_per_array) + num_rows_inserted;
%         ytick_spots(end) = max_yval;
%         ytick_labs{end} = num2str(sum(max_num_units_per_array)); % this y-tick will be "wrong" to MATLAB but meaningful to humans
        
        % Sort to be ascending so MATLAB is happy
        [ytick_spots, idx] = sort(ytick_spots);
        ytick_labs = ytick_labs(idx);
        
        % Format plot
%         title(sprintf('%s', KDE(m).Sessions(sessn).ShortName), 'Color', col)
%         title(sprintf('%s', KDE(m).Sessions(sessn).ShortName))
        if i == 1
            xlabel('Time from cue on (ms)')
            ylabel('Array (# units)')
        end
        linkaxes(axes, 'y');
        xtick_idx = [51 151 251];
        xticks(xtick_idx)
        xticklabels(kde_x_vals(xtick_idx))
        ylim([0 sum(max_num_units_per_array) + 6])
        yticks(ytick_spots)
        yticklabels(ytick_labs)
        ytickangle(90)
        
        yl = ylim;
        plot(repelem(find(kde_x_vals == 0),2), [yl(1) yl(2)], 'w--')
        plot(repelem(find(kde_x_vals == 100),2), [yl(1) yl(2)], 'w--')
        
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',26, ...
            'fontname', 'Helvetica', 'fontweight', 'normal', ...
            'XColor', 'black', 'YColor', 'black', 'Color', 'none')
    end
%     set(gcf, 'Renderer', 'Painters');
%     set(gcf, 'Color', 'none');
%     saveas(gcf, fullfile(save_dir, sprintf('pv_KDEs_heatmaps_mk%s.pdf', KDE(m).Code)))
    
end

%% Plot avg timecourses (cat/dog diff) with days overlaid (learning)

% Params
areas_to_plot = {'anterior', 'middle', 'posterior'};
% areas_to_plot = {'te'};
% boundary = 'static'; % static or dynamic
% Not convinced correction is needed in this case.
% chi_sq_base_alpha = 0.05;
% chi_sq_alpha = chi_sq_base_alpha / numel(chi_sq_inds);
chi_sq_alpha = 0.05;
% colors = cbrewer('qual', 'Dark2', 3);
% colors = colors([3 2 1], :); % baseline, pre, post
colors = [0.5 0.5 0.5; mlc(1:2)];
range_normalize = false;
minimum_latency = 70;  % discount any latencies lower than this as noise

% Pre allocate
pre_data = cell(length(KDE), length(areas_to_plot));
post_data = cell(length(KDE), length(areas_to_plot));
chi_sq_inds = 1:25:251; % inds in kde_x_vals to test...corresponds to -200:50:500 in real time.


for m = 1:length(KDE)
    
    % Get sessions to use
    sessions_to_plot = rSessionsByMonk{m};
    
    % Set up figure
    figure2('Position', [200 200 1200 330]) % for arrays
    hold on
    
    for a = 1:length(areas_to_plot)
        area = areas_to_plot{a};
                        
        % Just arrays
        subplot(1, 3, a)
        hold on
        
        % Get data and plot
        for i = 1:length(sessions_to_plot)
            sessn = sessions_to_plot(i);
            
            % Choose correct units
            switch area
                case 'te'
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area));
                otherwise
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, area));
            end
            
            % Choose color to plot with
            if ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Base', 'once'))
                col_ind = 1;
            elseif ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre', 'once'))
                col_ind = 2;
            elseif ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post', 'once'))
                col_ind = 3;
            end
            
            % Get real x-axis values
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(1).CueOnAllCues.KDEXVals;
            
            % Collect data
            heatmap_mat = zeros(length(units_to_plot), length(kde_x_vals));
            for j = 1:length(units_to_plot)
                unit = units_to_plot(j);
                
                if isfield(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues, 'KDEYVals') && ...
                    isnan(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals)
                    continue
                end
            
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
            
            % Calculate mean
            heatmap_mean = mean(heatmap_mat); % 1 x num xvals
            if range_normalize
                heatmap_mean = normalize(heatmap_mean, 'range');
            end
            
            % Set xvals to use
%             xvals = 1:length(heatmap_mean);
            
            % Plot the data
            plot(kde_x_vals, heatmap_mean, ...
                'LineWidth', 1.5, ...
                'DisplayName', Monkeys(m).Sessions(sessn).ShortName,...
                'Color', colors(col_ind, :))
            
            % Format plot
%             title(sprintf('Monkey %s, %s', KDE(m).Code, area), 'Interpreter', 'none')
            title(sprintf('%s', area), 'Interpreter', 'none')
            xtick_idx = [51 151 251];
            xticks(kde_x_vals(xtick_idx))
            xticklabels(kde_x_vals(xtick_idx))
            if strcmp(area, 'te')
                xlabel('Time from cue on')
                ylabel('Fraction units')
            end
            
            % Store pre / post data for stats testing
            if ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post', 'once'))
                post_data{m,a} = heatmap_mat;
            elseif ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre', 'once'))
                pre_data{m,a} = heatmap_mat;
            end
        end
        
        % Do stats testing at each requested timepoint
        pre = pre_data{m,a};
        post = post_data{m,a};
        yl = get(gca, 'YLim');
        for ii = chi_sq_inds
            [h,p] = prop_test([sum(pre(:,ii)) sum(post(:,ii))], [size(pre,1) size(post,1)], false, chi_sq_alpha);
%             fprintf('%s, %s, %d \n', Monkeys(m).Sessions(sessn).ShortName, area, p)
            if h
                scatter(ii, yl(2)*0.97, 'ko', 'filled')
            end
        end
        
        % Janky test for population latency
        pre_diff_from_zero_bools = zeros(length(kde_x_vals), 1);
        post_diff_from_zero_bools = zeros(length(kde_x_vals), 1);
        minimum_latency_bool = (kde_x_vals < minimum_latency)';
        for ii = 1:length(kde_x_vals)
            [pre_diff_from_zero_bools(ii),~] = prop_test([sum(pre(:,ii)) 0], [size(pre,1) size(pre,1)], false, chi_sq_alpha);
            [post_diff_from_zero_bools(ii),~] = prop_test([sum(post(:,ii)) 0], [size(post,1) size(post,1)], false, chi_sq_alpha);
        end
        
        pre_popn_lat = kde_x_vals(find(pre_diff_from_zero_bools & ~minimum_latency_bool, 1));
        post_popn_lat = kde_x_vals(find(post_diff_from_zero_bools & ~minimum_latency_bool , 1));
        fprintf("Monkey %s, session %s, array %s \n\t POPN pre latency: %d  \n\t POPN post latency: %d \n", ...
            Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName, area, pre_popn_lat, post_popn_lat)
        
        if ~isempty(pre_popn_lat)
            plot(repelem(pre_popn_lat, 2), [0 yl(2)], '--', 'Color', mlc(1))
        end
        if ~isempty(post_popn_lat)
            plot(repelem(post_popn_lat, 2), [0 yl(2)], '--', 'Color', mlc(2))
        end
        
        % Format the plot some more
        if a == 1
            ylabel('Proportion signf. units')
            xlabel('Time from cue on')
        end
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',26, ...
            'fontname', 'Helvetica', ...
            'XColor', 'black', 'YColor', 'black')
        
    end
    % Format the plot some more
%     sgtitle(sprintf('Monkey %s', KDE(m).Code))
end

%% Latency of category coding, at a single-unit level
areas_to_plot = {'anterior', 'middle', 'posterior', 'te'};
minimum_latency = 70;  % discount any latencies lower than this as noise

for m = 1:length(KDE)

    sessions_to_plot = rSessionsByMonk{m};
    
    figure2('Position', [400 400 200*length(areas_to_plot) 200])
    hold on
    
    for a = 1:length(areas_to_plot)
        area = areas_to_plot{a};
                        
        % Just arrays
        subplot(1, length(areas_to_plot), a)
        hold on
    
        all_latencies = {};
        for i = 1:length(sessions_to_plot)
            sessn = sessions_to_plot(i);
            
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(1).CueOnAllCues.KDEXVals;

            % Choose correct units
            switch area
                case 'te'
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area));
                otherwise
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, area));
            end

            % Get data 
            latencies = nan(length(units_to_plot),1);
            for j = 1:length(units_to_plot)
                unit = units_to_plot(j);
                unit_signf_mat = zeros(length(kde_x_vals), 1);

                if isfield(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues, 'KDEYVals') && ...
                    isnan(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals)
                    continue
                end

                switch boundary
                    case 'static'
                        exceedid = 'ExceedBD_FixedBound';
                        unit_signf_mat(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid)) = 1;
                    case 'dynamic'
                        exceedid = 'ExceedBD_DynamicBound';
                        unit_signf_mat(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid)) = 1;
                    otherwise
                        error('Unexpected value for boundary')
                end
                
                early_idx = find(kde_x_vals < minimum_latency);
                unit_signf_mat(early_idx) = 0;
                latency = find(unit_signf_mat, 1);
                if ~isempty(latency)
                    latencies(j) = latency;
                else
                    continue
                end
            end
            
            % Plotting color and label
            switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
                case 'Base'
                    col = [0.5 0.5 0.5];
                case 'Pre'
                    col = mlc(1);
                    all_latencies{1} = latencies;
                case 'Post'
                    col = mlc(2);
                    all_latencies{2} = latencies;
            end
            
            
            % Histogram of latencies for this array / session
            histogram(latencies, 'BinEdges', 70:20:400, 'FaceColor', mlc(i), 'Normalization', 'probability')
            
            % Draw mean / std on plots
            mu = nanmean(latencies);
            sigma = nanstd(latencies);
            yl = ylim;
            marker_yval = yl(2) * 1.01 + 0.002 * (i-1);
            plot([mu-sigma mu+sigma], repelem(marker_yval, 2), '-', 'Color', mlc(i), 'LineWidth', 1.5)
            scatter(mu, marker_yval , 100, mlc(i), 'filled', 'v',  'MarkerEdgeColor', 'k')
            xlabel("Latency")
            if a  == 1
                ylabel("Probability")
            end
            
        end
        
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',26, ...
            'fontname', 'Helvetica', ...
            'XColor', 'black', 'YColor', 'black')
        
        [p, h] = ranksum(all_latencies{1}, all_latencies{2});
        disp([nanmedian(all_latencies{1}) nanmedian(all_latencies{2})])
        fprintf("Monkey %s, session %s, array %s, pre vs post latencies, pval %0.4g \n", ...
            Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName, area, p)
        
    end
end

%% Duration of signf cat/dog diffs (pre/post overlay)

% Parameters
ranksum_alpha = 0.05;
% areas_to_plot = {'anterior', 'middle', 'posterior', 'te'};
areas_to_plot = {'te'};
figs_path = '/Users/jonahpearl/Documents/BJR group/Catdog paper';

for a = 1:length(areas_to_plot)
    area = areas_to_plot{a};
    
    % Plot both mks together
%     figure2('Position', [200 200 1000 1000])
    figure2('Position', [200 200 300 300])
    hold on
    
    for m = 1:length(KDE)
        sessions_to_plot = rSessionsByMonk{m};
        
        % Arrange subplots
        subplot(2,1,m)
        hold on

        for i = 1:length(sessions_to_plot)
            sessn = sessions_to_plot(i);
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(1).CueOnAllCues.KDEXVals;
            
            switch area
                case 'te'
                    units_to_plot = find(strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, 'te'));
                otherwise
                    units_to_plot = find(strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Location}, area));
            end
            runs = cell(length(units_to_plot), 2); % runs of significant KDE difference. First col is start time in ms of run, second col is length of run.

            for j = 1:length(units_to_plot)
                unit = units_to_plot(j);
                
                if isfield(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues, 'KDEYVals') && ...
                    isnan(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals)
                    continue
                end

                % Get data
                switch boundary
                    case 'static'
                        exceedid = 'ExceedBD_FixedBound';
                    case 'dynamic'
                        exceedid = 'ExceedBD_DynamicBound';
                    otherwise
                        error('Unexpected value for boundary')
                end
                exceed_bd_inds = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid);
                actual_xvals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals;
                max_val = length(actual_xvals);

                % Catch units with no signf time points
                if isempty(exceed_bd_inds)
                    runs{j,1} = 0;
                    runs{j,2} = 0;
                    continue
                end

                % Convert to string of zeros and ones for textscan
                exceed_bd_bools = zeros(1, max_val);
                exceed_bd_bools(exceed_bd_inds) = 1;
                s = sprintf('%d', exceed_bd_bools);

                % Find runs with textscan, then get their lengths
                runs{j,1} = (find(diff(uint8(s)))+1)'; % inds of first 1 in each run of 1's
                s = textscan(s, '%s', 'delimiter', '0', 'multipleDelimsAsOne',1); % all runs of 1s
                s = s{:}; % unpack
                runs{j,2} = cellfun('length', s); % length of each run
            end

            % Concatenate all the run lengths across all units into one
            % vector. Works even when one unit has multiple runs.
            data = vertcat(runs{:,2});

            % Remove zeros -- we're only interested in units with a significant
            % cat/dog difference
            data(data == 0) = [];
            
            % Account for step size in kde x vals
            step = mean(diff(kde_x_vals));  % should only ever be 1 or 2, atm
            data = data * step;

            % Plot the data as a histogram        
            histogram(data, 'BinEdges', [0:25:500 Inf], ... % 500
                'Normalization', 'probability', ...
                'DisplayName', Monkeys(m).Sessions(sessn).ShortName, ...
                'FaceColor', mlc(i))

            % Add mean/std above the histogram
%             yl = ylim();
            mu = mean(data);
            sigma = std(data);
            y_marker_val = 0.7 + 0.05 * (i-1);
            scatter(mu, y_marker_val, 20, '^', 'MarkerFaceColor', mlc(i), 'MarkerEdgeColor', mlc(i))
            plot([mu - sigma, mu + sigma], repelem(y_marker_val, 2), ...
                'LineWidth', 1, 'Color', mlc(i), 'MarkerFaceColor', mlc(i))
            
            % Store the data
            run_id = sprintf('KDE_EBD_Runs_%s', area); 
            KDE(m).Sessions(sessn).(run_id)= runs;
            if strcmp('Post', regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once'))
                post_data = data;
            end

            % Format the plot
%             title(sprintf('Lengths of signf. runs, %s, session %s, %s',...
%                 KDE(m).Name, KDE(m).Sessions(sessn).ShortName, area), 'Interpreter', 'none')
            title(sprintf('%s, %s', Monkeys(m).Name, area), 'interpreter', 'none')
%             ylim([0 0.3])
            yticks(0:0.2:0.4)
            xlim([0 450])
            
    %         legend
            if i == 1 && m == 1
%                 xlabel('Run lengths of cat/dog significanct differences (ms)')
                xlabel('Run length (ms)')
    %             ylabel({'Probability across all runs','(Not norm''d to num unts)'})
                ylabel('Probability')
            end
            set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
                'XMinorTick', 'off', 'YMinorTick', 'off',...
                'fontsize',14, ...
                'fontname', 'Helvetica', 'fontweight', 'normal', ...
                'XColor', 'black', 'YColor', 'black')
        end
        
        % Run stats testing
        for i = 1:length(sessions_to_plot)
            sessn = sessions_to_plot(i);
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(1).CueOnAllCues.KDEXVals;
            
            if strcmp('Post', regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once'))
                continue
            end
            step = mean(diff(kde_x_vals));  % should only ever be 1 or 2, atm
            run_id = sprintf('KDE_EBD_Runs_%s', area); 
            runs = KDE(m).Sessions(sessn).(run_id);
            data = vertcat(runs{:,2});
            data = data * step;
            if isempty(data)
                fprintf('%s, session %d, area %s: no signf timepoints \n', KDE(m).Name, sessn, area, h, p)
                continue
            end
            data(data == 0) = [];  % remove zeros
            [p,h] = ranksum(data, post_data, 'Alpha', ranksum_alpha);
            fprintf('%s, session %d, area %s, 1-sided ranksum session vs post: h = %d, p = %0.3f \n', KDE(m).Name, sessn, area, h, p)
        end
    end
    
    % Save figure
    save_path = fullfile(figs_path, sprintf('KDE_diffs_duration_%s.png', area));
%     saveas(gcf, save_path)
end

% Format the plot more
% sgtitle('Session, median +/- iqr')
set(gcf, 'Color', 'white')

%% testing assym kernels
% 
% x = -100:400;
% y = asymmetric_gaussian_kernel(x);
% figure
% plot(x, y)
% 
% spike_times = [-50 10 80 90 100 110 112 113 150];
% spike_times = [spike_times (spike_times + 10)];
% y1 = ksdensity(spike_times, x, 'Kernel', @asymmetric_gaussian_kernel, "Bandwidth", 15);
% y2 = ksdensity(spike_times, x, "Kernel", @gaussian_kernel, "Bandwidth", 15);
% disp(trapz(y1))
% disp(trapz(y2))
% 
% figure
% hold on
% plot(x, y1)
% plot(x, y2)
% legend(["custom", "bw=15"])


% x = -100:100;
% figure
% plot(x, asymmetric_response_kernel(x))

% x = 0:100;
% figure
% plot(x, synaptic_kernel(x))


%% Functions

function y = gaussian_kernel(x)
    sigma = 1; % Standard deviation of the Gaussian kernel
    y = exp(-0.5 * (x / sigma).^2) / (sigma * sqrt(2 * pi));
end


function y = asymmetric_gaussian_kernel(x)
    
    % Connor lab uses 15/5, so define a ratio of 1:3 here and then pass bw
    % = 15 to the ksdensity function.
    sigma1 = 1;
    sigma2 = 1/3;
    
    % Allocate output array
    y = zeros(size(x));
    
    % Indices for positive and negative x
    positive = x >= 0;
    negative = x < 0;
    
    % Calculate kernel values
    y(positive) = exp(-0.5 * (x(positive) / sigma1).^2) / (sigma1 * sqrt(2 * pi));
    y(negative) = exp(-0.5 * (x(negative) / sigma2).^2) / (sigma1 * sqrt(2 * pi));
    y = y * 1.5;   % use a scaling factor to get area underneath to still be 1 
end

% 
% function y = asymmetric_response_kernel(t)
%     tau_g = 5;  % Growth time constant in ms
%     tau_d = 15; % Decay time constant in ms
%     
%     y = zeros(size(t));
%     growth_side = t < 0;
%     decay_side = t >= 0;
%     
%     % Applying growth formula to negative t (t approaches 0 from the left)
%     y(growth_side) = exp(t(growth_side) / tau_g);
%     
%     % Applying decay formula to positive t (t moves away from 0)
%     y(decay_side) = exp(-t(decay_side) / tau_d);
% end
% 
% function y = synaptic_kernel(t)
%     tau_g = 5;  % Growth time constant in ms
%     tau_d = 15; % Decay time constant in ms
% 
%     % Ensure t is non-negative for this formula to make sense
%     y = zeros(size(t));
%     positive = t >= 0;
%     y(positive) = (1 - exp(-t(positive) / tau_g)) .* exp(-t(positive) / tau_d);
% end

