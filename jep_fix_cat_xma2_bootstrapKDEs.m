% looks for significant differences in cat/dog responses by bootstrapping

%% Updates
% updated may 6 to use real trial count (previously only counted trials
% with spikes, which was most but not necessarily all trials).

%% Load data
bw = 20;
path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = sprintf('MaxMarta_xma2_dualKDEs_bw%d.mat', bw);
Data = load(fullfile(path1, fname));
[status, KDE] = stitch_monkeyStruct_from_parts(Data);
clear Data
clearvars -except KDE bw

%% Add CueInfo, get number of cat/dog trials, then remove CueInfo
path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = 'MaxMarta_xma2_raw_behav.mat';
catgs = {1:260, 261:520};

% load data
Monkeys = load(fullfile(path1, fname));
Monkeys = Monkeys.Monkeys;

for m = 1:length(KDE)
    for i = 1:length(KDE(m).Sessions)
        
        % error check
        if ~strcmp(Monkeys(m).Sessions(i).SessionName, KDE(m).Sessions(i).SessionName)
            error('Session names do not match')
        else % move into KDE struct
            num_total_cat_trials = sum([Monkeys(m).Sessions(i).CueInfo(catgs{1}).NumApp]);
            num_total_dog_trials = sum([Monkeys(m).Sessions(i).CueInfo(catgs{2}).NumApp]);
            KDE(m).Sessions(i).NumCatgTrials = [num_total_cat_trials num_total_dog_trials];
        end
    end
end
clearvars -except KDE bw

%% Store real cat/dog diffs

diffid = sprintf('TrueDiff_BW%d', bw);

for m = 1:length(KDE)
    % for m = 1
    sessions_to_plot = 1:length(KDE(m).Sessions);
    %     sessions_to_plot = [7 9];
    KDE(m).Name = KDE(m).Name;
    
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        units_to_plot = 1:length(KDE(m).Sessions(sessn).UnitInfo);
        
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            
            if isempty(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterXVals)
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid) = NaN;
                continue
            end
            
            % get data
            raster_x_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterXVals; % normalized spike times
            raster_y_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterYVals; % useful for keeping trials together in the shuffle
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals; % normalized time
            catg_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.CatgVals; % cat or dog for each spike
            num_total_cat_trials = KDE(m).Sessions(sessn).NumCatgTrials(1);
            num_total_dog_trials = KDE(m).Sessions(sessn).NumCatgTrials(2);
            
            % Get real cat-dog diff.
            % scale like so: (density / ms) x (1000 ms / s) x (num spikes) /
            % (num trials)  = spikes / second / trial.
            % See Demo_KDE_scaling if you don't believe me :D
            catg1_scaled = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg1_raw * 1000 * sum(catg_vals==1) / num_total_cat_trials;
            catg2_scaled = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg2_raw * 1000 * sum(catg_vals==2) / num_total_dog_trials;
            true_diff = catg1_scaled - catg2_scaled;
            
            % store
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid) = true_diff;
        end
    end
end

%% (REALLY LONG -- like 12 hours) Bootstrap

n_bootstrap = 50;
rng(10)
catgs = {1:260, 261:520};
% debug_ratios = zeros(1, n_bootstrap);
bootid = sprintf('BootstrappedDiffs_BW%d', bw);

for m = 1:length(KDE)
    % for m = 2
    %     sessions_to_plot = 1:length(KDE(m).Sessions);
    %     sessions_to_plot = 1;
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_plot = [1 2 7 9]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_plot = [1:3 5:9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_plot = [1 2 6 7]; % This is base 01, base02, pre and post. skip base03,04,05 and sub03 (low trial count)
        sessions_to_plot = 1:7;
    end
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        units_to_plot = find(ismember({KDE(m).Sessions(sessn).UnitInfo.Location}, {'anterior', 'middle', 'posterior'}));
        %         units_to_plot = 1:length(KDE(m).Sessions(sessn).UnitInfo);
        %         units_to_plot = 102;
        
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            
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
            bootstrapped_diffs = zeros(n_bootstrap, length(kde_x_vals));
            shuffled_catgs = catg_vals; % to be shuffled
            n_uq = length(uq);
            
            % bootstrap
            for n = 1:n_bootstrap
                
                
                % OLD WAY
                % Shuffle categories while keeping trials intact, but
                % cats/dogs UNBALANCED (ie, if cats have 1000 trials with spikes and
                % dogs have 500, it will stay like that).
                %                 img_shuffle = rounded_uq(randperm(n_uq)); % 1 x num trials, with shuffled image numbers
                %                 for t = 1:length(uq)
                %                     if img_shuffle(t) <= 260
                %                         shuffled_catgs(trial_inds{t}) = 1;
                %                     elseif img_shuffle(t) > 261
                %                         shuffled_catgs(trial_inds{t}) = 2;
                %                     end
                %                 end
                
                % NEW WAY
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
                
                % debugging the shuffle
                %                 ratio = sum(shuffled_catgs==1) / sum(shuffled_catgs==2);
                %                 fprintf('ratio %0.2f \n', ratio)
                %                 debug_ratios(n) = ratio;
                
                
                % fit KDEs
                shuffled_f_catg1 = ksdensity(raster_x_vals(shuffled_catgs==1), kde_x_vals, 'Bandwidth', bw);
                shuffled_f_catg2 = ksdensity(raster_x_vals(shuffled_catgs==2), kde_x_vals, 'Bandwidth', bw);
                
                
                % OLD WAY
                % scale KDEs
                %                 shuffled_f_catg1 = shuffled_f_catg1 * 1000 * sum(shuffled_catgs==1) / sum(rounded_img_num <= 260);
                %                 shuffled_f_catg2 = shuffled_f_catg2 * 1000 * sum(shuffled_catgs==2) / sum(rounded_img_num >= 261);
                
                % NEW WAY
                % scale KDEs
                shuffled_f_catg1 = shuffled_f_catg1 * 1000 * sum(shuffled_catgs==1) / num_total_cat_trials;
                shuffled_f_catg2 = shuffled_f_catg2 * 1000 * sum(shuffled_catgs==2) / num_total_dog_trials;
                
                % get diff
                bootstrapped_diffs(n,:) = shuffled_f_catg1 - shuffled_f_catg2;
            end
            % store
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid) = bootstrapped_diffs;
            fprintf('Done with %s, session %d, unit %d \n', KDE(m).Name, sessn, unit)
        end
    end
end

%% Save data
path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = sprintf('MaxMarta_fix_cat_xma2_bootstrappedKDEs_BW%d.mat', bw);
[status, vnames] = split_monkeyStruct_in_parts(KDE);
save(fullfile(path1, fname), vnames{:});

%% Re load data, add Area labels
clearvars
close all
bw = 20;
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

%% Get short names
marta_xls = './RecordingMarta.xlsx';
max_xls = './RecordingMax.xlsx';

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

%% (moving/timecourse boundary) Get xvals where true diffs are larger than bootstrapped popn

pct = 99;
diffid = sprintf('TrueDiff_BW%d', bw);
bootid = sprintf('BootstrappedDiffs_BW%d', bw);

for m = 1:length(KDE)
    % for m = 1
    sessions_to_plot = 1:length(KDE(m).Sessions);
    %     sessions_to_plot = [7 9];
    KDE(m).Name = KDE(m).Name;
    
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        units_to_plot = 1:length(KDE(m).Sessions(sessn).UnitInfo);
        
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals;
            true_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid);
            boot_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid);
            boot_diffs = abs(boot_diffs); % just look at magnitude, not sign
            percentile_diffs = prctile(boot_diffs, pct, 1); % 1 x (num xvals) vector of the 95th percentile of bootstrapped data
            exceedid = 'ExceedBD_DynamicBound';
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid) = find(abs(true_diffs) > percentile_diffs); % normalized times in ms where true diffs exceeds bootstrapped diffs
            
        end
    end
end

%% (fixed boundary) Get xvals where true diffs are larger than bootstrapped popn

pct = 99;
diffid = sprintf('TrueDiff_BW%d', bw);
bootid = sprintf('BootstrappedDiffs_BW%d', bw);

for m = 1:length(KDE)
    % for m = 1
   if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_plot = [1 2 7 9]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_plot = [1:3 5:9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_plot = [1 2 6 7]; % This is base 01, base02, pre and post. skip base03,04,05 and sub03 (low trial count)
        sessions_to_plot = 1:7;
   end
    
    KDE(m).Name = KDE(m).Name;
    
    
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
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid) = find(abs(true_diffs) > percentile_diff); % normalized times in ms where true diffs exceeds bootstrapped diffs
            
        end
    end
end

%% ***** Boolean timecourse plots *****

%% Plot dynamic boundary example

m = 2;
sessn = 1;
unit = 40;
bw = 20;
pct = 99;

diffid = sprintf('TrueDiff_BW%d', bw);
bootid = sprintf('BootstrappedDiffs_BW%d', bw);
true_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid);
boot_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid);

percentile_diffs = prctile(abs(boot_diffs), pct, 1); % 1 x (num xvals) vector of the 95th percentile of bootstrapped data
t = find(abs(true_diffs) > percentile_diffs);

figure
hold on
plot(abs(true_diffs), 'g-', 'LineWidth', 2)
plot(abs(boot_diffs)', '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1)
plot(percentile_diffs, 'r--', 'LineWidth', 1.5)
yl = get(gca, 'YLim');
scatter(t, repelem(yl(2)*0.98, numel(t)), 'go', 'filled')
xlabel('Time from cue on (ms)')
xticks(0:200:700)
xticklabels(-200:200:500)
ylabel('Abs. cat/dog difference (sp/sec/tr)')
title(sprintf('%s, session %s, unit %d', KDE(m).Name, KDE(m).Sessions(sessn).ShortName, unit), 'Interpreter', 'none')
make_custom_patch_legend({'g', [0.4 0.4 0.4], 'r'}, {'Abs. Diff', 'Abs. Boot', 'Abs. 99th percentile boot'})
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',18, ...
    'fontname', 'Helvetica', 'fontweight', 'normal', ...
    'XColor', 'black', 'YColor', 'black')

%% Plot static boundary example

m = 2;
sessn = 1;
unit = 40;
bw = 20;
pct = 99;

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
yticks([0 9])
ylim([0 9])
% title(sprintf('%s, session %s, unit %d', KDE(m).Name, KDE(m).Sessions(sessn).ShortName, unit), 'Interpreter', 'none')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',18, ...
    'fontname', 'Helvetica', 'fontweight', 'normal', ...
    'XColor', 'black', 'YColor', 'black')
make_custom_patch_legend({'r', 'b'}, {'Cats', 'Dogs'}, 'FontSize', 10, 'Location', 'northwest')

subplot(2,1,2)
hold on
plot(kde_x_vals, abs(true_diffs), 'g-', 'LineWidth', 2)
plot(kde_x_vals, abs(boot_diffs)', '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1)
plot([-200 500], repelem(percentile_diff, 2,1), 'r--', 'LineWidth', 1.5)
yl = get(gca, 'YLim');
scatter(kde_x_vals(t), repelem(yl(2)*0.98, numel(t)), 'go', 'filled')
xlabel('Time from cue on (ms)')
xlim([-200 500])
xticks(kde_x_vals(1):200:kde_x_vals(end))
ylabel('Spikes / second / trial')
yticks([0 2])
% title(sprintf('%s, session %s, unit %d', KDE(m).Name, KDE(m).Sessions(sessn).ShortName, unit), 'Interpreter', 'none')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',18, ...
    'fontname', 'Helvetica', 'fontweight', 'normal', ...
    'XColor', 'black', 'YColor', 'black')
make_custom_patch_legend({'g', [0.4 0.4 0.4], 'r'}, {'Abs. Diff', 'Abs. Shuff. Diff', '99th percentile (static)'}, 'FontSize', 10, 'Location', 'northwest')
set(gcf, 'Color', 'white')

%% Plot results as heatmap

area_to_plot = 'te';
boundary = 'static'; % static or dynamic
nrows_insert = 3; % size of red line separating arrays
sort_by = 'after_zero_latency'; % after_zero_latency or duration
sort_empty_val = 1000; % controls if units with no signf diff go above or below others

for m = 1:length(KDE)
    
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_plot = [1 2 7 9]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_plot = [1:3 5:9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_plot = [1 2 6 7]; % This is base 01, base02, pre and post. skip base03,04,05 and sub03 (low trial count)
        sessions_to_plot = 1:7;
    end
    
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
            ytick_labs{j} = sprintf('%s (%d)', uq{j}, numel(unit_inds));
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
        title(sprintf('%s', KDE(m).Sessions(sessn).ShortName))
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
            'fontsize',18, ...
            'fontname', 'Helvetica', 'fontweight', 'normal', ...
            'XColor', 'black', 'YColor', 'black')
    end
    set(gcf, 'Color', 'white')
end

%% Plot avg timecourses (cat/dog diff) with days overlaid (learning)

% params
% areas_to_plot = {'te', 'anterior', 'middle', 'posterior'};
areas_to_plot = {'te'};
boundary = 'static'; % static or dynamic
% Not convinced correction is needed in this case.
% chi_sq_base_alpha = 0.05;
% chi_sq_alpha = chi_sq_base_alpha / numel(chi_sq_inds);
chi_sq_alpha = 0.05;

% pre allocate
pre_data = cell(length(KDE), length(areas_to_plot));
post_data = cell(length(KDE), length(areas_to_plot));
chi_sq_inds = 200:50:700; % inds in kde_x_vals to test...corresponds to -200:50:500 in real time.



for m = 1:length(KDE)
    
    % get sessions to use
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
        %                 sessions_to_plot = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_plot = [1 2 7 9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
        %                 sessions_to_plot = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
        sessions_to_plot = [1 2 6 7];
    end
    
    
    % prep figure
    %     figure('Position', [100 100 1100 900]) % for te + arrays
    figure('Position', [400 400 500 300]) % just te
    hold on
    
    for a = 1:length(areas_to_plot)
        area = areas_to_plot{a};
        
        %         % configure subplots
        %         switch area
        %             case 'te'
        %                 subplot(2,3,1:3)
        %             otherwise
        %                 subplot(2,3,2+a)
        %         end
        % %         subplot(2, 2, a)
        %         hold on
        
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
            %             plot(movmean(heatmap_mean,20), 'DisplayName', KDE(m).Sessions(sessn).ShortName)
            plot(heatmap_mean, 'DisplayName', KDE(m).Sessions(sessn).ShortName, 'LineWidth', 1.5)
            
            % format plot
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
        
        % format plot
        %         legend('Location', 'eastoutside')
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',18, ...
            'fontname', 'Helvetica', ...
            'XColor', 'black', 'YColor', 'black')
    end
    %     sgtitle(sprintf('%s (Chi-sq, p < %0.2f with Bonferroni)', KDE(m).Name, chi_sq_base_alpha),  'Interpreter', 'none')
    %     sgtitle(sprintf('%s (Chi-sq, p < %0.2f)', KDE(m).Name, chi_sq_alpha),  'Interpreter', 'none')
end

%% Plot avg timecourses (cat/dog diff) with arrays overlaid (kinetics)

areas_to_plot = {'anterior', 'middle', 'posterior'};
boundary = 'static'; % static or dynamic
% colors = cbrewer('qual', 'Set2', 3);
colors = cbrewer('qual', 'Paired', 6);
light_colors = colors([1 3 5], :)*0.9;
dark_colors = colors([2 4 6], :);

% overlay post timecourses of different arrays to show kinetics
figure('Position', [200 200 1000 400])
hold on
for m = 1:length(KDE)
    
    % get sessions to use
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
        %                 sessions_to_plot = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_plot = [1 2 7 9];
        %         sessions_to_plot = [9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
        %                 sessions_to_plot = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
        sessions_to_plot = [1 2 6 7];
        %         sessions_to_plot = [7];
    end
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        
        % subplots of monks x sessions
        subplot(length(KDE), length(sessions_to_plot), length(sessions_to_plot)*(m-1) + mod(i-1, length(sessions_to_plot)) + 1)
        hold on
        
        for a = 1:length(areas_to_plot)
            area = areas_to_plot{a};
            
            % get correct units
            switch area
                case 'te'
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area));
                otherwise
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, area));
            end
            
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
            
            %             plot(normalize(mean(heatmap_mat), 'range'), 'Color', colors(a,:), ...
            %                 'DisplayName', areas_to_plot{a}, 'LineWidth', 1.5)
            plot(mean(heatmap_mat), 'Color', dark_colors(a,:), ...
                'DisplayName', areas_to_plot{a}, 'LineWidth', 1.5)
            
        end
        %         legend
        title(sprintf('%s', KDE(m).Sessions(sessn).ShortName), 'Interpreter', 'none')
        xticks(0:200:600)
        xticklabels(kde_x_vals(1:200:end))
        if m == 1 && i == 3
            xlabel('Time from cue on')
            ylabel('Fraction units')
        end
        ylim([0 0.7])
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',18, ...
            'fontname', 'Helvetica',  ...
            'XColor', 'black', 'YColor', 'black')
    end
end

%% ***** Magnitude of diff timecourse plots *****

%% Plot magnitude of diff as similar heatmap

area_to_plot = 'te';
nrows_insert = 3; % size of red line separating arrays
sort_by = 'after_zero_latency'; % after_zero_latency or duration
sort_empty_val = 1000; % controls if units with no signf diff go above or below others
diffid = sprintf('TrueDiff_BW%d', bw);

for m = 1:length(KDE)
    
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
        %         sessions_to_use = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_plot = [1 2 7 9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
        %         sessions_to_use = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
        sessions_to_plot = [1 2 6 7];
    end
    
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
        
        % take only te units
        area_bool = strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area_to_plot);
        units_to_plot = find(area_bool);
        
        % sort units by array
        array_list = {KDE(m).Sessions(sessn).UnitInfo(area_bool).Location};
        [array_list, idx] = sort(array_list);
        units_to_plot = units_to_plot(idx);
        
        % get data (already sorted by array)
        % normalize here within each unit, if desired
        heatmap_mat = zeros(length(units_to_plot), length(kde_x_vals));
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            data = abs(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid));
            exceedid = 'ExceedBD_FixedBound';
            
            % if not enough spikes, or bootstrap analysis is never
            % significant, leave the row as zeros.
            if isnan(data)
                continue
            elseif isempty(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid))
                continue
            end
            %             mean_baseline = mean(data(100:200));
            %             std_baseline = std(data(100:200));
            %             heatmap_mat(j, :) = (data - mean_baseline) / std_baseline;
            heatmap_mat(j, :) = normalize(data, 'range');
            %             heatmap_mat(j, :) = data;
        end
        
        percentile_99 = prctile(reshape(heatmap_mat, 1, numel(heatmap_mat)), 99);
        
        % sort data in each array, put lines between arrays
        uq = unique(array_list);
        for j = 1:length(uq)
            
            % find units in array
            unit_inds = find(strcmp(array_list, uq{j}));
            full_data = heatmap_mat(unit_inds, :);
            
            % sort like before
            sortid = sprintf('HeatmapSortInds_%s', sort_by);
            idx = KDE(m).Sessions(sessn).(sortid){j};
            heatmap_mat(unit_inds,:) = full_data(idx,:);
            
            % add data for the lines
            array_end_ind = max(unit_inds);
            heatmap_mat = insertrows(heatmap_mat, repmat(zeros(1, length(kde_x_vals))-1, nrows_insert,1), array_end_ind);
            array_list = insertrows(array_list', repelem({'NA'}, nrows_insert, 1), array_end_ind)';
        end
        
        % add yticklabels for arrays
        ytick_spots = zeros(1, length(uq)+1);
        ytick_labs = cell(1, length(uq)+1);
        for j = 1:length(uq)
            unit_inds = find(strcmp(array_list, uq{j}));
            array_end_ind = max(unit_inds);
            ytick_spots(j) = round(median(unit_inds));
            ytick_labs{j} = sprintf('%s (%d)', uq{j}, numel(unit_inds));
        end
        
        % set colormap
        cm = parula(100);
        cm = [[0.75 0.75 0.75]; cm];
        colormap(cm);
        caxis([-0.1 1]);
        %         caxis([-0.1 3]);
        %         caxis([-0.1 percentile_99])
        colorbar
        
        % plot data
        imagesc(heatmap_mat)
        
        % add ytick at top
        ytick_spots(end) = max_num_units;
        ytick_labs{end} = num2str(max_num_units);
        
        % sort to be ascending so MATLAB is happy
        [ytick_spots, idx] = sort(ytick_spots);
        ytick_labs = ytick_labs(idx);
        
        % format plot
        title(sprintf('%s', KDE(m).Sessions(sessn).ShortName))
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
            'fontsize',18, ...
            'fontname', 'Helvetica', 'fontweight', 'normal', ...
            'XColor', 'black', 'YColor', 'black')
    end
    set(gcf, 'Color', 'white')
end

%% Get cat/dog diffs' latencies and max-times/max-vals thereof

bw = 20;
boundary = 'static'; % static or dynamic
diffid = sprintf('TrueDiff_BW%d', bw);
num_consec_needed = 20;

for m = 1:length(KDE)
    %     sessions_to_plot = 1:length(KDE(m).Sessions);
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
        %         sessions_to_plot = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_plot = [1 2 7 9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
        %         sessions_to_plot = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
        sessions_to_plot = [1 2 6 7];
    end
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        %         units_to_plot = 1:length(KDE(m).Sessions(sessn).UnitInfo);
        units_to_plot = find(ismember({KDE(m).Sessions(sessn).UnitInfo.Location}, {'anterior', 'middle', 'posterior'}));
        
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            
            % choose type of analysis
            switch boundary
                case 'static'
                    exceedid = 'ExceedBD_FixedBound';
                    latid = 'ExceedBD_FixedBound_Latency';
                case 'dynamic'
                    exceedid = 'ExceedBD_DynamicBound';
                    latid = 'ExceedBD_DynamicBound_Latency';
                otherwise
                    error('Unexpected value for boundary')
            end
            
            % get data
            data = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid);
            
            % latency to first diff exceeding bootstrapped diffs
            if ~issorted(data)
                error('Exceed BD inds not sorted')
            elseif isempty(data)
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(latid) = NaN;
            else
                diffs = diff(data);
                inds = strfind(diffs, repelem(1, num_consec_needed));
                if isempty(inds)
                    KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(latid) = NaN;
                else
                    KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(latid) = data(inds(1));
                end
            end
            
            
            % latency to largest diff (doesn't depend on bootstraps)
            true_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid);
            [max_val , max_t] = max(true_diffs);
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.Diff_MaxVal = max_val;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.Diff_MaxTime = max_t;
        end
    end
end

%% Barplot of max diffs by array
arrays = {'anterior', 'middle', 'posterior'};
matlab_colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3;
    0          0        0    ];
exceedid = 'ExceedBD_FixedBound';
 
figure

for m = 1:length(KDE)
    %     sessions_to_plot = 1:length(KDE(m).Sessions);
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
        %         sessions_to_plot = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_plot = [1 2 7 9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
        %         sessions_to_plot = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
        sessions_to_plot = [1 2 6 7];
    end
    
    barmat_means = zeros(length(sessions_to_plot), length(arrays));
    barmat_stds = zeros(length(sessions_to_plot), length(arrays));
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        
        for a = 1:length(arrays)
            array = arrays{a};
            units = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, array));
            concat = zeros(1,numel(units));
            for j = 1:length(units)
                if isempty(KDE(m).Sessions(sessn).UnitInfo(units(j)).CueOnAllCues.(exceedid))
                    concat(j) = nan;
                else
                    concat(j) = KDE(m).Sessions(sessn).UnitInfo(units(j)).CueOnAllCues.Diff_MaxVal;
                end
                barmat_means(i,a) = mean(concat, 'omitnan');
                barmat_stds(i,a) = std(concat, 'omitnan');
            end
        end
    end
    
    subplot(2,1,m)
    hold on
    b = bar(barmat_means);
    pause(0.5)
    xoff = vertcat(b.XOffset);
    xoff = repmat(xoff, 1, 4);
    xvals = repmat(1:length(sessions_to_plot), numel(arrays), 1);
    errorbar(xvals'+xoff', barmat_means, barmat_stds, 'ko')
    if m == 1
        ylim([-1.75 4])
    elseif m == 2
        ylim([-1 3])
        xlabel('Session')
        ylabel('Cat/dog difference (spikes / sec / trial)')
        make_custom_patch_legend(matlab_colors(1:3,:), arrays, 'FontSize', 12)
    end
    xticks(1:4)
    xticklabels({KDE(m).Sessions(sessions_to_plot).ShortName})
    xtickangle(45)
    yticks([0 3])
    plot([0 4], [0 0], 'k-', 'HandleVisibility', 'off')
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',18, ...
        'fontname', 'Helvetica', 'fontweight', 'normal', ...
        'XColor', 'black', 'YColor', 'black')
    
end
set(gcf, 'Color', 'white')

%% ***** Add visual responsiveness data *****

%% Load EXCIT / SUPPRN latency data, add Area labels to it

path1 = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
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

%% (note) Window results key and window results summary key:

%{
Raw results (each time point):
0: suppressed
1: nothing
2: excited
%}

%{
Summary (three 0's in a row or three 2's in a row or both):
0: suppressed
1: nothing
2: exicted
3: suppressed then excited
4: excited then suppressed
%}

%% Plot avg timecourse (EXCIT) with arrays overlaid (kinetics)

areas_to_plot = {'anterior', 'middle', 'posterior'};
boundary = 'static'; % static or dynamic
colors = cbrewer('qual', 'Set1', 3);
excit_supprn_timecourse_xvals = -50:10:430;

% overlay post timecourses of different arrays to show kinetics
figure('Position', [200 200 1200 400])
hold on
for m = 1:length(KDE)
    
    % get sessions to use
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
        sessions_to_plot = [1 2 3 5 6 7 9]; % skip base04 (outlier) and post01 (low trial count)
%         sessions_to_plot = [1 2 7 9];
%         sessions_to_plot = [9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
        sessions_to_plot = [1 2 3 4 5 6 7]; % skip base04 and sub03 (low trial count)
%         sessions_to_plot = [1 2 6 7];
        %         sessions_to_plot = [7];
    end
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        
        % subplots of monks x sessions
        subplot(length(KDE), length(sessions_to_plot), length(sessions_to_plot)*(m-1) + mod(i-1, length(sessions_to_plot)) + 1)
        hold on
        
        for a = 1:length(areas_to_plot)
            area = areas_to_plot{a};
            
            % get correct units
            switch area
                case 'te'
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area));
                otherwise
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, area));
            end
            
            % get data
            heatmap_mat = zeros(length(units_to_plot), length(excit_supprn_timecourse_xvals));
            for j = 1:length(units_to_plot)
                unit = units_to_plot(j);
                heatmap_mat(j, [Monkeys(m).Sessions(sessn).UnitInfo(unit).WindowResults]==2) = 1;
            end
            
            %             plot(normalize(mean(heatmap_mat), 'range'), 'Color', colors(a,:), ...
            %                 'DisplayName', areas_to_plot{a}, 'LineWidth', 1.5)
            plot(mean(heatmap_mat), 'Color', colors(a,:), ...
                'DisplayName', areas_to_plot{a}, 'LineWidth', 1.5)
            
        end
        %         legend
        title(sprintf('%s', KDE(m).Sessions(sessn).ShortName), 'Interpreter', 'none')
        xticks(0:200:600)
        xticklabels(kde_x_vals(1:200:end))
        if m == 1 && i == 1
            xlabel('Time from cue on')
            ylabel('Fraction units')
        end
        ylim([0 0.7])
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',18, ...
            'fontname', 'Helvetica',  ...
            'XColor', 'black', 'YColor', 'black')
    end
end

%% Plot cat/dog diff and EXCIT timecourse on same plots

areas_to_plot = {'anterior', 'middle', 'posterior'};
boundary = 'static'; % static or dynamic
colors = cbrewer('qual', 'Paired', 6);
light_colors = colors([1 3 5], :)*0.9;
dark_colors = colors([2 4 6], :);
excit_supprn_timecourse_xvals = -50:10:430;

% overlay post timecourses of different arrays to show kinetics
figure('Position', [200 200 380 470])
hold on
for m = 1:length(KDE)
    
    % get sessions to use
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
        %         sessions_to_plot = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
        %         sessions_to_plot = [1 2 7 9];
        sessions_to_plot = [9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
        %         sessions_to_plot = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
        %         sessions_to_plot = [1 2 6 7];
        sessions_to_plot = [7];
    end
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        
        % subplots of monks x sessions
        subplot(length(KDE), length(sessions_to_plot), length(sessions_to_plot)*(m-1) + mod(i-1, length(sessions_to_plot)) + 1)
        hold on
        
        for a = 1:length(areas_to_plot)
            area = areas_to_plot{a};
            
            % get correct units
            switch area
                case 'te'
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area));
                otherwise
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, area));
            end
            
            
            % EXCIT plot
            % =====
            % get data
            heatmap_mat = zeros(length(units_to_plot), length(excit_supprn_timecourse_xvals));
            for j = 1:length(units_to_plot)
                unit = units_to_plot(j);
                heatmap_mat(j, [Monkeys(m).Sessions(sessn).UnitInfo(unit).WindowResults]==2) = 1;
            end
            %             plot(normalize(mean(heatmap_mat), 'range'), 'Color', colors(a,:), ...
            %                 'DisplayName', areas_to_plot{a}, 'LineWidth', 1.5)
            plot(excit_supprn_timecourse_xvals, mean(heatmap_mat), 'Color', light_colors(a,:), ...
                'DisplayName', areas_to_plot{a}, 'LineWidth', 1.5)
            
            % CAT/DOG plot
            % =====
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
            %             plot(normalize(mean(heatmap_mat), 'range'), 'Color', colors(a,:), ...
            %                 'DisplayName', areas_to_plot{a}, 'LineWidth', 1.5)
            plot(kde_x_vals, mean(heatmap_mat), 'Color', dark_colors(a,:), ...
                'DisplayName', areas_to_plot{a}, 'LineWidth', 1.5)
            
            
        end
        %         legend
        %         make_custom_patch_legend(dark_colors, areas_to_plot)
        title(sprintf('%s', KDE(m).Sessions(sessn).ShortName), 'Interpreter', 'none')
        if m == 1 && i == 1
            xlabel('Time from cue on')
            ylabel('Fraction units')
        end
        ylim([0 1])
        xlim([-50 450])
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XGrid', 'on', 'YGrid', 'off',...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',18, ...
            'fontname', 'Helvetica',  ...
            'XColor', 'black', 'YColor', 'black')
    end
end

%% ***** Old Code *****

%% (messy) Plot histograms of cat/dog and excit. latencies, with arrays overlaid (kinetics)

areas_to_plot = {'anterior', 'middle', 'posterior'};
boundary = 'static'; % static or dynamic
colors = cbrewer('qual', 'Set2', 3);
hist_bins = 0:25:700;

% overlay post timecourses of different arrays to show kinetics
figure('Position', [200 200 1200 400])
hold on

for m = 1:length(KDE)
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
        %         sessions_to_plot = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_plot = [1 2 7 9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
        %         sessions_to_plot = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
        sessions_to_plot = [1 2 6 7];
    end
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        subplot(length(KDE), length(sessions_to_plot), length(sessions_to_plot)*(m-1) + mod(i-1, length(sessions_to_plot)) + 1)
        hold on
        
        for a = 1:length(areas_to_plot)
            area = areas_to_plot{a};
            
            % get correct units
            switch area
                case 'te'
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area));
                otherwise
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, area));
            end
            
            % CAT/DOG DIFFS
            % =============
            % nb, these are in relative time (ie 0 to 700, but the real
            % time is -200 to 500).
            % This gets corrected with hist_bins (0:25:700) and by
            % adjusting the xticks.
            % get correct field names
            switch boundary
                case 'static'
                    exceedid = 'ExceedBD_FixedBound';
                    latid = 'ExceedBD_FixedBound_Latency';
                case 'dynamic'
                    exceedid = 'ExceedBD_DynamicBound';
                    latid = 'ExceedBD_DynamicBound_Latency';
                otherwise
                    error('Unexpected value for boundary')
            end
            
            % get data
            data = zeros(1, length(units_to_plot));
            for j = 1:length(units_to_plot)
                unit = units_to_plot(j);
                data(j) = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(latid);
            end
            hist_lines = histcounts(data, 'BinEdges', hist_bins);
            % plot
            plot(hist_lines, 'Color', colors(a,:)*1.6 - 0.6, 'LineWidth', 1.5, 'DisplayName', area)
            
            % EXCITATORY RESPONSE
            % =============
            % nb, these go -50 to 450, not -200 to 500, and are in
            % real-time already.
            % This gets adjusted with hist_bins (-200:25:500) and then
            % xticks.
            % find excitatory timestamps for marked units
            summary_nums = [Monkeys(m).Sessions(i).UnitInfo(units_to_plot).WindowResults_summary];
            excit_times = [Monkeys(m).Sessions(i).UnitInfo(units_to_plot(summary_nums == 2)).WindowResults_time];
            se_times= vertcat(Monkeys(m).Sessions(i).UnitInfo(units_to_plot(summary_nums == 3)).WindowResults_time);
            es_times = vertcat(Monkeys(m).Sessions(i).UnitInfo(units_to_plot(summary_nums == 4)).WindowResults_time);
            if isempty(se_times) && isempty(es_times)
                % do nothing
            elseif isempty(se_times)
                excit_times = [excit_times es_times(:,1)'];
            elseif isempty(es_times)
                excit_times = [excit_times se_times(:,2)'];
            else
                excit_times = [excit_times se_times(:,2)' es_times(:,1)'];
            end
            
            % plot data
            hist_line = histcounts(excit_times, 'BinEdges', hist_bins-200); % adjust hist bins to real time
            %             plot(hist_line, 'Color', colors(a,:), 'LineWidth', 1.5)
            
        end
        seventh = length(hist_lines)/7;
        xticks([floor(2*seventh) floor(4*seventh) floor(6*seventh)])
        xticklabels([0 200 400])
        if m == 1 && i == 1
            xlabel(sprintf('First string of %d signf. cat/dog diffs. (ms after cue on)', num_consec_needed))
            ylabel('Number of units')
        end
        if m == 1
            ylim([0 20])
            yticks([0 20])
        elseif m == 2
            ylim([0 15])
            yticks([0 15])
        end
        title(sprintf('%s', KDE(m).Sessions(sessn).ShortName))
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',18, ...
            'fontname', 'Helvetica',  ...
            'XColor', 'black', 'YColor', 'black')
    end
end

%% (uninformative) Bar graph of excit vs cat/dog latencies, by array

areas_to_plot = {'anterior', 'middle', 'posterior'};
boundary = 'static'; % static or dynamic
colors = cbrewer('qual', 'Set2', 3);
hist_bins = 0:25:700; % relative time for cat/dog diffs; need to adjsut for excit/supprn data.

matlab_colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3;
    0          0        0    ];


for m = 1:length(KDE)
    
    % prep figure
    figure('Position', [200 200 600 300])
    hold on
    
    % get sessions to plot
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
        %         sessions_to_plot = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
        %         sessions_to_plot = [1 2 7 9];
        sessions_to_plot = 9;
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
        %         sessions_to_plot = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
        %         sessions_to_plot = [1 2 6 7];
        sessions_to_plot = 7;
    end
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        %         subplot(length(KDE), length(sessions_to_plot), length(sessions_to_plot)*(m-1) + mod(i-1, length(sessions_to_plot)) + 1)
        %         subplot(2,1,m)
        %         hold on
        
        % pre allocate
        bar_mat_means = zeros(length(areas_to_plot), 2); % (num arrays) x 2. Col 1 excit, col 2 cat/dog.
        bar_mat_stds = zeros(length(areas_to_plot), 2);
        
        for a = 1:length(areas_to_plot)
            area = areas_to_plot{a};
            
            % get correct units
            switch area
                case 'te'
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, area));
                otherwise
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, area));
            end
            
            % CAT/DOG DIFFS
            % =============
            % nb, these are in relative time (ie 0 to 700, but the real
            % time is -200 to 500).
            % This gets corrected with hist_bins (0:25:700) and by
            % adjusting the xticks.
            % get correct field names
            switch boundary
                case 'static'
                    exceedid = 'ExceedBD_FixedBound';
                    latid = 'ExceedBD_FixedBound_Latency';
                case 'dynamic'
                    exceedid = 'ExceedBD_DynamicBound';
                    latid = 'ExceedBD_DynamicBound_Latency';
                otherwise
                    error('Unexpected value for boundary')
            end
            
            % get data
            data = zeros(1, length(units_to_plot));
            for j = 1:length(units_to_plot)
                unit = units_to_plot(j);
                data(j) = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(latid);
            end
            bar_mat_means(a, 2) = mean(data, 'omitnan');
            bar_mat_stds(a, 2) = std(data, 'omitnan');
            
            % EXCITATORY RESPONSE
            % =============
            % nb, these go -50 to 450, not -200 to 500, and are in
            % real-time already.
            % This gets adjusted with hist_bins (-200:25:500) and then
            % xticks.
            % find excitatory timestamps for marked units
            summary_nums = [Monkeys(m).Sessions(i).UnitInfo(units_to_plot).WindowResults_summary];
            excit_times = [Monkeys(m).Sessions(i).UnitInfo(units_to_plot(summary_nums == 2)).WindowResults_time];
            se_times= vertcat(Monkeys(m).Sessions(i).UnitInfo(units_to_plot(summary_nums == 3)).WindowResults_time);
            es_times = vertcat(Monkeys(m).Sessions(i).UnitInfo(units_to_plot(summary_nums == 4)).WindowResults_time);
            if isempty(se_times) && isempty(es_times)
                % do nothing
            elseif isempty(se_times)
                excit_times = [excit_times es_times(:,1)'];
            elseif isempty(es_times)
                excit_times = [excit_times se_times(:,2)'];
            else
                excit_times = [excit_times se_times(:,2)' es_times(:,1)'];
            end
            bar_mat_means(a, 1) = mean(excit_times, 'omitnan');
            bar_mat_stds(a, 1) = std(excit_times, 'omitnan');
            
        end
        b = bar(bar_mat_means, 'grouped');
        pause(0.5)
        xoffs = vertcat(b.XOffset);
        xvals = [1 2 3; 1 2 3];
        xvals = xvals + xoffs;
        errorbar(xvals', bar_mat_means, bar_mat_stds, 'o', 'Color', 'k')
        xticks(1:length(areas_to_plot))
        xticklabels(areas_to_plot)
        ylabel('Latency (ms after cue on)')
        ylim([0 600])
        title(sprintf('%s', KDE(m).Sessions(sessn).ShortName))
        make_custom_patch_legend(matlab_colors([1 2],:), {'Excitatory latency', 'Cat/dog difference latency'})
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',18, ...
            'fontname', 'Helvetica',  ...
            'XColor', 'black', 'YColor', 'black')
        pause(0.5)
    end
end

%% (old) Look at start times and number of runs

% the change here is more subtle...just looks like absolute number of long
% runs has gone up. I doubt the time profile has actually changed.
% Should probably break down heatmaps by array for base1, base2, pre and
% post. Could test proportion of units showing a cat/dog diff for atleast
% XYZ ms -- maybe 20, considering that was the bandwidth?

for m = 1:length(KDE)
    sessions_to_plot = 1:length(KDE(m).Sessions);
    %     sessions_to_plot = [1];
    figure('Position', [200 200 1000 1000])
    hold on
    sq = ceil(sqrt(length(sessions_to_plot)));
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        runs =  KDE(m).Sessions(sessn).KDE_EBD_Runs;
        
        % concat
        data = vertcat(runs{:,1});
        
        % plot histograms
        subplot(sq, sq, i)
        hold on
        title(sprintf('%s, \n median +/- iqr, %0.2f +/- %0.2f', KDE(m).Sessions(sessn).ShortName, median(data), iqr(data)))
        histogram(data, 'BinEdges', 0:25:700)
    end
end

%% (old) Look at distribution of lengths of time dogs/cats differ

% Distributions are at least bi-modal, if not tri or quatra.
% So ANOVAs are out for comparing. Post learning has longer for sure, but
% effect is weak -- probably looking at non-parametric test.

figure('Position', [200 200 1000 1000])
hold on
ranksum_alpha = 0.01;

for m = 1:length(KDE)
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
        %         sessions_to_use = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_plot = [1 2 7 9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
        %         sessions_to_use = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
        sessions_to_plot = [1 2 6 7];
    end
    
    %     figure % for each mk separately, all sessns
    %     hold on
    sq = ceil(sqrt(length(sessions_to_plot)));
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        subplot(length(KDE), length(sessions_to_plot), length(sessions_to_plot)*(m-1) + mod(i-1, length(sessions_to_plot)) + 1) % both mks together, selected sessions
        hold on
        %         units_to_plot = 1:length(KDE(m).Sessions(sessn).UnitInfo);
        units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, 'te'));
        runs = cell(length(units_to_plot), 2); % runs of significant KDE difference. First col is start time in ms of run, second col is length of run.
        
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            
            % get data
            exceed_bd_inds = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.ExceedBD;
            actual_xvals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals;
            max_val = length(actual_xvals);
            
            % catch issues
            if isempty(exceed_bd_inds)
                runs{j,1} = 0;
                runs{j,2} = 0;
                continue
            end
            
            % convert to string of zeros and ones
            exceed_bd_bools = zeros(1, max_val);
            exceed_bd_bools(exceed_bd_inds) = 1;
            s = sprintf('%d', exceed_bd_bools);
            
            % here goes some matlab magic
            runs{j,1} = (find(diff(uint8(s)))+1)'; % inds of first 1 in each run of 1's
            s = textscan(s, '%s', 'delimiter', '0', 'multipleDelimsAsOne',1); % all runs of 1s
            s = s{:}; % unpack
            runs{j,2} = cellfun('length', s); % length of each run
        end
        
        % concat
        data = vertcat(runs{:,2});
        
        % plot histograms
        %         subplot(sq, sq, i) % for each mk separately, all sessns
        %         hold on
        title(sprintf('%s, \n %0.2f +/- %0.2f', KDE(m).Sessions(sessn).ShortName, median(data), iqr(data))) % both mks together, selected sessions
        histogram(data, 'BinEdges', 0:5:500, 'Normalization', 'probability')
        %         set(gca, 'YScale', 'log')
        %         set(gca, 'XScale', 'log')
        
        % store data
        KDE(m).Sessions(sessn).KDE_EBD_Runs = runs;
        if strcmp('Post', regexp(KDE(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once'))
            post_data = data;
        end
        
        % format plot
        ylim([0 0.3])
        if i == 1 && m == 1
            xlabel('Length of run')
            ylabel({'Probability across all runs','(Not norm''d to num unts)'})
        end
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',18, ...
            'fontname', 'Helvetica', 'fontweight', 'normal', ...
            'XColor', 'black', 'YColor', 'black')
    end
    
    % stats testing
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        if strcmp('Post', regexp(KDE(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once'))
            continue
        end
        runs = KDE(m).Sessions(sessn).KDE_EBD_Runs;
        data = vertcat(runs{:,2});
        
        [p,h] = ranksum(data, post_data, 'Alpha', ranksum_alpha);
        fprintf('%s, session %d: h = %d, p = %0.3f \n', KDE(m).Name, sessn, h, p)
    end
end
sgtitle('Session, median +/- iqr')
set(gcf, 'Color', 'white')


