% uses KDEs to estimate visual and catg. latencies

%% Load temporal activation/suppression data
clearvars
close all
path1 = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = 'MaxMarta_xma2_temporalActvnSuppn.mat';
Data = load(fullfile(path1, fname)); 
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Load raw KDEs, transfer, remove
bw = 20;
path1 = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = sprintf('MaxMarta_xma2_dualKDEs_bw%d.mat', bw);
Data = load(fullfile(path1, fname));
[status, KDE] = stitch_monkeyStruct_from_parts(Data);
clear Data
clearvars -except Monkeys KDE bw

for m = 1:length(Monkeys)
    for i = 1:length(Monkeys(m).Sessions)
        for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
            Monkeys(m).Sessions(i).UnitInfo(j).CueOnAllCues = KDE(m).Sessions(i).UnitInfo(j).CueOnAllCues;
        end
    end
end
clear KDE

%% Load bootstrapped catg data, get bootstrap diffs, transfer, remove
bw = 20;
path1 = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = sprintf('MaxMarta_fix_cat_xma2_bootstrappedKDEs_BW%d.mat', bw);
Data = load(fullfile(path1, fname));
[status, KDE] = stitch_monkeyStruct_from_parts(Data);
clear Data

pct = 99;
diffid = sprintf('TrueDiff_BW%d', bw);
bootid = sprintf('BootstrappedDiffs_BW%d', bw);

for m = 1:length(KDE)
       
    % skip Marta session 4 (bad data)
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
        sessions_to_plot = [1:3 5:9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
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
            
            % get 99th percentile of diffs
            percentile_diff = prctile(boot_diffs, pct);
            
            % field name
            exceedid = 'ExceedBD_FixedBound';
            
            % store
            Monkeys(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid) = find(abs(true_diffs) > percentile_diff); % normalized times in ms where true diffs exceeds bootstrapped diffs
            Monkeys(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid) = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid);
            
        end
    end
end

clear KDE


% moving boundary version
% pct = 99;
% diffid = sprintf('TrueDiff_BW%d', bw);
% bootid = sprintf('BootstrappedDiffs_BW%d', bw);
% 
% for m = 1:length(Monkeys)
%     % for m = 1
%     sessions_to_plot = 1:length(Monkeys(m).Sessions);
%     %     sessions_to_plot = [7 9];
%     Monkeys(m).Name = Monkeys(m).Name;
%     
%     
%     for i = 1:length(sessions_to_plot)
%         sessn = sessions_to_plot(i);
%         units_to_plot = 1:length(Monkeys(m).Sessions(sessn).UnitInfo);
%         
%         for j = 1:length(units_to_plot)
%             unit = units_to_plot(j);
%             
%             kde_x_vals = Monkeys(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals;
%             true_diffs = Monkeys(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid);
%             boot_diffs = Monkeys(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid);
%             boot_diffs = abs(boot_diffs); % just look at magnitude, not sign
%             percentile_diffs = prctile(boot_diffs, pct, 1); % 1 x (num xvals) vector of the 95th percentile of bootstrapped data
%             exceedid = 'ExceedBD_DynamicBound';
%             Monkeys(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid) = find(abs(true_diffs) > percentile_diffs); % normalized times in ms where true diffs exceeds bootstrapped diffs
%             
%         end
%     end
% end

%% Get short names, add area labels
marta_xls = '/Users/jonahpearl/Documents/MATLAB/matsumoto/RecordingMarta.xlsx';
max_xls = '/Users/jonahpearl/Documents/MATLAB/matsumoto/RecordingMax.xlsx';

for m = 1:length(Monkeys)
    
    % get short session names
    date_strs = {Monkeys(m).Sessions.DateStr};
    if regexp(Monkeys(m).Name, 'Marta\w*')
        short_names = get_short_names(date_strs, marta_xls, '\w*cats\w*');
    else
        short_names = get_short_names(date_strs, max_xls, '\w*cats\w*');
    end
    for i = 1:length(Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).ShortName = short_names{i};
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

%% Add CueInfo, get number of cat/dog trials, then remove CueInfo
path1 = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = 'MaxMarta_xma2_raw_behav.mat';
catgs = {1:260, 261:520};

% load data
Data = load(fullfile(path1, fname));
Data = Data.Monkeys;

for m = 1:length(Monkeys)
    for i = 1:length(Monkeys(m).Sessions)
        
        % error check
        if ~strcmp(Data(m).Sessions(i).SessionName, Monkeys(m).Sessions(i).SessionName)
            error('Session names do not match')
        else % move into KDE struct
            num_total_cat_trials = sum([Data(m).Sessions(i).CueInfo(catgs{1}).NumApp]);
            num_total_dog_trials = sum([Data(m).Sessions(i).CueInfo(catgs{2}).NumApp]);
            Monkeys(m).Sessions(i).NumCatgTrials = [num_total_cat_trials num_total_dog_trials];
        end
    end
end
clearvars -except Monkeys bw

%% Save compiled data
path1 = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = sprintf('MaxMarta_fix_cat_xma2_latencies_and_KDEs_BW%d.mat', bw);
[status, vnames] = split_monkeyStruct_in_parts(Monkeys);
save(fullfile(path1, fname), vnames{:});

%% Reload compiled data
clearvars
close all
bw = 20;
path1 = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = sprintf('MaxMarta_fix_cat_xma2_latencies_and_KDEs_BW%d.mat', bw);
Data = load(fullfile(path1, fname)); 
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Plot KDEs / catdog-diff timecourses of chosen units

% Window results summary key:
% 0:suppressed, 1:nothing, 2:exicted, 3:suppressed then excited, 4:excited then suppressed

% KDE Params
latency_bounds = [75 175]; % only plot units that respond first in this window
shape_to_plot = 'exc'; % exc (ie, exc and exc then supp), supp_only, supp_exc (ie, supp then exc), noise

% Catdog diff params
catdog_diff_curves_to_show = 'none'; % all, none, or signf (only signf cat/dog diffs)                                       % all (show all units' KDEs and catdog diffs)
num_diff_timepoints_for_signf = 40; % not necessarily consecutive bc that's a pain to code for now
bw = 20; % kernel bandwidth

% Plotting params
savePath = '/Users/jonahpearl/Documents/BJR group/Catdog paper';

% Set up figure
figure2('Position', [200 200 500 500])
hold on

% Get x values for plotting
kde_x_vals = Monkeys(1).Sessions(1).UnitInfo(1).CueOnAllCues.KDEXVals;

% Field name
diffid = sprintf('TrueDiff_BW%d', bw);

for m = 1:length(Monkeys)
    
    % Get sessions to use
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_plot = [1 2 3 5 6 7 9 10]; % skip base04 (outlier) and post01 (low trial count)
%         sessions_to_plot = [1 2 7 9];
        sessions_to_plot = 7;
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_plot = [1 2 3 5 6 7 9]; % skip base04 and sub03 (low trial count)
%         sessions_to_plot = [1 2 6 7];
        sessions_to_plot = 6;
    end
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        
        % Subplots of monks x sessions
        subplot(length(Monkeys), length(sessions_to_plot), length(sessions_to_plot)*(m-1) + mod(i-1, length(sessions_to_plot)) + 1)
        hold on
        
        % Get data
        window_results = [Monkeys(m).Sessions(sessn).UnitInfo.WindowResults_summary];
        lats = {Monkeys(m).Sessions(sessn).UnitInfo.WindowResults_time};
        
        % Subset based on params
        switch shape_to_plot
            case 'exc'
                shape_bool = (window_results == 2) | (window_results == 4); % excited first
                lats = cellfun(@(x) x(1), lats); % take first value if there are two
                within_lat_bool = (lats >= latency_bounds(1)) & (lats <= latency_bounds(2));
                
            case 'supp_only'
                shape_bool = (window_results == 0);
                lats = cellfun(@(x) x(1), lats);
                within_lat_bool = (lats >= latency_bounds(1)) & (lats <= latency_bounds(2));
                
            case 'supp_exc'
                shape_bool = (window_results == 3); % excited second
                lats = cellfun(@(x) x(end), lats); % take second value if there are two
                within_lat_bool = (lats >= latency_bounds(1)) & (lats <= latency_bounds(2));
                
            case 'noise'
                shape_bool = (window_results == 1); % noise
                lats = cellfun(@(x) x(1), lats);
                within_lat_bool = ones(1, length(window_results)); % take all
        end
        
        % Subset further only to te units
        area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, 'te');
        units_to_plot = find(shape_bool & within_lat_bool & area_bool);
        
        % Plot
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            
            % Try to get data, skip if no data (ie too few spikes)
            try
                kd = Monkeys(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_raw; % raw KDE
                catdog_diff = abs(Monkeys(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid)); % raw cat/dog diff
                exceed_bd = Monkeys(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.ExceedBD_FixedBound; % timepoints where catdog diff is signf
            catch
                continue
            end
            
            % Scale to spikes / second / trial:
            % (density / ms) x (1000 ms / s) x (num spikes) / (num trials)
            num_spikes = length(Monkeys(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterXVals);
            num_trials = sum(Monkeys(m).Sessions(sessn).NumCatgTrials);
            scaled_kd = kd * 1000 * num_spikes / num_trials;
            
            % Plot kde
%             yyaxis left
            plot(kde_x_vals, normalize(kd, 'range'), '-', 'Color', mlc(1), 'LineWidth', 0.5)
%             plot(kde_x_vals, scaled_kd, '-', 'Color', mlc(1), 'LineWidth', 0.5)
            
            % Plot catdog diffs if requested
%             yyaxis right
            if strcmp(catdog_diff_curves_to_show, 'all')
                plot(kde_x_vals, catdog_diff, '-', 'Color', mlc(2), 'LineWidth', 0.5)
            elseif strcmp(catdog_diff_curves_to_show, 'signf') && (length(exceed_bd) > num_diff_timepoints_for_signf)
                plot(kde_x_vals, catdog_diff, '-', 'Color', mlc(2), 'LineWidth', 0.5)
%                 plot(kde_x_vals, normalize(catdog_diff, 'range'), '-', 'Color', mlc(2), 'LineWidth', 0.5)
            elseif strcmp(catdog_diff_curves_to_show, 'none')
                % do nothing
            end
        end
        
        % Format plot
%         yyaxis left
        yticks([0 0.5 1])
%         yyaxis right
%         yl = get(gca, 'YLim');
%         yticks([0 yl(2)/2 yl(2)])
        title(sprintf('%s (Mean exc. lat. %0.0f ms , n = %d)',...
            Monkeys(m).Sessions(sessn).ShortName, ...
            mean(lats(units_to_plot), 'omitnan'), length(units_to_plot)),...
            'Interpreter', 'none')
        sgtitle(sprintf('Response type %s, latency within interval %d to %d', shape_to_plot, latency_bounds(1), latency_bounds(2)), ...
             'Interpreter', 'none')
         formatSVMPlot(gca, gcf)
         saveas(gcf, fullfile(savePath, sprintf('%s_%d_to_%d', ...
             shape_to_plot, latency_bounds(1), latency_bounds(2))), 'epsc')
    end
end

%% Functions
function formatSVMPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',28, 'YGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end