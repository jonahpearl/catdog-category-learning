% figs for the cat/dog manuscript, 18 Oct 2023

% to run the bootstrapping, see: /Volumes/Alex's Mac Backup/Documents/MATLAB/matsumoto/jep_fix_cat_xma2_bootstrapKDEs.m
%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs'; % pv --> "passive viewing"

% Pick KDE bandwidth parameter and load
bw = 20;
fname = sprintf('MaxMarta_fix_cat_xma2_bootstrappedKDEs_BW%d.mat', bw);
Data = load(fullfile(EXT_HD, pv_path, fname));
[status, KDE] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Misc prep

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
marta_xls = fullfile(EXT_HD, '/RecordingMarta.xlsx');
max_xls = fullfile(EXT_HD, '/RecordingMax.xlsx');

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

xtickcolors = [mlc(1:2); 0.5 0.5 0.5];
for m = 1:length(KDE)
    if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
%         KDE(m).Sessions_to_use = [1 2 3 5 6 7 9];
%         KDE(m).XTickLabs= {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
%         KDE(m).Sessions_to_use = [1 2 7 9];
%         KDE(m).XTickLabs= {'Base 1', 'Base 2', 'Pre', 'Post'};
        KDE(m).Sessions_to_use = [7 9];
        KDE(m).XTickLabs= {'Pre', 'Post'};

%         KDE(m).Sessions_to_use = [7 9];
%         KDE(m).XTickLabs= {'Pre', 'Post'};
        KDE(m).Code = 'R';
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
%         KDE(m).Sessions_to_use = [1 2 3 4 5 6 7];
%         KDE(m).XTickLabs = {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
%         KDE(m).Sessions_to_use = [1 2 6 7];
%         KDE(m).XTickLabs = {'Base 1', 'Base 2', 'Pre', 'Post'};
        KDE(m).Sessions_to_use = [6 7];
        KDE(m).XTickLabs = {'Pre', 'Post'};
        KDE(m).Code = 'X';
    end
end

%% Fig 4B params

rArrayLocs = {'anterior', 'middle', 'posterior'};
chi_sq_alpha = 0.05;
colors = [0.5 0.5 0.5; mlc(1:2)];
range_normalize = false;
rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)
boundary = 'dynamic'; % static or dynamic boundary for bootstrapping

%% Fig 4 prep: get xvals where true diffs are larger than bootstrapped popn (dynamic boundary) 

diffid = sprintf('TrueDiff_BW%d', bw);
bootid = sprintf('BootstrappedDiffs_BW%d', bw);
pct = 99;
for m = 1:length(KDE)
    sessions_to_plot = rSessionsByMonk{m};
    
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
            percentile_diffs = prctile(boot_diffs, pct, 1); % 1 x (num xvals) vector of the 95th percentile of bootstrapped data
            exceedid = 'ExceedBD_DynamicBound';
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid) = find(abs(true_diffs) > percentile_diffs); % normalized times in ms where true diffs exceeds bootstrapped diffs
            
        end
    end
end

%% Fig 4B: plot results as heatmap

nrows_insert = 3; % size of black line separating arrays
sort_by = 'after_zero_latency'; % after_zero_latency or duration
sort_empty_val = 1000; % controls if units with no signf diff go above or below others
divider_color = [0.5 0.5 0.5 ];
save_dir = '/Users/jonahpearl/Documents/BJR group/Catdog_paper/Feb_2023_addtl_figs/kde_diff_heatmaps';

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
        switch regexp(KDE(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
            case 'Base'
                col = xtickcolors(1,:);
            case 'Pre'
                col = xtickcolors(2,:);
            case 'Post'
                col = xtickcolors(3,:);
        end
        
        % Take only te units
        area_bool = strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, 'te');
        units_to_plot = find(area_bool);
        
        % Sort units by array
        array_list = {KDE(m).Sessions(sessn).UnitInfo(area_bool).Location};
        [array_list, idx] = sort(array_list);
        units_to_plot = units_to_plot(idx);
        
        % Get data (already sorted by array)
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
        xticks(0:200:700)
        xticklabels(-200:200:500)
%         xticks([0 400])
%         xticklabels([0 400])
        ylim([0 sum(max_num_units_per_array) + 6])
        yticks(ytick_spots)
        yticklabels(ytick_labs)
        ytickangle(90)
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

%% Fig 4C: plot avg timecourses (cat/dog diff) with days overlaid 

% Pre allocate
pre_data = cell(length(KDE), length(rArrayLocs));
post_data = cell(length(KDE), length(rArrayLocs));
chi_sq_inds = 200:50:700; % inds in kde_x_vals to test...corresponds to -200:50:500 in real time.

yl_by_monk = {[0 0.7], [0 0.4]};

for m = 1:length(KDE)
    
    % Get sessions to use
    sessions_to_plot = KDE(m).Sessions_to_use;
%     if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_plot = [7 9];
%     elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_plot = [6 7];
%     end
    
    
    % Set up figure
    figure2('Position', [200 200 1200 330]) % for arrays
%     figure2('Position', [400 400 500 300]) % just te
    hold on
    
    for a = 1:length(rArrayLocs)
        area = rArrayLocs{a};
        
                % Plot all in one figure2:
%                 switch area
%                     case 'te'
%                         subplot(2,3,1:3)
%                     otherwise
%                         subplot(2,3,2+a)
%                 end
%                 hold on
                

                %TE alone: don't do any subplots
                
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
            if ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Base', 'once'))
                col_ind = 1;
            elseif ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Pre', 'once'))
                col_ind = 2;
            elseif ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Post', 'once'))
                col_ind = 3;
            end
            
            % Get real x-axis values
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(1).CueOnAllCues.KDEXVals;
            
            % Collect data
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
            
            % Calculate mean
            heatmap_mean = mean(heatmap_mat); % 1 x num xvals
            if range_normalize
                heatmap_mean = normalize(heatmap_mean, 'range');
            end
            
            % Set xvals to use
            xvals = 1:length(heatmap_mean);
            
            
            % Plot the data
            plot(xvals, heatmap_mean, ...
                'LineWidth', 1.5, ...
                'DisplayName', KDE(m).Sessions(sessn).ShortName,...
                'Color', colors(col_ind, :))
            
            % Format plot
%             title(sprintf('Monkey %s, %s', KDE(m).Code, area), 'Interpreter', 'none')
            title(sprintf('%s', area), 'Interpreter', 'none')
            xticks(0:200:600)
            xticklabels(kde_x_vals(1:200:end))
            
            if strcmp(area, 'te')
                xlabel('Time from cue on')
                ylabel('Fraction units')
            end
            
            % Store pre / post data for stats testing
            if ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Post', 'once'))
                post_data{m,a} = heatmap_mat;
            elseif ~isempty(regexp(KDE(m).Sessions(sessn).ShortName, 'Pre', 'once'))
                pre_data{m,a} = heatmap_mat;
            end
        end
        
        % Do stats testing
        pre = pre_data{m,a};
        post = post_data{m,a};
        yl = get(gca, 'YLim');
        for ii = chi_sq_inds
            [h,p] = prop_test([sum(pre(:,ii)) sum(post(:,ii))], [size(pre,1) size(post,1)], false, chi_sq_alpha);
            fprintf('%s, %s, %d \n', KDE(m).Sessions(sessn).ShortName, area, p)
            if h
                scatter(ii, yl(2)*0.97, 'ko', 'filled')
            end
        end

        % Format the plot some more
        if a == 1
            ylabel('Proportion signf. units')
            xlabel('Time from cue on')
        end
        ylim(yl_by_monk{m})
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',...
            'fontsize',26, ...
            'fontname', 'Helvetica', ...
            'XColor', 'black', 'YColor', 'black')
        
    end
    % Format the plot some more
%     sgtitle(sprintf('Monkey %s', KDE(m).Code))
end

%% Fig 4D: duration of signf cat/dog runs (pre/post overlay)

% Parameters
ranksum_alpha = 0.05;
areas_to_plot = {'anterior', 'middle', 'posterior', 'te'};
figs_path = '/Users/jonahpearl/Documents/BJR group/Catdog paper';

for a = 1:length(areas_to_plot)
    area = areas_to_plot{a};
    
    % Plot both mks together
%     figure2('Position', [200 200 1000 1000])
    figure2('Position', [200 200 300 300])
    hold on
    
    for m = 1:length(KDE)
        sessions_to_plot = KDE(m).Sessions_to_use;
        % Arrange subplots
        subplot(2,1,m)
        hold on

        for i = 1:length(sessions_to_plot)
            sessn = sessions_to_plot(i);
            switch area
                case 'te'
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Area}, 'te'));
                otherwise
                    units_to_plot = find(strcmp({KDE(m).Sessions(sessn).UnitInfo.Location}, area));
            end
            runs = cell(length(units_to_plot), 2); % runs of significant KDE difference. First col is start time in ms of run, second col is length of run.

            for j = 1:length(units_to_plot)
                unit = units_to_plot(j);

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

            % Plot the data as a histogram        
            histogram(data, 'BinEdges', [0:25:500 Inf], ... % 500
                'Normalization', 'probability', ...
                'DisplayName', KDE(m).Sessions(sessn).ShortName, ...
                'FaceColor', mlc(i))

            % Add mean/std above the histogram
            yl = ylim();
            mu = mean(data);
            sigma = std(data);
            scatter(mu, yl(2)+0.05, 100, '^', 'MarkerFaceColor', mlc(i), 'MarkerEdgeColor', mlc(i))
            plot([mu - sigma, mu+sigma], repelem(yl(2)+0.05,2), ...
                'LineWidth', 2, 'Color', mlc(i), 'MarkerFaceColor', mlc(i))
            
            % Store the data
            run_id = sprintf('KDE_EBD_Runs_%s', area); 
            KDE(m).Sessions(sessn).(run_id)= runs;
            if strcmp('Post', regexp(KDE(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once'))
                post_data = data;
            end

            % Format the plot
%             title(sprintf('Lengths of signf. runs, %s, session %s, %s',...
%                 KDE(m).Name, KDE(m).Sessions(sessn).ShortName, area), 'Interpreter', 'none')
            title(sprintf('%s, %s', KDE(m).Name, area), 'interpreter', 'none')
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
            if strcmp('Post', regexp(KDE(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once'))
                continue
            end
            run_id = sprintf('KDE_EBD_Runs_%s', area); 
            runs = KDE(m).Sessions(sessn).(run_id);
            data = vertcat(runs{:,2});
            if isempty(data)
                fprintf('%s, session %d, area %s: no signf timepoints \n', KDE(m).Name, sessn, area, h, p)
                continue
            end
            data(data == 0) = [];  % remove zeros
            [p,h] = ranksum(data, post_data, 'Alpha', ranksum_alpha, 'tail', 'left');
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