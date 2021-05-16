%analyzes population-level, non-visual-related statistics of neural recordings, eg number of units, locations of units. Also VR.

%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';
figureSavePath = '/Users/jonahpearl/Documents/BJR group/Catdog paper/';
fname = 'MaxMarta_xma2_ni.mat';
Data = load(fullfile(EXT_HD, pv_path, fname));
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

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

%% Pool posteriors together and get num of arrays
num_arrays = [];
for m = 1:length(Monkeys)
    num_arrays_sessn = [];
    for i = 1:length(Monkeys(m).Sessions)
        arrays = unique({Monkeys(m).Sessions(i).UnitInfo.Location});
        num_arrays_sessn(i) = numel(arrays);
        for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
            if regexp(Monkeys(m).Sessions(i).UnitInfo(j).Location, 'posterior\d')
                Monkeys(m).Sessions(i).UnitInfo(j).Location = 'posterior';
            end
        end
    end
    num_arrays(m) = max(num_arrays_sessn);
end

%% (Marta) Get Utah locations

% from Narihisa: elecnum maps the "electrode number" in the NEV files onto
% the real position on the Utah array.

% "Real position" in the Utah array means:
%  x  ..  28  18  x
%  96 ..  27  17  8
%  .. ..  ..  ..  ..
%  90 ..  21  11  2
%  89 ..  20  10  1
%  x  ..  19  9   x

% elecnum("electrode number") = real position.
% EG, elecnum("2") = 1.
elecnum = [17,1,18,2,19,3,20,4,21,5,22,6,23,7,24,8,25,9,26,10,27,11,28,12,29,13,30,14,31,15,32,16,49,33,50,34,51,35,52,36,53,37,54,38,55,39,56,40,57,41,58,42,59,43,60,44,61,45,62,46,63,47,64,48,81,65,82,66,83,67,84,68,85,69,86,70,87,71,88,72,89,73,90,74,91,75,92,76,93,77,94,78,95,79,96,80];
elecnum_Marta = horzcat(elecnum,elecnum,elecnum);
elecnum_A = horzcat(elecnum,elecnum(1:32),elecnum,elecnum(33:64));  %NSP1 (anterior, posterior 1:32, teo, posterior 33:64)
elecnum_B = horzcat(elecnum,elecnum(65:96)); %NSP2 (middle, posterior 65:96)

for m = 1:length(Monkeys)
    for i = 1:length(Monkeys(m).Sessions)
        for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
            if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2') 
               Monkeys(m).Sessions(i).UnitInfo(j).Utah_Location = elecnum_Marta(Monkeys(m).Sessions(i).UnitInfo(j).Electrode); 
                
            elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
                if strcmp('a', Monkeys(m).Sessions(i).UnitInfo(j).NEV)
                    Monkeys(m).Sessions(i).UnitInfo(j).Utah_Location = elecnum_A(Monkeys(m).Sessions(i).UnitInfo(j).Electrode);
                elseif strcmp('b', Monkeys(m).Sessions(i).UnitInfo(j).NEV)
                   Monkeys(m).Sessions(i).UnitInfo(j).Utah_Location = elecnum_B(Monkeys(m).Sessions(i).UnitInfo(j).Electrode); 
                end
                
            end
        end
    end
end

%% Set sessions to use and xtick names and colors
xtickcolors = cbrewer('qual', 'Dark2', 3);
xtickcolors = xtickcolors([3 2 1], :);
for m = 1:length(Monkeys)
    if strcmp(Monkeys(m).Name, 'Marta_fix_cat_xma2')
%         Monkeys(m).Sessions_to_use = [1 2 3 5 6 7 9];
%         Monkeys(m).XTickLabs= {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};
        Monkeys(m).Sessions_to_use = [7 9];
        Monkeys(m).XTickLabs= {'Pre', 'Post'};

        Monkeys(m).Code = 'R';
    elseif strcmp(Monkeys(m).Name, 'Max_fix_cat_xma2')
        Monkeys(m).Sessions_to_use = [6 7];
        Monkeys(m).XTickLabs = {'Pre', 'Post'};
%         Monkeys(m).Sessions_to_use = 1:7;
%         Monkeys(m).XTickLabs = {'Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5', 'Pre', 'Post'};

        Monkeys(m).Code = 'X';
    end
end

%% Collect spike counts
catgs = {1:260, 261:520};

% intervals_to_test = {[-150 -50], [75 175], [175 275], [275 375]}; % time window from cue on
intervals_to_test = {[-175 0], [-150 -50], [175 350], [75 175], [175 275]};
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
        fprintf('Done with %s session %d \n', Monkeys(m).Name, sessn)
    end
end

%% Get visual responsiveness

% test_intervals = {[75 175], [175 275], [275 375], [75 375]};
% baseline_intervals = {[-150 -50], [-150 -50], [-150 -50], [-250 50]};
% test_intervals = {[75 175], [175 350], [175 275]};
% baseline_intervals = {[-150 -50], [-175 0], [-150 -50]};


by = 'all'; % other option is by image
switch by
    case 'all'
        groupings = []; % overall vis resp
    case 'category'
        groupings = [1 2]; % cat vs dog vis resp
    case 'image'
        groupings = 1:520; % vis resp for each image (count of this gives rough image-level selectivity)
end

for m = 1:length(Monkeys)
    for i = 1:length(Monkeys(m).Sessions)
        
        
        for k = 1:length(test_intervals)
            
            test_int = test_intervals{k};
            baseline_int= baseline_intervals{k};
            xid = get_good_interval_name2(test_int, 'full', 'X');
            yid = get_good_interval_name2(baseline_int, 'full', 'X');
            
            % get data
            X1 = Monkeys(m).Sessions(i).(xid); % test int
            X2 = Monkeys(m).Sessions(i).(yid); % baseline int

            switch by
                case 'all'
                    tid = get_good_interval_name2(test_int, 'full', 'VisResp_all_test');
                    bid = get_good_interval_name2(baseline_int, '', '');
                    vrID = strcat(tid,bid);

                    for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
                        Monkeys(m).Sessions(i).UnitInfo(j).(vrID) = signrank(X2(:,j) - X1(:,j), 0, 'tail', 'left'); % tests that median(baseline-test) < 0
                    end
                    
                case 'category'
                    for p = 1:length(groupings)
                        tid = get_good_interval_name2(test_int, 'full', sprintf('VisResp_test_catg_%d',p));
                        bid = get_good_interval_name2(baseline_int, '', '');
                        vrID = strcat(tid,bid);
                        idx = find(Monkeys(m).Sessions(i).Session_Y_catg == p);

                        for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
                            Monkeys(m).Sessions(i).UnitInfo(j).(vrID) = ...
                                signrank(reshape(X2(idx,j),1,numel(X2(idx,j))) - reshape(X1(idx,j),1,numel(X1(idx,j))),...
                                0, 'tail', 'left'); % tests that median(baseline-test) < 0
                        end
                    end
                    
                case 'image'
                    
                    % session 8 only has 3 presentations per img, skip it
                    if i == 8
                        continue
                    end
        
                    img_vr_mat = zeros(length(groupings), length(Monkeys(m).Sessions(i).UnitInfo));
                    
                    for p = 1:length(groupings)
                        tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img');
                        bid = get_good_interval_name2(baseline_int, '', '');
                        vrID = strcat(tid,bid);    
                        idx = find(Monkeys(m).Sessions(i).Session_Y_imageID == p);
                        
                        idx = idx(1:5); % use first five presentations
                        for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
                            img_vr_mat(p,j) = ...
                                signrank(reshape(X2(idx,j),1,numel(X2(idx,j))) - reshape(X1(idx,j),1,numel(X1(idx,j))),...
                                0, 'tail', 'left'); % tests that median(baseline-test) < 0
                        end
                    end
                    Monkeys(m).Sessions(i).(vrID) = img_vr_mat;
            end
                
                    
        end
    end
end

%% Bar plot of unit array / amount across days

% test_interval = {[75 175]};
% baseline_interval = {[-150 -50]};
test_interval = {[175 350]};
baseline_interval = {[-175 0]};
% test_interval = {[175 275]};
% baseline_interval = {[-150 -50]};
tid = get_good_interval_name2(test_int, 'full', 'VisResp_all_test');
bid = get_good_interval_name2(baseline_int, '', '');
vrID = strcat(tid,bid);
vr_alpha = 0.001;


figure
hold on

for m = 1:length(Monkeys)
    
    % Create subplot
    subplot(2,1,m)
    hold on
    
    % Plotting params
    sessions_to_use = Monkeys(m).Sessions_to_use;
    xticklabs = Monkeys(m).XTickLabs;
    
    % To make sure array data plots in the same order across mks (unique()
    % should be sorted alphabetically so this is probably redundant)
    if regexp(Monkeys(m).Name, 'Marta\w*', 'once')
        all_arrays = {'anterior', 'middle', 'posterior'};
    elseif regexp(Monkeys(m).Name, 'Max\w*', 'once')
%         all_arrays = {'anterior', 'middle', 'posterior', 'teo'};
        all_arrays = {'anterior', 'middle', 'posterior'};
    end
    
    % Pre allocate
    nUnits = zeros(length(sessions_to_use), length(all_arrays));
    nVR_Units = zeros(length(sessions_to_use), length(all_arrays));
    
    % Get data to plot
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        arrays = unique({Monkeys(m).Sessions(sessn).UnitInfo.Location});
        VR_bools = [Monkeys(m).Sessions(sessn).UnitInfo.(vrID)] < vr_alpha;
        for j = 1:length(arrays)
            idx = find(strcmp(all_arrays, arrays{j}));
            array_bools = strcmp(arrays{j}, {Monkeys(m).Sessions(sessn).UnitInfo.Location});
            nUnits(i,idx) = sum(array_bools);
            nVR_Units(i,idx) = sum(array_bools & VR_bools);
        end
    end
    
    % Tell MATLAB to repeat only the first three default colors
    cols = colororder;
    colororder(cols(1:3,:))
    
    % Plot all neurons with empty bars 
    b = bar(nUnits, 'grouped');
    for iB = 1:length(b)
        b(iB).EdgeColor = b(iB).FaceColor;
        b(iB).FaceColor = [1 1 1];
        b(iB).LineWidth = 2;
    end
    
    % Plot filled-in bars with num VR units
    b = bar(nVR_Units, 'grouped');
    for iB = 1:length(b)
        b(iB).EdgeColor = b(iB).FaceColor;
    end
    
    % Format the plot
    xticks(1:length(nUnits))
%     xticklabels(get_xtick_labs_colored(xticklabs, short_names, xtickcolors))
    xticklabels(xticklabs)
%     legend(all_arrays)
    ylim([0 150])
    yticks(0:75:150)
    ylabel('# units')
    set(gca, 'FontSize', 36)
end

%% Plot proportion VR units per day, across arrays


test_intervals = {[75 175], [175 350]};
baseline_intervals = {[-150 -50], [-175 0]};

% by all imgs

for m = 1:length(Monkeys)
    num_units = [];
    legend_strs = {};
    for i = 1:length(Monkeys(m).Sessions)
        for k = 1:length(test_intervals)
            test_int = test_intervals{k};
            baseline_int= baseline_intervals{k};
            tid = get_good_interval_name2(test_int, 'full', 'VisResp_all_test');
            bid = get_good_interval_name2(baseline_int, '', '');
            vrID = strcat(tid,bid);
            num_units(i,k) = sum([Monkeys(m).Sessions(i).UnitInfo.(vrID)]<0.001) ./ length([Monkeys(m).Sessions(i).UnitInfo.(vrID)]);
        end
    end
    
    for k = 1:length(test_intervals)
        legend_strs{k} = sprintf('%d to %d', test_intervals{k}(1),test_intervals{k}(2));
    end
    figure
    hold on
    bar(num_units, 'grouped')
    legend(legend_strs);
    xlabel('Session')
    xticks(1:length(Monkeys(m).Sessions))
    ylabel('Proportion units VR')
    ylim([0.25 0.75])
    title(sprintf('%s, VR units (all stims)', Monkeys(m).Name), 'Interpreter', 'none')
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',....
        'XGrid', 'on','YGrid', 'on', 'fontsize',20, ...
        'fontname', 'Helvetica','FontWeight','Bold', ...
        'XColor', 'black', 'YColor', 'black')
    set(gcf,'color','w');
end

%% Plot proportion VR units per day, separate arrays


test_intervals = {[75 175], [175 350]};
baseline_intervals = {[-150 -50], [-175 0]};
normalize = true;
alpha = 0.001;

for m = 1:length(Monkeys)
    figure
    hold on
    
    if regexp(Monkeys(m).Name, 'Marta\w*', 'once')
        all_arrays = {'anterior', 'middle', 'posterior'};
    elseif regexp(Monkeys(m).Name, 'Max\w*', 'once')
        all_arrays = {'anterior', 'middle', 'posterior', 'teo'};
    end
    
    for k = 1:length(test_intervals)
        test_int = test_intervals{k};
        baseline_int= baseline_intervals{k};
        tid = get_good_interval_name2(test_int, 'full', 'VisResp_all_test');
        bid = get_good_interval_name2(baseline_int, '', '');
        vrID = strcat(tid,bid);
        
        num_units = [];        
        for i = 1:length(Monkeys(m).Sessions)
            arrays = unique({Monkeys(m).Sessions(i).UnitInfo.Location});
            contrast_bools = [Monkeys(m).Sessions(i).UnitInfo.(vrID)] < alpha;
            for a = 1:length(arrays)
                idx = find(strcmp(all_arrays, arrays{a}));
                array_bools = strcmp(all_arrays{idx}, {Monkeys(m).Sessions(i).UnitInfo.Location});
                if normalize
                    norm_factor = sum(array_bools);
                else
                    norm_factor = 1;
                end
                num_units(i,idx) = sum(contrast_bools & array_bools) ./ norm_factor;
            end
        end
        
        subplot(1, length(test_intervals), k)
        bar(num_units, 'grouped')
        
        xlabel('Session')
        xticks(1:length(Monkeys(m).Sessions))
        xticklabels({Monkeys(m).Sessions.ShortName})
        xtickangle(90)
        if normalize
            ylim([0 1])
            ylabel('Proportion units VR')
        else
           ylim([0 150]) 
           ylabel('Number units VR')
        end
        title(sprintf('Interval %d to %d',  test_int(1), test_int(2)), 'Interpreter', 'none')
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], ...
            'XMinorTick', 'off', 'YMinorTick', 'off',....
            'XGrid', 'on','YGrid', 'on', 'fontsize',20, ...
            'fontname', 'Helvetica','FontWeight','Bold', ...
            'XColor', 'black', 'YColor', 'black')
        set(gcf,'color','w'); 
    end
    legend(all_arrays, 'Position', [0.2 0.9 0.5 0.075], 'Orientation', 'horizontal');
    sgtitle(sprintf('%s, VR units (all stims)', Monkeys(m).Name))
end

