

%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';

% Load spike data
Data = load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_ni.mat'));
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

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

%% Collect spike counts for all neurons

% intervals_to_test = {[-150 -50], [-200 0], [-250 50], [0 75], [75 175], [175 275], [275 375], [75 275], [75 375]}; % time window from cue on
% intervals_to_test = {[75 175], [175 275], [275 375], [175 375]}; % time window from cue on
% intervals_to_test = {[-150 -50], [-175 0], [75 175], [175 350]};
intervals_to_test = {[-150 -50], [175 275]};
catgs = {1:260, 261:520};

for m = 1:length(Monkeys)
    sessions_to_use = 1:length(Monkeys(m).Sessions);
        
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        units_to_use = 1:length(Monkeys(m).Sessions(sessn).UnitInfo);
        
        for p = 1:length(intervals_to_test)
            interval_to_test = intervals_to_test{p};
            
            % Collect all times that an image came on
            catg1 = Monkeys(m).Sessions(sessn).CueInfo(catgs{1});
            catg2 = Monkeys(m).Sessions(sessn).CueInfo(catgs{2});
            catg1_timeson = vertcat(catg1.Times_on);
            catg2_timeson = vertcat(catg2.Times_on);

            
            % Pre-allocate spike times
            X = zeros(length(catg1_timeson)+length(catg2_timeson), length(units_to_use));

            for j = 1:length(units_to_use)

                % Get and sort this unit's spike times
                unit = units_to_use(j);
                allspike_times = sort(Monkeys(m).Sessions(sessn).UnitInfo(unit).Spike_times);

                % Collect spike counts for catg1
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

                 % Collect spike counts for catg2
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


                % Concat into X
                X(1:length(catg1_sc),j) = catg1_sc;
                X(length(catg1_sc)+1:end,j) = catg2_sc;


            end
            % Store
            testID = get_good_interval_name2(interval_to_test, 'full', 'X');
            Monkeys(m).Sessions(sessn).(testID) = X;
            
        end
    end
end

%% Get overall visual responsiveness

% test_intervals = {[75 175], [175 275], [275 375], [75 375]};
% baseline_intervals = {[-150 -50], [-150 -50], [-150 -50], [-250 50]};
% test_intervals = {[175 350], [75 175]};
% baseline_intervals = {[-175 0], [-150 -50]};
test_intervals = {[175 275]};
baseline_intervals = {[-150 -50]};

for m = 1:length(Monkeys)
    for i = 1:length(Monkeys(m).Sessions)
        for k = 1:length(test_intervals)
            
            % Find the requested spike count fields
            test_int = test_intervals{k};
            baseline_int= baseline_intervals{k};
            testID = get_good_interval_name2(test_int, 'full', 'X');
            baseID = get_good_interval_name2(baseline_int, 'full', 'X');
            testSC = Monkeys(m).Sessions(i).(testID); % test int
            baselineSC = Monkeys(m).Sessions(i).(baseID); % baseline int

            % Field names for results
            tid = get_good_interval_name2(test_int, 'full', 'VisResp_all_TTEST');
            bid = get_good_interval_name2(baseline_int, '', '');
            fullID = strcat(tid,bid);

            % Test each unit
            for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
                [~, Monkeys(m).Sessions(i).UnitInfo(j).(fullID)] = ttest2(baselineSC(:,j), testSC(:,j), 'tail', 'left'); 
                
                % Old way with non-parametric testing
%                 Monkeys(m).Sessions(i).UnitInfo(j).(contrast_id) = signrank(X2(:,j) - X1(:,j), 0, 'tail', 'left'); % tests that median(baseline-test) < 0
            end
        end
    end
end

%% Get visual responsiveness to each image

% test_intervals = {[75 175], [175 275], [275 375], [75 375]};
% baseline_intervals = {[-150 -50], [-150 -50], [-150 -50], [-250 50]};
% test_intervals = {[175 350], [75 175]};
% baseline_intervals = {[-175 0], [-150 -50]};
test_intervals = {[175 275]};
baseline_intervals = {[-150 -50]};

groupings = 1:520; % vis resp for each image (count of this gives rough image-level selectivity)


for m = 1:length(Monkeys)
    for i = 1:length(Monkeys(m).Sessions)
        for k = 1:length(test_intervals)
            
            % Find the requested spike count fields
            test_int = test_intervals{k};
            baseline_int= baseline_intervals{k};
            testID = get_good_interval_name2(test_int, 'full', 'X');
            baseID = get_good_interval_name2(baseline_int, 'full', 'X');
            testSC = Monkeys(m).Sessions(i).(testID); % test int
            baselineSC = Monkeys(m).Sessions(i).(baseID); % baseline int
            
            % Pre-allocate, (num imgs) x (num units)
            img_vr_mat = zeros(length(groupings), length(Monkeys(m).Sessions(i).UnitInfo));
            
            % For each image...
            for p = 1:length(groupings)
                
                % Get field name for results
%                 tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img');
                tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img_TTEST');
                bid = get_good_interval_name2(baseline_int, '', '');
                fullID = strcat(tid,bid);   
                
                % Find relevant trials (trials when that image was shown)
                idx = find(Monkeys(m).Sessions(i).Session_Y_imageID == p);
                
                % Test each unit on those trials
                for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
                    [~,img_vr_mat(p,j)] = ttest2(baselineSC(idx,j), testSC(idx,j), 'tail', 'left'); 
                    
                    % Old way of non-parametric testing
%                     img_vr_mat(p,j) = signrank(baselineSC(idx,j), testSC(idx,j), 'tail', 'left'); % tests that median(baseline-test) < 0
%                     signrank(reshape(baselineSC(idx,j),1,numel(baselineSC(idx,j))) - reshape(testSC(idx,j),1,numel(testSC(idx,j))),...
%                         0, 'tail', 'left')
                end
            end
            
            % Store results
            Monkeys(m).Sessions(i).(fullID) = img_vr_mat;
        end
    end
end

%% Save data
fname = 'MaxMarta_VR_img_TTest_v2.mat'; % v2 has 175-275
t = whos('Monkeys');
if t.bytes > 2e9 % split large struct and save all the vars
    [status, vnames] = split_monkeyStruct_in_parts(Monkeys);
    save(fullfile(EXT_HD, pv_path, fname), vnames{:})
else
    save(fullfile(EXT_HD, pv_path, fname), 'Monkeys')
end




