% check whether the neurons that respond to more images post-training also
% have reduced spike rates, to resolve the paradoxical unchanged GLM coeff
% result.

% hypothesis is that if tuning curves broadened, maybe the overall rate
% also went down (ie literally taking a probability density and broadening
% it), which could explain why we see units "visually responsive" to
% greater numbers of stimuli but no change in the GLM coefficients post
% training.


%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';
figureSavePath = '/Users/jonahpearl/Documents/BJR group/Catdog_paper/Feb_2023_addtl_figs';
% fname = 'MaxMarta_VR_img_TTest_v2.mat';
fname = 'MaxMarta_VR_img_TTest_jun2021.mat'; % changed ttest2 to ttest (because they're paired!)
Data = load(fullfile(EXT_HD, pv_path, fname));
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Set parameters

% GLM Parameters
rSessionsByMonk = {[7 9], [6 7]};
% rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};


% Other Parameters
spikeCountPath = 'XMA2/Spike_count_mats';
TE_LOCS = {'anterior', 'middle', 'posterior'};
catg2_ind1 = 261;
loc_in_trials_fname_base = '%s_imgLocInTrials.mat';
trial_num_fname_base = '%s_trialNum.mat';

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

%% Look at spike rates to images to which units were visually responsive

test_intervals = {[175 275]};
baseline_intervals = {[-150 -50]};
area = 'te';
vr_alpha = 0.05;

for m = 1:length(Monkeys)
    sessions_to_use = rSessionsByMonk{m};
    
    for p = 1:length(test_intervals)
        test_int = test_intervals{p};
        baseline_int = baseline_intervals{p};
%         tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img'); % using non-parametric testing
        tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img_TTEST'); % with a t-test
        xid = get_good_interval_name2(test_int, 'full', 'X');
        vr_spike_rate_id = get_good_interval_name2(test_int, '', 'VRSpikeInfo');
        bid = get_good_interval_name2(baseline_int, '', '');
        vr_id = strcat(tid,bid);
        
        for i = 1:length(sessions_to_use)
            sessn = sessions_to_use(i);

            % Get data
            data_mat = Monkeys(m).Sessions(sessn).(vr_id); % stimuli x units
            if strcmp(area, 'te')
                area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, area);
            else
                area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Location}, area);
            end
            data_mat = data_mat(:, area_bool); % imgs x units
            [responsive_paired_imgs, responsive_paired_units] = find(data_mat < vr_alpha);
            
            
            % for each img in list, grab corresponding trials
            ii = 1;  % counter for struct we're going to store data in
            uq_imgs = unique(responsive_paired_imgs);
            for iImg = 1:length(uq_imgs)
                this_img = uq_imgs(iImg);
                if this_img >= catg2_ind1
                    this_catg = 2;
                else
                    this_catg = 1;
                end
                these_trials = Monkeys(m).Sessions(sessn).CueInfo(this_img).TrialNum;
                uq_units = unique(responsive_paired_units(responsive_paired_imgs == this_img));
                for iUnit = 1:length(uq_units)
                    % for each unit, grab the spike rates from those trials
                    this_unit = uq_units(iUnit);
                    scs = Monkeys(m).Sessions(sessn).(xid)(these_trials, this_unit);
                    Monkeys(m).Sessions(sessn).(vr_spike_rate_id)(ii) = struct('Unit', this_unit,...
                        'Image', this_img, 'Category', this_catg, ...
                        'Trials', these_trials, 'SpikeCounts', scs);
                    ii = ii + 1;
                end
            end
        end
    end
end

%% Plot results
for m = 1:length(Monkeys)
    sessions_to_use = rSessionsByMonk{m};
    
    for p = 1:length(test_intervals)
        test_int = test_intervals{p};
        vr_spike_rate_id = get_good_interval_name2(test_int, '', 'VRSpikeInfo');
        
        catg_spike_counts = {};
        vr_counts = {};
        
        for i = 1:length(sessions_to_use)
            sessn = sessions_to_use(i);
            
            % concat all spike counts by session and category
            for catg = [1 2]
                idx = find([Monkeys(m).Sessions(sessn).(vr_spike_rate_id).Category] == catg);
                scs = vertcat(Monkeys(m).Sessions(sessn).(vr_spike_rate_id)(idx).SpikeCounts);
                catg_spike_counts{i, catg} = scs;
            end
            
            % show how many total images each neuron in this set responded to
            uq_units = unique([Monkeys(m).Sessions(sessn).(vr_spike_rate_id).Unit]);
            vr_count = zeros(length(uq_units),1);
            for iUnit = 1:length(uq_units)
                vr_count(iUnit) = sum([Monkeys(m).Sessions(sessn).(vr_spike_rate_id).Unit] == uq_units(iUnit));
            end
            vr_counts{i} = vr_count;
        end
        
        % just concat dogs and cats to start
        figure
        hold on
        all_pre = vertcat(catg_spike_counts{1, 1}, catg_spike_counts{1, 2});
        all_post = vertcat(catg_spike_counts{2, 1}, catg_spike_counts{2, 2});
        histogram(all_pre, 'DisplayName', 'Pre', 'Normalization', 'probability')
        histogram(all_post, 'DisplayName', 'Post', 'Normalization', 'probability')
        title(Monkeys(m).Name, 'Interpreter', 'none')
        xlabel('Spike counts (across trials) for responsive unit/img pairs')
        legend
        saveas(gcf, fullfile(figureSavePath, sprintf('%s_spikeCounts_vr_pairs.pdf', Monkeys(m).Name)))
        
        % try separating them --> hists look ~identical.
%         catg_names = {'Cats', 'Dogs'};
%         for catg = [1 2]
%             figure
%             hold on
%             histogram(catg_spike_counts{1, catg}, 'DisplayName', 'Pre', 'Normalization', 'probability')
%             histogram(catg_spike_counts{2, catg}, 'DisplayName', 'Post', 'Normalization', 'probability')
%             title(sprintf('%s, %s', Monkeys(m).Name, catg_names{catg}), 'Interpreter', 'none')
%         end

        figure
        hold on
        histogram(vr_counts{1}, 1:15:450, 'DisplayName', 'Pre', 'Normalization', 'probability')
        histogram(vr_counts{2}, 1:15:450, 'DisplayName', 'Post', 'Normalization', 'probability')
        title(Monkeys(m).Name, 'Interpreter', 'none')
        xlabel('Num imgs to which units responded')
        legend
        saveas(gcf, fullfile(figureSavePath, sprintf('%s_numImgsVR.pdf', Monkeys(m).Name)))
    end
end