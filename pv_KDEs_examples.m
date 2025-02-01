%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';

% Pick KDE bandwidth parameter and load
% bw = 20;
% fname = sprintf('MaxMarta_fix_cat_xma2_bootstrappedKDEs_BW%d.mat', bw);
bw = 15;
% fname = sprintf('MaxMarta_fix_cat_xma2_bootstrappedKDEs_ASYMM_BW%d_Marta_fix_cat_xma2', bandwidth);
fname = sprintf('MaxMarta_fix_cat_xma2_bootstrappedKDEs_ASYMM_BW%d_Max_fix_cat_xma2', bw);
Data = load(fullfile(EXT_HD, pv_path, fname));
[status, KDE] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Tmp hot fix
% just working with marta's data rn, delete Max
% KDE = KDE(1);

% Max data only
KDE = KDE(2);
rSessionsByMonk = {[6 7]};


%% Add Area labels, get short session names
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

%% Get xvals where true diffs are larger than bootstrapped popn (dynamic boundary) 

pct = 99;
diffid = sprintf('TrueDiff_BW%d', bw);
bootid = sprintf('BootstrappedDiffs_BW%d', bw);
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

%% Print out all neurons with large cat/dog differences in this analysis

% rSessionsByMonk = {[7 9], [6 7]};

pct = 99;
diffid = sprintf('TrueDiff_BW%d', bw);
bootid = sprintf('BootstrappedDiffs_BW%d', bw);

for m = 1:length(KDE)
    
    KDE(m).Name = KDE(m).Name;
    sessions_to_plot = rSessionsByMonk{m};
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        KDE(m).Sessions(sessn).Potential_example_units = [];
        
        %         units_to_plot = 1:length(KDE(m).Sessions(sessn).UnitInfo);
        units_to_plot = find(ismember({KDE(m).Sessions(sessn).UnitInfo.Location}, {'anterior', 'middle', 'posterior'}));
%         units_to_plot = find(ismember({KDE(m).Sessions(sessn).UnitInfo.Location}, {'anterior'}));
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            
            if isfield(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues, 'KDEYVals') && ...
                isnan(KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals)
                continue
            end
            
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals;
%             true_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid);
            catg_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.CatgVals;
            num_total_cat_trials = KDE(m).Sessions(sessn).NumCatgTrials(1);
            num_total_dog_trials = KDE(m).Sessions(sessn).NumCatgTrials(2);
            catg1_scaled = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg1_raw * 1000 * sum(catg_vals==1) / num_total_cat_trials;
            catg2_scaled = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg2_raw * 1000 * sum(catg_vals==2) / num_total_dog_trials;
            true_diffs = catg1_scaled - catg2_scaled;
            
            boot_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(bootid);
            boot_diffs = abs(boot_diffs); % just look at magnitude, not sign
            
            % reshape to prepare for fixed boundary calculation
            boot_diffs = reshape(boot_diffs, 1, numel(boot_diffs));
            
            percentile_diff = prctile(boot_diffs, pct);
            exceedid = 'ExceedBD_FixedBound';
            times = find(abs(true_diffs) > percentile_diff); % normalized times in ms where true diffs exceeds bootstrapped diffs
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(exceedid) = times;
            if numel(times) > 15 && sum(times > 200)>15
                fprintf('Monkey %s, session %d, unit %d \n', KDE(m).Name(1:3), sessn, unit)
                KDE(m).Sessions(sessn).Potential_example_units = [KDE(m).Sessions(sessn).Potential_example_units unit];
            end
            
        end
    end
end

%% Show all potential example units

for m = 1:length(KDE)
    sessions_to_plot = rSessionsByMonk{m};
    
    for i = 1:length(sessions_to_plot)    
        figure2
        hold on 
        
        sessn = sessions_to_plot(i);
        units = KDE(m).Sessions(sessn).Potential_example_units;
        
        for j = 1:length(units)
            unit = units(j);
            
            subplot(8, 8,j)
            hold on
            
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals;
%             true_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid);
            catg_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.CatgVals;
            num_total_cat_trials = KDE(m).Sessions(sessn).NumCatgTrials(1);
            num_total_dog_trials = KDE(m).Sessions(sessn).NumCatgTrials(2);
            catg1_scaled = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg1_raw * 1000 * sum(catg_vals==1) / num_total_cat_trials;
            catg2_scaled = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg2_raw * 1000 * sum(catg_vals==2) / num_total_dog_trials;
            
            plot(kde_x_vals, catg1_scaled, 'r-')
            plot(kde_x_vals, catg2_scaled, 'b-')
            title(unit)
%             formatSVMPlot(gca, gcf)
        end
        sgtitle(KDE(m).Sessions(sessn).ShortName)
    end
end


%% (revision) Plot units for paper

% Marta
% rSessionsByMonk = {[9], [7]};  % post-training sessions
% units_by_monk = {[5, 64, 76, 114], []};

% Max
units_by_monk = {[285 13 206 164]};
sessions_to_plot = [7];

for m = 1:length(KDE)
    
    for i = 1:length(sessions_to_plot)
        
        figure2
        hold on 
        
        sessn = sessions_to_plot(i);
%         units = KDE(m).Sessions(sessn).Potential_example_units;
        units = units_by_monk{m};
        
        for j = 1:length(units)
            unit = units(j);
            
            subplot(1, length(units),j)
            hold on
            
            kde_x_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals;
%             true_diffs = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.(diffid);
            catg_vals = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.CatgVals;
            num_total_cat_trials = KDE(m).Sessions(sessn).NumCatgTrials(1);
            num_total_dog_trials = KDE(m).Sessions(sessn).NumCatgTrials(2);
            catg1_scaled = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg1_raw * 1000 * sum(catg_vals==1) / num_total_cat_trials;
            catg2_scaled = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg2_raw * 1000 * sum(catg_vals==2) / num_total_dog_trials;
            
            plot(kde_x_vals, catg1_scaled, 'r-')
            plot(kde_x_vals, catg2_scaled, 'b-')
%             title(unit)
            formatSVMPlot(gca, gcf)
            if j == 1
                ylabel("Spikes / sec")
            end
        end
    end
end


% %% (old) Plot units for paper
% 
% % Example units [monk, sessn, unit]
% 
% 
% % Sorted by array and rough spikr rate within array
% % unitIDs = {[1 7 24], [1,7,6], [2,7,170], [1,7,2]}; arr = 'posterior'; % posterior
% % unitIDs = {[2,6,187], [2 6 224], [1,9,105], [1 7 112]}; arr = 'middle'; % middle
% % unitIDs = {[2 7 10], [1 9 78], [1 9 76], [1 7 80]}; arr = 'anterior'; % anterior [m,sessn,unit]
% 
% % Sorted by array and rough time of cat/dog diff
% % (early, late, both)
% unitIDs_post = {[1,7,6], [2,7,170], [2,7,284]}; 
% unitIDs_mid = {[1 7 112], [2,6,187], [1,9,105]}; 
% unitIDs_ant = {[1 7 80], [1 9 96], [2 7 12]}; % prefer [1 9 78] for extended here, but want one of Max's
% pct = 99;
% 
% concat_unitIDs = {unitIDs_post unitIDs_mid unitIDs_ant};
% arrayNames = {'posterior', 'middle', 'anterior'};
% 
% for iArr = 1:length(arrayNames)
%     unitIDs = concat_unitIDs{iArr};
%     arr = arrayNames{iArr};
%     figure2('Position', [400 400 400 600])
%     for i = 1:length(unitIDs)
%         subplot(length(unitIDs),1,i)
%         hold on
%         plotCatDogKDEs(KDE, unitIDs{i}(1), unitIDs{i}(2), unitIDs{i}(3), bw, pct)
%         if i == length(unitIDs)
%             xticks(-200:200:400)
%             xlabel('Time from cue on (ms)')
%             ylabel('Spikes / sec')
%         else
%             xticks(0)
%         end
%     end
%     saveas(gcf, sprintf('../exampleKDEs_%s', arr), 'epsc')
% end


%% Funcions
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


