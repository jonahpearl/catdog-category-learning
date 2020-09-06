%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';

% Pick KDE bandwidth parameter and load
bw = 20;
fname = sprintf('MaxMarta_fix_cat_xma2_bootstrappedKDEs_BW%d.mat', bw);
Data = load(fullfile(EXT_HD, pv_path, fname));
[status, KDE] = stitch_monkeyStruct_from_parts(Data);
clear Data

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

%% Print out all neurons with signf cat/dog differences in this analysis

pct = 99;
diffid = sprintf('TrueDiff_BW%d', bw);
bootid = sprintf('BootstrappedDiffs_BW%d', bw);

% for m = 1:length(KDE)
for m = 2
   if strcmp(KDE(m).Name, 'Marta_fix_cat_xma2')
%         sessions_to_plot = [1 2 7 9]; % skip base04 (outlier) and post01 (low trial count)
        sessions_to_plot = [1:3 5:9];
    elseif strcmp(KDE(m).Name, 'Max_fix_cat_xma2')
%         sessions_to_plot = [1 2 6 7]; % This is base 01, base02, pre and post. skip base03,04,05 and sub03 (low trial count)
        sessions_to_plot = 1:7;
   end
    
    KDE(m).Name = KDE(m).Name;
    sessions_to_plot = 7;
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        %         units_to_plot = 1:length(KDE(m).Sessions(sessn).UnitInfo);
%         units_to_plot = find(ismember({KDE(m).Sessions(sessn).UnitInfo.Location}, {'anterior', 'middle', 'posterior'}));
        units_to_plot = find(ismember({KDE(m).Sessions(sessn).UnitInfo.Location}, {'anterior'}));
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
                fprintf('Monkey %s, session %d, unit %d \n', KDE(m).Name(1:3), sessn, unit)
            end
            
        end
    end
end

%% Look through those units
m = 2;
sessn = 7;
units = [5 10 12 13 14 15 20]; 
% Martha sessn 7: 2, 6, 24 (later), 58, 90 (later), 112
% Martha sessn 9: 105, 114.
% Martha sessn 9 anterior late: 81, 85, 96
% Max sessn 6: 187, 300, 317, 
% Max sessn 7: 10, 170
% Max sessn 7 posterior both: 284, 285
% % Martha, sessn 1, 168 10 120 102 

bw = 20; 
pct = 99;
figure2
for j = 1:length(units)
    unit = units(j);
    subplot(length(units),1,j)
    hold on
    
    plotCatDogKDEs(KDE, m, sessn, unit, bw, pct)
end

%% Plot units for paper

% Example units [monk, sessn, unit]


% Sorted by array and rough spikr rate within array
% unitIDs = {[1 7 24], [1,7,6], [2,7,170], [1,7,2]}; arr = 'posterior'; % posterior
% unitIDs = {[2,6,187], [2 6 224], [1,9,105], [1 7 112]}; arr = 'middle'; % middle
% unitIDs = {[2 7 10], [1 9 78], [1 9 76], [1 7 80]}; arr = 'anterior'; % anterior [m,sessn,unit]

% Sorted by array and rough time of cat/dog diff
% (early, late, both)
unitIDs_post = {[1,7,6], [2,7,170], [2,7,284]}; 
unitIDs_mid = {[1 7 112], [2,6,187], [1,9,105]}; 
unitIDs_ant = {[1 7 80], [1 9 96], [2 7 12]}; % prefer [1 9 78] for extended here, but want one of Max's

concat_unitIDs = {unitIDs_post unitIDs_mid unitIDs_ant};
arrayNames = {'posterior', 'middle', 'anterior'};

for iArr = 1:length(arrayNames)
    unitIDs = concat_unitIDs{iArr};
    arr = arrayNames{iArr};
    figure2('Position', [400 400 400 600])
    for i = 1:length(unitIDs)
        subplot(length(unitIDs),1,i)
        hold on
        plotCatDogKDEs(KDE, unitIDs{i}(1), unitIDs{i}(2), unitIDs{i}(3), bw, pct)
        if i == length(unitIDs)
            xticks(-200:200:400)
            xlabel('Time from cue on (ms)')
            ylabel('Spikes / sec')
        else
            xticks(0)
        end
    end
    saveas(gcf, sprintf('../exampleKDEs_%s', arr), 'epsc')
end

%% Functions
function plotCatDogKDEs(KDE, m, sessn, unit, bw, pct)
    UnitKDE = KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues; 
    srs = UnitKDE.KDEYVals_SRs;
    f1 = 1000 * srs(2) / KDE(m).Sessions(sessn).NumCatgTrials(1) * UnitKDE.KDEYVals_catg1_raw;
    f2 = 1000 * srs(3) / KDE(m).Sessions(sessn).NumCatgTrials(2) * UnitKDE.KDEYVals_catg2_raw;

    diffid = sprintf('TrueDiff_BW%d', bw);
    bootid = sprintf('BootstrappedDiffs_BW%d', bw);
    true_diffs = UnitKDE.(diffid);
    boot_diffs = UnitKDE.(bootid);
    kde_x_vals = UnitKDE.KDEXVals;
    
    boot_diffs_reshape = reshape(boot_diffs, 1, numel(boot_diffs));
    percentile_diff = prctile(abs(boot_diffs_reshape), pct); % 1 x (num xvals) vector of the 95th percentile of bootstrapped data
    t = find(abs(true_diffs) > percentile_diff);

    plot(kde_x_vals, f1, 'r-', 'LineWidth', 2)
    plot(kde_x_vals, f2, 'b-', 'LineWidth', 2)
    
    yl = ylim;
    plot(kde_x_vals(t), yl(2)*1.06, 'ko','MarkerFaceColor', 'k', 'MarkerSize', 4)
    ylim([0 yl(2)*1.12])
    % xlabel('Time from cue on (ms)')
    % ylabel('Spikes / second / trial')
    xlim([-200 500])
    shortName = regexp(KDE(m).Sessions(sessn).ShortName, '(Pre|Post)', 'match', 'once');
%     title(sprintf('%s, %s, # %d (%s)', KDE(m).Name(1:3), KDE(m).Sessions(sessn).ShortName, unit, KDE(m).Sessions(sessn).UnitInfo(unit).Location), 'Interpreter', 'none')
%     title(sprintf('%s, %s, %s', upper(KDE(m).Name(3)), shortName, KDE(m).Sessions(sessn).UnitInfo(unit).Location))
%     title(sprintf('%s', KDE(m).Sessions(sessn).UnitInfo(unit).Location))
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',18, ...
        'fontname', 'Helvetica', 'fontweight', 'normal', ...
        'XColor', 'black', 'YColor', 'black')
    h = gca;
    h.Title.FontWeight = 'normal';
end