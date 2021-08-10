% analyze visual responsiveness of units to the cat/dog images

%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';
figureSavePath = '/Users/jonahpearl/Documents/BJR group/Catdog paper/';
% fname = 'MaxMarta_VR_img_TTest_v2.mat';
fname = 'MaxMarta_VR_img_TTest_jun2021.mat'; % changed ttest2 to ttest (because they're paired!)
Data = load(fullfile(EXT_HD, pv_path, fname));
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

%% Run major-axis regression to assess slopes, store residuals
% MAR alpha is not dynamic for figures -- if you want to use a different
% MAR alpha (for slope CI's), come back here and re-run the analysis first.

vr_alpha = 0.05;
mar_alpha = 0.05;
alpha_string_scale_factor = 100;
% test_intervals = {[75 175], [175 350]};
% baseline_intervals = {[-150 -50], [-175 0]};
% test_intervals = {[175 350]};
% baseline_intervals = {[-175 0]};
test_intervals = {[175 275]};
baseline_intervals = {[-150 -50]};
make_plots = false;
areas_to_plot = {'te', 'anterior', 'middle', 'posterior'};

for m = 1:length(Monkeys)
    sessions_to_use = 1:length(Monkeys(m).Sessions);
    
    for p = 1:length(test_intervals)
        test_int = test_intervals{p};
        baseline_int = baseline_intervals{p};
%         tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img'); % using non-parametric testing
        tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img_TTEST'); % with a t-test
        bid = get_good_interval_name2(baseline_int, '', '');
        vr_id = strcat(tid,bid);
        
        for a = 1:length(areas_to_plot)
            area = areas_to_plot{a};
            if make_plots
                figure('Position', [400 400 1300 1000])
                hold on
            end
            for i = 1:length(sessions_to_use)
                sessn = sessions_to_use(i);
                
                % Get data
                data_mat = Monkeys(m).Sessions(sessn).(vr_id); % stimuli x units
                if strcmp(area, 'te')
                    area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, area);
                else
                    area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Location}, area);
                end
                data_mat = data_mat(:, area_bool);
                cat_props = sum(data_mat(1:260,:) < vr_alpha) / 260;
                dog_props = sum(data_mat(261:520,:) < vr_alpha) / 260;
                
                % Calculate correlation
                x = corrcoef(cat_props, dog_props);
                
                % Get angle of best-fit line with major axis regression
                % method. I checked this function against a reference, it's
                % correct.
                [coeffs,coeff_ints,~,~] = maregress(dog_props, cat_props, mar_alpha);
                theta = rad2deg(atan(coeffs(2)));
                theta_bounds = rad2deg(atan(coeff_ints(2,:)));
                intercept = coeffs(1);
                intercept_bounds = coeff_ints(1,:);
                
                % Print for easy reference
                fprintf('%s, %s, %s, interval %d to %d, correlation %0.2f, slope %0.2f \n', ...
                    Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName,...
                    areas_to_plot{a}, test_int(1), test_int(2),...
                    x(1,2), tan(theta/180*pi))
                
                % Calculate distance from points to regression line:
                % D = abs(ax0+by0+c) / sqrt(a^2+b^2), where point is
                % (x0,y0), and line is ax+by+c = 0.
                % We have y=mx+c form, so a=coeffs(2), b = -1, and c =
                % coeffs(1).
                perpendicular_residuals = abs(coeffs(2)*dog_props + -1*cat_props + coeffs(1)) / sqrt(coeffs(2)^2 + 1);
                
                % Store data
                propn_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_CatDogPropns_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                slope_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Slope_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                intercept_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Intercept_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                resid_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Resids_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                rval_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_RValue_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                Monkeys(m).Sessions(sessn).(propn_id) = [cat_props; dog_props]; % 2 x units
                Monkeys(m).Sessions(sessn).(resid_id) = perpendicular_residuals; % 1 x units
                Monkeys(m).Sessions(sessn).(slope_id) = [theta_bounds(1) theta theta_bounds(2)]; % 1 x 3
                Monkeys(m).Sessions(sessn).(intercept_id) = [intercept_bounds(1) intercept intercept_bounds(2)]; % 1 x 3
                Monkeys(m).Sessions(sessn).(rval_id) = x(1,2); % pearson correlation coeff.
                
                % Plot if requested
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
            end
            if make_plots
                sgtitle(sprintf('%s, area %s, %d to %d, image VR propns,(alpha = %0.2f)', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), vr_alpha), 'Interpreter', 'none')
            end
        end
    end
end

%% Plot selected MAR results with scatterplots

% interval_to_plot = [175 350];
interval_to_plot = [175 275];
vr_alpha = 0.05;
alpha_string_scale_factor = 100;
area_to_plot = 'te';
color_by_waveform_class = false;
colors = cbrewer('qual', 'Set1', 5);
f1 = figure('Position', [400 400 860 800]); % scatter plots
hold on
f2 = figure('Position', [400 400 400 800]); % histograms
f3 = figure('Position', [400 400 400 800]); % sums

hold on

for m = 1:length(Monkeys)
    sessions_to_use = Monkeys(m).Sessions_to_use;
    diffs = cell(1,length(sessions_to_use));
    sums = cell(1,length(sessions_to_use));
    resids = cell(1,length(sessions_to_use));
    
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        
        % Get data
        propn_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_CatDogPropns_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        t = Monkeys(m).Sessions(sessn).(propn_id);
        
        % Raw num of cat and dogs to which each neuron responds
        cat_props = t(1,:);
        dog_props = t(2,:);
        
        % Difference is like a selectivity measure
        diffs{i} = dog_props - cat_props;
        
        % Sum, to see if neurons respond to more images overall
        sums{i} = (dog_props + cat_props)*520;
        
        % Slope of fit line indicates population bias
        slope_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_Slope_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        slope = tan(deg2rad(Monkeys(m).Sessions(sessn).(slope_id)(2)));
        
        % Intercept is ~meaningless, mostly 0
        intercept_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_Intercept_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        intercept = Monkeys(m).Sessions(sessn).(intercept_id)(2);
        xvals = 0:0.01:1;
        
        % Residuals is similar to difference, since the fit lines are so
        % close to 45 degrees. Reflects variance of neural responsiveness
        % around its mean.
        resid_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Resids_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        perpendicular_resids = Monkeys(m).Sessions(sessn).(resid_id);
        resids{i} = perpendicular_resids;
        
        % ID session name and set some params accordingly.
        switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
            case 'Base'
                col = epoch_colors(1,:);
                thk = 1;
            case 'Pre'
                col = epoch_colors(2,:);
                thk = 2;
                pre_diff_idx = i;
            case 'Post'
                col = epoch_colors(3,:);
                thk = 2;
                post_diff_idx = i;
        end
                
        % Make scatter plot
        set(0, 'CurrentFigure', f1)
        subplot(length(Monkeys),length(sessions_to_use),...
            length(sessions_to_use)*(m-1) + mod(i-1, length(sessions_to_use)) + 1)
        hold on
        scatter(dog_props, cat_props, 'bo')
%         binscatter(dog_props, cat_props, 20)
%         scatterDiagHist(dog_props, cat_props, -0.2:0.025:0.2, 'bo')
        if i == 1 && m == 1
            xlabel('Fraction dogs causing exc.')
            ylabel('Fraction cats causing exc.')
        end
        xlim([0 1])
        ylim([0 1])
        xticks([0 0.5 1])
        yticks([0 0.5 1])
%         title(sprintf('%s', Monkeys(m).Sessions(sessn).ShortName), 'FontWeight', 'normal')
        title(Monkeys(m).XTickLabs{i},'FontWeight', 'normal')
        formatPlot(gca, f1)
        axis square
        
        % Color in example points
%         if i == 1 && m == 1
%             scatter(dog_props(14), cat_props(14), 'go', 'filled')
%             scatter(dog_props(2), cat_props(2), 'go', 'filled')
%         end

        % Plot reference lines
        set(0, 'CurrentFigure', f1)
        plot(xvals, xvals, 'k--', 'LineWidth', 2)
        plot(xvals, intercept + xvals*slope, 'r--', 'LineWidth', 2)
        
        
        
        % Make histograms to add to diag of scatters
        set(0, 'CurrentFigure', f2)
%         subplot(length(Monkeys),length(sessions_to_use),...
%             length(sessions_to_use)*(m-1) + mod(i-1, length(sessions_to_use)) + 1)
        subplot(length(Monkeys), 1, m)
        hold on
        histogram(abs(dog_props - cat_props), 'BinEdges', 0:0.025:0.3)
        xlabel('Selectivity |cat-dog|')
%         histogram(perpendicular_resids, 'BinEdges', 0:0.025:0.2)
        ax = gca;
        formatPlot(ax, f2)
        ax.YAxis.Visible = 'off';
        ax.YGrid = 'off';
        ax.XGrid = 'off';
%         ax.XLabel.FontSize = 40;
        ax.YAxis.Scale = 'log';
        

    end
    
%     % Informative sgtitle if desired
%     sgtitle(sprintf('Area %s, %d to %d, image VR propns,(alpha = %0.2f)',...
%         area_to_plot, interval_to_plot(1), interval_to_plot(2), vr_alpha), 'Interpreter', 'none')

    % Run variance tests on diffs
    pre_diffs = diffs{pre_diff_idx};
    post_diffs = diffs{post_diff_idx};
    
%     [h,pval] = vartest2(pre_diffs, post_diffs);
%     fprintf('DIFFS Pre var: %0.3g, post var: %0.3g. H = %d, p = %d \n', ...
%             var(pre_diffs), var(post_diffs), h, pval)
    [h,pval] = ttest2(abs(pre_diffs), abs(post_diffs));
    fprintf('ABS DIFFS pre vs post ttest2: h = %d, pval = %0.2f \n', ...
            h, pval)
    fprintf('Mean+std ABS DIFF pre: %0.2g +/- %0.2g, post: %0.2g + %02.g \n', ...
        mean(abs(pre_diffs)), std(abs(pre_diffs)), mean(abs(post_diffs)), std(abs(post_diffs)))
        
    % Test sums (did units responsd to more images overall?)
%     [h,pval] = ttest2(sums{pre_diff_idx}, sums{post_diff_idx});
%     fprintf('SUMS, pre: %0.0f, post: %0.0f, pre vs post ttest2: h = %d, pval = %0.2f \n', ...
%             mean(sums{pre_diff_idx}), mean(sums{post_diff_idx}), h, pval)
    
    
    % Ditto on residuals
%     [h,pval] = vartest2(resids{pre_diff_idx}, resids{post_diff_idx});
%     fprintf('RESIDUALS Pre var: %0.3g, post var: %0.3g. H = %d, p = %d \n', ...
%             var(resids{pre_diff_idx}), var(resids{post_diff_idx}), h, pval)

    % plot sums
%     set(0, 'CurrentFigure', f3)
%     subplot(length(Monkeys), 1, m)
%     hold on
%     errorbar(1:length(sessions_to_use), cellfun(@mean, sums), cellfun(@std, sums)./cellfun(@length,sums))
%     formatPlot(ax, f3)
%     title('Sums of cat and dog responsiveness')
end

% saveas(f1, fullfile(figureSavePath, sprintf('VR_scatter_%s', propn_id)), 'epsc')
% saveas(f2, fullfile(figureSavePath, sprintf('VR_hist_%s', propn_id)), 'epsc')

%% Draw each histogram separately

for m = 1:length(Monkeys)
    sessions_to_use = Monkeys(m).Sessions_to_use;
    diffs = cell(1,length(sessions_to_use));
    
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        
        % Get data
        figure
        propn_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_CatDogPropns_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        slope_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_Slope_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        intercept_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_Intercept_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        t = Monkeys(m).Sessions(sessn).(propn_id);
        cat_props = t(1,:);
        dog_props = t(2,:);
        histogram(dog_props - cat_props, 'BinEdges', -0.2:0.025:0.2)
        ax = gca;
        formatPlot(ax, gcf)
        ax.YAxis.Visible = 'off';
        ax.YGrid = 'off';
        ax.XGrid = 'off';
        saveas(gcf, fullfile(figureSavePath, 'VR_hists', sprintf('VR_hist_m%d_sessn%d_%s', m, sessn, propn_id)), 'epsc')
    end
end

%% Plot CDFs of MAR residuals, and do stats

vr_alpha = 0.05;
alpha_string_scale_factor = 100;
% ranksum_alpha = 0.025;
ranksum_alpha = 0.05;
test_intervals = {[175 275]};

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
                
                % ***********************************
                % use selectivity index (cats - dogs)
                % ***********************************
                propn_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_CatDogPropns_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
                t = Monkeys(m).Sessions(sessn).(propn_id);
                cat_props = t(1,:);
                dog_props = t(2,:);
                diffs = abs(cat_props - dog_props);
                [f, x] = ecdf(diffs);
                xlabel_str = 'Selectivity (|prop_cat - prop_dog|)';
                
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
                        pre_data = diffs;
                    case 'Post'
%                         col = 'g';
                        col = epoch_colors(3,:);
                        thk = 2;
                        post_data = diffs;
                    case 'Sub'
                        col = [0.8 0.8 0.8];
                        thk = 1;
                 end

                % ***********************************
                % use residuals from major axis regression
                % ***********************************
%                 propn_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_CatDogPropns_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
%                 resid_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_Resids_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
%                 t = Monkeys(m).Sessions(sessn).(propn_id);
%                 cat_props = t(1,:);
%                 dog_props = t(2,:);
%                 perpendicular_residuals = Monkeys(m).Sessions(sessn).(resid_id);
%                 
%                 % control: restrict to only neurons with low num VR imgs
% %                 units_to_use = (cat_props/2 + dog_props/2) < 0.25;
% %                 perpendicular_residuals = perpendicular_residuals(units_to_use);
%                 
%                 [f, x] = ecdf(perpendicular_residuals);
%                 xlabel_str = 'Residuals from major axis regression';
                
%                 switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
%                     case 'Base'
% %                         col = [0.4 0.4 0.4];
%                         col = epoch_colors(1,:);
%                         thk = 1;
%                     case 'Pre'
%                         %                         col = 'r';
%                         col = [0.4 0.4 0.4];
%                         col = epoch_colors(2,:);
%                         thk = 2;
%                         pre_data = perpendicular_residuals;
%                     case 'Post'
% %                         col = 'g';
%                         col = epoch_colors(3,:);
%                         thk = 2;
%                         post_data = perpendicular_residuals;
%                     case 'Sub'
%                         col = [0.8 0.8 0.8];
%                         thk = 1;
%                 end

                % ***********************************
                % End differences between diffs and resids
                % ***********************************
                
                plot(x, f, 'Color', col, 'LineWidth', thk, 'DisplayName', Monkeys(m).Sessions(sessn).ShortName)
            end
            
            [pval,h] = ranksum(pre_data, post_data, 'Alpha', ranksum_alpha, 'Tail', 'left');
            fprintf('%s, %s, %d to %d, pre vs. post ranksum (one-sided non-parametric), hyp %d, pval of %0.4f \n', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), h, pval)
            
            [h, pval] = ttest2(pre_data, post_data, 'tail', 'left');
            fprintf('%s, %s, %d to %d, pre vs. post ttest2 (one-sided), hyp %d, pval of %0.4f \n', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), h, pval)
            
            %             title(sprintf('%s, %s, %d to %d, VR alpha %0.2f', Monkeys(m).Name, areas_to_plot{a}, test_int(1), test_int(2), vr_alpha), 'Interpreter', 'none')
            %             legend
            xlabel(xlabel_str)
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

%% Plot slopes of selected MARs

vr_alpha = 0.05;
alpha_string_scale_factor = 100;
% test_intervals = {[175 350]};
test_intervals = {[175 275]};
area_to_plot = 'te';

% make fig
figure2('Position', [400 400 300 800]);
hold on
for m = 1:length(Monkeys)
    subplot(2,1,m)
    hold on
    
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
                        scatter(mean([pre_i_val post_i_val]), t(1)*1.075, '*k', 'HandleVisibility', 'off')
                        plot([pre_i_val post_i_val], repelem(t(1)*1.05, 1, 2) , '-k', 'LineWidth', 2, 'HandleVisibility', 'off')
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
%     xticklabels(get_xtick_labs_colored(xticklabs, short_names, epoch_colors))
    xtickangle(45)
    xticklabels(xticklabs)
    xlim([0.5 0.5 + length(sessions_to_use)])
    ylabel('Slope (degrees)')
    ylim([39 52])
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
% saveas(gcf, fullfile(figureSavePath, sprintf('VR_slopes_%s', m, sessn, slope_id)), 'epsc')

%% Plot correlation of selected MARs

vr_alpha = 0.05;
alpha_string_scale_factor = 100;
% test_intervals = {[175 350]};
test_intervals = {[175 275]};
area_to_plot = 'te';

% make fig
figure2('Position', [400 400 300 800]);
hold on
for m = 1:length(Monkeys)
    subplot(2,1,m)
    hold on
    
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
        rval_id = get_good_interval_name2(test_int, 'full', sprintf('VisPropMAR_RValue_%s_Alpha%d', areas_to_plot{a}, vr_alpha*alpha_string_scale_factor));
        
        % pre allocate
        rvals = zeros(1,length(sessions_to_use));
        
        % collect data
        for i = 1:length(sessions_to_use)
            sessn = sessions_to_use(i);
            rvals(i) = Monkeys(m).Sessions(sessn).(rval_id);
            if ~isempty(regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post', 'once'))
                post_data = rvals(i);
            end
        end
        
        % plot data
        plot(baseline_i_vals, rvals(baseline_i_vals), 'LineWidth', 2, ...
            'DisplayName', sprintf('%s, %d to %d', Monkeys(m).Name, test_int(1), test_int(2)),...
            'Color', matlab_colors(1,:))
        errorbar([pre_i_val post_i_val], rvals([pre_i_val post_i_val]), diff(theta_bounds(:, [pre_i_val post_i_val]))/2, 'LineWidth', 2, ...
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
                        scatter(mean([pre_i_val post_i_val]), t(1)*1.075, '*k', 'HandleVisibility', 'off')
                        plot([pre_i_val post_i_val], repelem(t(1)*1.05, 1, 2) , '-k', 'LineWidth', 2, 'HandleVisibility', 'off')
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
%     xticklabels(get_xtick_labs_colored(xticklabs, short_names, epoch_colors))
    xtickangle(45)
    xticklabels(xticklabs)
    xlim([0.5 0.5 + length(sessions_to_use)])
    ylabel('Slope (degrees)')
    ylim([39 52])
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
% saveas(gcf, fullfile(figureSavePath, sprintf('VR_slopes_%s', m, sessn, slope_id)), 'epsc')

%% Functions
function formatPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',20, 'YGrid', 'on', 'XGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end