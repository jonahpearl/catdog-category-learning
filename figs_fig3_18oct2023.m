% figs for the cat/dog manuscript, 18 Oct 2023

%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';
figureSavePath = '/Users/jonahpearl/Documents/BJR group/Catdog paper/';
fname = 'MaxMarta_VR_img_TTest_jun2021.mat'; % changed ttest2 to ttest (because they're paired!)
Data = load(fullfile(EXT_HD, pv_path, fname));
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Misc. prep 
epoch_colors = [0.5 0.5 0.5; mlc(1:2)];
matlab_colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3;
    0          0        0    ];

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

%% Fig 3 params

vr_alpha = 0.05;
% vr_alpha = 0.01;  % if want to try more VR's, need to go back and re-run
% the analysis.
mar_alpha = 0.05;
test_intervals = {[175 275]};
baseline_intervals = {[-150 -50]};

% for retrieving data from different alpha vals -- dont change
alpha_string_scale_factor = 100;

rSessionsByMonk = {[7 9] [6 7]}; % pre/post
% rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};

%% Stats for text
area = 'te';
chi_sq_alpha = 0.05;
for m = 1:length(Monkeys)
%     sessions_to_use = 1:length(Monkeys(m).Sessions);
     sessions_to_use = rSessionsByMonk{m};
    for iInt = 1:length(test_intervals)
        test_int = test_intervals{iInt};
        baseline_int = baseline_intervals{iInt};
        tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img_TTEST'); % with a t-test
        bid = get_good_interval_name2(baseline_int, '', '');
        vr_id = strcat(tid,bid);
        baseline_props = [];
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
            
            % Compare cat and dog props
%             fprintf('%s, %s, %d to %d, %2g images is mean betw-catg diff \n',...
%                 Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName,...
%                 test_int(1), test_int(2), mean(cat_props - dog_props)*260);
            
            
            % Find proportion of units responding signf. to at least some fraction
            % (10 or 25%) of the images
            n_units = size(data_mat,2);
            overall_props = sum(data_mat < vr_alpha) / 520;
            
            x_10 = sum(overall_props > 0.10);
%             fprintf('%s, %s, %d to %d, %2g%% of units respond to more than 10%% of images \n',...
%                 Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName,...
%                 test_int(1), test_int(2), round(x_10/n_units*100,1));
            
            x_25 = sum(overall_props > 0.25);
%             fprintf('%s, %s, %d to %d, %2g%% of units respond to more than 25%% of images \n',...
%                 Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName,...
%                 test_int(1), test_int(2), round(x_25/n_units*100,1));
            
            % Prepare for chi-sq test, which needs [X,N] for each condition
            switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
            case 'Pre'
                pre_props = {[x_10 n_units], [x_25 n_units]};
            case 'Post'
                post_props = {[x_10 n_units], [x_25 n_units]};
            case 'Base'
                baseline_props = [baseline_props x_10/n_units];
            end
        end
        
        
%         fprintf('%s, Baseline, %2g +/- %2g units respond to more than 10%% of imgs \n\n',...
%             Monkeys(m).Name, round(mean(baseline_props)*100,1), round(std(baseline_props)*100,1))
        
        % Do stats
        [h,pval, chi2stat, df] = prop_test(...
            [pre_props{1}(1) post_props{1}(1)],...
            [pre_props{1}(2) post_props{1}(2)],...
            false, chi_sq_alpha);
        fprintf('%s, PRE vs POST, %d to %d, chi-sq for 10%% level: p=%3g, chisq=%3g \n', Monkeys(m).Name, test_int(1), test_int(2), pval, chi2stat)
        
        [h,pval, chi2stat, df] = prop_test(...
            [pre_props{2}(1) post_props{2}(1)],...
            [pre_props{2}(2) post_props{2}(2)],...
            false, chi_sq_alpha);
        fprintf('%s, PRE vs POST, %d to %d, chi-sq for 25%% level: p=%3g, chisq=%3g \n', Monkeys(m).Name, test_int(1), test_int(2), pval, chi2stat)
            
    end
end
        
%% Fig 3B/C prep: run major-axis regression to assess slopes, store residuals
% MAR alpha is not dynamic for figures -- if you want to use a different
% MAR alpha (for slope CI's), come back here and re-run the analysis first.
make_plots = false;
areas_to_plot = {'te', 'anterior', 'middle', 'posterior'};

for m = 1:length(Monkeys)
    sessions_to_use = 1:length(Monkeys(m).Sessions);
    
    for iInt = 1:length(test_intervals)
        test_int = test_intervals{iInt};
        baseline_int = baseline_intervals{iInt};
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

%% Fig 3B/C: plot selected MAR results with scatterplots

area_to_plot = 'te';
f1 = figure('Position', [400 400 860 800]); % scatter plots
f2 = figure('Position', [400 400 400 800]); % histograms
hold on
interval_to_plot = [175 275];

for m = 1:length(Monkeys)
    sessions_to_use = rSessionsByMonk{m};
    diffs = cell(1,length(sessions_to_use));
    sums = cell(1,length(sessions_to_use));
    resids = cell(1,length(sessions_to_use));
    baseline_diffs = {};
    
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
                baseline_diffs = [baseline_diffs diffs{i}];
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
        scatter(dog_props, cat_props, 30,'bo') % size 70 for main figs, 30 for supps
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
%         title(Monkeys(m).XTickLabs{i},'FontWeight', 'normal')
        formatPlot(gca, f1)
        axis square

        % Plot reference lines
        set(0, 'CurrentFigure', f1)
        plot(xvals, xvals, 'k--', 'LineWidth', 2)
        plot(xvals, intercept + xvals*slope, 'r--', 'LineWidth', 2)
        
        % Make histograms to visualize stats
        bins = 0:0.01:0.2;
        set(0, 'CurrentFigure', f2)
        subplot(length(Monkeys), 1, m)
        hold on
        xlabel('Selectivity |cat-dog|')
        ylabel('Probability')
        histogram(abs(dog_props - cat_props), ...
                'Normalization', 'probability',...
                'BinEdges', bins,...
                'FaceColor', mlc(i))
        xl = get(gca,'xlim');
        xlim([-0.02 xl(2)])
        
        set(gca, 'FontSize', 20);
%         ax = gca;
%         ax.YAxis.Scale = 'log';  % for insets?
        
        mu = mean(abs(dog_props - cat_props));
        sigma = std(abs(dog_props - cat_props));
%         sem = std(abs(dog_props - cat_props)) / sqrt(length(dog_props));
        scatter(mu, 1 + i/20, 100, mlc(i), 'v', 'filled', 'MarkerEdgeColor', 'k')
        plot([mu-sigma mu+sigma], repelem(0.9 + i/20,2), '-', 'Color', mlc(i), 'LineWidth', 1.5)
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
    fprintf('ABS DIFFS pre vs post ttest2: h = %d, pval = %0.4f \n', ...
            h, pval)
        
    fprintf('Mean+std ABS DIFFS: \n\t pre: %0.2g%% +/- %0.2g%% \n\t post: %0.2g%% +/- %0.2g%% \n', ...
        mean(abs(pre_diffs))*100, std(abs(pre_diffs))*100,...
        mean(abs(post_diffs))*100, std(abs(post_diffs))*100)
        
%     baseline_means = cellfun(@(x) mean(abs(x)), baseline_diffs);
%     fprintf('Mean+std of baseline mean diffs: %0.2g%% +/- %0.2g%%\n',...
%         mean(baseline_means)*100, std(baseline_means)*100)
    
    % Test sums (did units responsd to more images overall?)
%     [h,pval] = ttest2(sums{pre_diff_idx}, sums{post_diff_idx});
%     fprintf('SUMS, pre: %0.0f, post: %0.0f, pre vs post ttest2: h = %d, pval = %0.2f \n', ...
%             mean(sums{pre_diff_idx}), mean(sums{post_diff_idx}), h, pval)
    
    
    % Ditto on residuals
%     [h,pval] = vartest2(resids{pre_diff_idx}, resids{post_diff_idx});
%     fprintf('RESIDUALS Pre var: %0.3g, post var: %0.3g. H = %d, p = %d \n', ...
%             var(resids{pre_diff_idx}), var(resids{post_diff_idx}), h, pval)

end

%% Fig 3D: plot slopes of MARs
area_to_plot = 'te';

% make fig
figure2('Position', [400 400 300 800]);
hold on
for m = 1:length(Monkeys)
    subplot(2,1,m)
    hold on
    
    sessions_to_use = rSessionsByMonk{m};
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
    
    for iInt = 1:length(test_intervals)
        test_int = test_intervals{iInt};
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
            'Color', 'k')
        errorbar([pre_i_val post_i_val], theta_means([pre_i_val post_i_val]), diff(theta_bounds(:, [pre_i_val post_i_val]))/2, 'LineWidth', 2, ...
            'DisplayName', sprintf('%s, %d to %d', Monkeys(m).Name, test_int(1), test_int(2)),...
            'Color', 'k')
        
        % stats testing
%         for i = 1:length(sessions_to_use)
%             sessn = sessions_to_use(i);
%             switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
%                 case {'Pre'}
%                     t = Monkeys(m).Sessions(sessn).(slope_id);
%                     test_bounds = t([1 3]);
%                     post_bounds = post_data([1 3]);
%                     % look for overlap in the confidence intervals.
%                     if isempty(intersect(round(test_bounds(1),1):0.1:round(test_bounds(2),1), round(post_bounds(1),1):0.1:round(post_bounds(2),1)))
%                         scatter(mean([pre_i_val post_i_val]), t(1)*1.075, '*k', 'HandleVisibility', 'off')
%                         plot([pre_i_val post_i_val], repelem(t(1)*1.05, 1, 2) , '-k', 'LineWidth', 2, 'HandleVisibility', 'off')
%                     end
%                 case 'Post'
%                     % dont test against itself
%                     continue
%             end
%         end
    end
    
    % format plot
    xticks(1:length(sessions_to_use))
%     xticklabs = Monkeys(m).XTickLabs;
%     xtickangle(45)
%     xticklabels(xticklabs)
    xlim([0.5 0.5 + length(sessions_to_use)])
    ylabel('Slope (degrees)')
    ylim([39 52])
    yticks(40:5:50)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
        'fontsize', 26, ...
        'fontname', 'Helvetica', 'fontweight', 'normal', ...
        'XColor', 'black', 'YColor', 'black')
    set(gcf, 'Color', 'w')
end


%% Functions

function formatPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',15, 'YGrid', 'on', 'XGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end
