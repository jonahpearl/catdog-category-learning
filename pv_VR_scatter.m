% analyze visual responsiveness of units to the cat/dog images

%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';
figureSavePath = '/Users/jonahpearl/Documents/BJR group/Catdog paper/';
fname = 'MaxMarta_VR_img_TTest.mat';
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

%% (Prop VR) Use major-axis regression to assess slopes, store residuals
% MAR alpha is not dynamic for figures -- if you want to use a different
% MAR alpha (for slope CI's), come back here and re-run the analysis first.

vr_alpha = 0.05;
mar_alpha = 0.05;
alpha_string_scale_factor = 100;
test_intervals = {[75 175], [175 350]};
baseline_intervals = {[-150 -50], [-175 0]};
% test_intervals = {[175 350]};
% baseline_intervals = {[-175 0]};
make_plots = false;


for m = 1:length(Monkeys)
    sessions_to_use = 1:length(Monkeys(m).Sessions);
    
    for p = 1:length(test_intervals)
        test_int = test_intervals{p};
        baseline_int = baseline_intervals{p};
%         tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img'); % using non-parametric testing
        tid = get_good_interval_name2(test_int, 'full', 'VisResp_test_img_TTEST'); % with a t-test
        bid = get_good_interval_name2(baseline_int, '', '');
        vr_id = strcat(tid,bid);
        
        areas_to_plot = unique({Monkeys(m).Sessions(sessions_to_use(1)).UnitInfo.Area});
        for a = 1:length(areas_to_plot)
            if make_plots
                figure('Position', [400 400 1300 1000])
                hold on
            end
            for i = 1:length(sessions_to_use)
                sessn = sessions_to_use(i);
                
                % Get data
                data_mat = Monkeys(m).Sessions(sessn).(vr_id); % stimuli x units
                area_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, areas_to_plot{a});
                data_mat = data_mat(:, area_bool);
                cat_props = sum(data_mat(1:260,:) < vr_alpha) / 260;
                dog_props = sum(data_mat(261:520,:) < vr_alpha) / 260;
                
                % Calculate correlation
                x = corrcoef(cat_props, dog_props);
                fprintf('%s, %s, %s, interval %d to %d, correlation %0.2f \n', ...
                    Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName, areas_to_plot{a}, test_int(1), test_int(2), x(1,2))
                
                
                % Get angle of best-fit line with major axis regression
                % method. I checked this function against a reference, it's
                % correct.
                [coeffs,coeff_ints,~,~] = maregress(dog_props, cat_props, mar_alpha);
                theta = rad2deg(atan(coeffs(2)));
                theta_bounds = rad2deg(atan(coeff_ints(2,:)));
                intercept = coeffs(1);
                intercept_bounds = coeff_ints(1,:);
                
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
                Monkeys(m).Sessions(sessn).(propn_id) = [cat_props; dog_props]; % 2 x units
                Monkeys(m).Sessions(sessn).(resid_id) = perpendicular_residuals; % 1 x units
                Monkeys(m).Sessions(sessn).(slope_id) = [theta_bounds(1) theta theta_bounds(2)]; % 1 x 3
                Monkeys(m).Sessions(sessn).(intercept_id) = [intercept_bounds(1) intercept intercept_bounds(2)]; % 1 x 3
                
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

interval_to_plot = [175 350];
vr_alpha = 0.05;
alpha_string_scale_factor = 100;
area_to_plot = 'te';
color_by_waveform_class = false;
colors = cbrewer('qual', 'Set1', 5);
f1 = figure('Position', [400 400 860 800]);
hold on
f2 = figure('Position', [400 400 860 800]);
% figure('Position', [400 400 1600 500])
hold on

for m = 1:length(Monkeys)
    sessions_to_use = Monkeys(m).Sessions_to_use;
   
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        
        % Get data
        propn_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_CatDogPropns_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        slope_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_Slope_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        intercept_id = get_good_interval_name2(interval_to_plot, 'full', sprintf('VisPropMAR_Intercept_%s_Alpha%d', area_to_plot, vr_alpha*alpha_string_scale_factor));
        t = Monkeys(m).Sessions(sessn).(propn_id);
        cat_props = t(1,:);
        dog_props = t(2,:);
        slope = tan(deg2rad(Monkeys(m).Sessions(sessn).(slope_id)(2)));
        intercept = Monkeys(m).Sessions(sessn).(intercept_id)(2);
        xvals = 0:0.01:1;
        
        % Get color / plotting params
        switch regexp(Monkeys(m).Sessions(sessn).ShortName, '([^0-9-]*)', 'match', 'once')
            case 'Base'
                col = epoch_colors(1,:);
                thk = 1;
            case 'Pre'
                col = epoch_colors(2,:);
                thk = 2;
                pre_data = perpendicular_residuals;
            case 'Post'
                col = epoch_colors(3,:);
                thk = 2;
                post_data = perpendicular_residuals;
        end
                
        % Make scatter plot
        set(0, 'CurrentFigure', f1)
        subplot(length(Monkeys),length(sessions_to_use),...
            length(sessions_to_use)*(m-1) + mod(i-1, length(sessions_to_use)) + 1)
        hold on
        scatter(dog_props, cat_props, 'bo')
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
        subplot(length(Monkeys),length(sessions_to_use),...
            length(sessions_to_use)*(m-1) + mod(i-1, length(sessions_to_use)) + 1)
        histogram(dog_props - cat_props, 'BinEdges', -0.2:0.025:0.2)
        ax = gca;
        formatPlot(ax, f2)
        ax.YAxis.Visible = 'off';
        ax.YGrid = 'off';
        ax.XGrid = 'off';
        
        
    end
    %     sgtitle(sprintf('Area %s, %d to %d, image VR propns,(alpha = %0.2f)',...
    %         area_to_plot, interval_to_plot(1), interval_to_plot(2), vr_alpha), 'Interpreter', 'none')
end

saveas(f1, fullfile(figureSavePath, sprintf('VR_scatter_%s.svg', propn_id)))
saveas(f2, fullfile(figureSavePath, sprintf('VR_hist_%s.svg', propn_id)))

%% Functions
function formatPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',22, 'YGrid', 'on', 'XGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end