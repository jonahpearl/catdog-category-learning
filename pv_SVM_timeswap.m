%% Load behavioral data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
svmRecordPath = 'XMA2/Monkey_structs/SVM_Records.mat';
pv_path = 'XMA2/Monkey_structs';
fullSVMPath = fullfile(EXT_HD, pv_path, 'SVM_results_%g.mat');

% Load behavioral data
load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_behav_and_metaNI.mat')) % behavior and neural summaries, but w/o spike times

%% Load the compiled timeswap data

in = load('/Users/jonahpearl/Documents/BJR group/Catdog_paper/timeswap_SVMs/compiled_results.mat');
AllResults = in.AllResults;

%% Other params

rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)
start_times = [-100:5:450];

%% Plot heatmaps

cax_by_monk = {[0.5, 0.8], [0.5, 0.65]};
for m = 1:length(AllResults)
    kfls = AllResults(m).Compiled_KFLs;
    sessions_to_use = rSessionsByMonk{m};

    figure('Position', [580, 640, 950, 400])
%     figure
    for i = 1:length(sessions_to_use)

        mean_kfl = mean(kfls(:,:,:,i), 3);
        to_plot = 1 - mean_kfl;
%         to_plot = permute(to_plot, [2,1,3]);

        subplot(1,2,i)
        hold on
        imagesc(to_plot, 'Interpolation', 'nearest')
        plot(0:111, 0:111, 'k--')
        
        set(gca, 'YDir','reverse')
        colormap('viridis')
        caxis(cax_by_monk{m})
%         colorbar()
        title(sprintf('%s', i))
        xlabel('Training bin')
        ylabel('Testing bin')
        xticks(1:20:111)
        xticklabels(start_times(1:20:end))
        yticks(1:20:111)
        yticklabels(start_times(1:20:end))
        formatSVMPlot(gca, gcf)
    end
end

%%
formatSVMPlot(gca, gcf)

%% Functions
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
