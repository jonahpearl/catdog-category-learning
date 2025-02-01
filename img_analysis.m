% analyzes the cat/dog images themselves

%% Re-load Images struct
clearvars
% close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
figureSavePath = '/Users/jonahpearl/Documents/BJR group/Catdog paper/';
imgPath = 'Images';
fname = 'Catdog_260_Images_RgbSilh.mat';
load(fullfile(EXT_HD, imgPath, fname), 'Images');

%% Compare cat/dog colors

% Strategy: use silhouettes to extract foreground pixels. Then histograms.

% Parameters
% channels = {'red', 'green', 'blue'};
channels = {'Hue', 'Saturation', 'Value'};

% Pre-allocate
cat_cols = zeros(260*(256^2), 3); % (num imgs * num pixels) x num channels (rgb)
dog_cols = zeros(260*(256^2), 3);

% Iterators
iic = 1;
iid = 1;

for i = 1:length(Images)
    numpix = numel(Images(i).Silhouette);
    
    % Pre-allocate for the image
    rgb = zeros(numpix, 3);
    
    % Reshape into 1 x 256^2 vectors
    silh = reshape(Images(i).Silhouette, 1, numpix);
    for c = 1:3
       rgb(:,c) = reshape(Images(i).RGBImage(:,:,c), 1, numpix); 
    end
    
    % Turn rgb into hsv
    hsv = rgb2hsv(rgb ./ 255);
    
    % Find foreground pixels
    silh_inds = find(silh == 0);
    
    % Store foreground color values
    for c = 1:3
        if Images(i).Category == 1    
            cat_cols(iic:(iic+numel(silh_inds)-1), c) = hsv(silh_inds,c);
        elseif Images(i).Category == 2
            dog_cols(iid:(iid+numel(silh_inds)-1), c) = hsv(silh_inds,c);
        end
    end
    
    % Iterate
    if Images(i).Category == 1    
        iic = iic + numel(silh_inds);
    elseif Images(i).Category == 2
        iid = iid + numel(silh_inds);
    end
end

% Remove excess zeros
cat_cols(iic:end, :) = [];
dog_cols(iid:end, :) = [];


% Plot the histograms
figure2
for c = 1:3
    subplot(3,1,c)
    hold on
    histogram(cat_cols(:,c), 'NumBins', 50, 'Normalization', 'probability')
    histogram(dog_cols(:,c), 'NumBins', 50, 'Normalization', 'probability')
    xlabel(sprintf('%s', channels{c}))
    ylabel('Probability')
    if c == 1
        legend({'Cats', 'Dogs'})
    end
    formatPlot(gca, gcf)
    set(gcf,'renderer','Painters')
end

%% Functions
function formatPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize', 32, 'YGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end

