% Analyzes basic visual properties of the cat/dog stimuli

%% Read in all images
% clearvars
% close all
% 
% stimset_path = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/Images/260-260/';
% imgs_to_ignore = [".",  "..", ".DS_Store"];
% imgs = dir(stimset_path);
% imgs = string({imgs.name});
% Images = struct();
% Images(length(imgs)-3).FilePath = [];
% Images(length(imgs)-3).RGBImage = [];
% ii = 1;
% for i = 1:length(imgs)
%     img_name = imgs(i);
%     if any(strcmp(img_name, imgs_to_ignore))
%         continue
%     end
%     img_path = strcat(stimset_path, img_name);
%     img_in = imread(img_path); 
%     Images(ii).FilePath = img_path;
%     
%     if ismatrix(img_in) == 1 % catch grayscale imgs
%         img_in = cat(3, img_in, img_in, img_in);
%     end
%         
%     Images(ii).RGBImage = img_in;
%     if ~isempty(regexp(img_name, 'cat00', 'once'))
%         Images(ii).Category = 1;
%         Images(ii).CategoryStr = "cat";
%     elseif ~isempty(regexp(img_name, 'cat01', 'once'))
%         Images(ii).Category = 2;
%         Images(ii).CategoryStr = "dog";
%     end
%     ii = ii + 1;
% end

%% Get silhouettes and total non-background pixels
% 
% noOfColors = 4;
% 
% for i = 1:length(Images)
%     img = Images(i).RGBImage;
%     
%     % set up input
%     img_size = size(img);
%     r = img(:,:,1);
%     g = img(:,:,2);
%     b = img(:,:,3);
%     inputImg = zeros((img_size(1) * img_size(2)), 3);
%     inputImg(:,1) = r(:);
%     inputImg(:,2) = g(:);
%     inputImg(:,3) = b(:);
%     inputImg = double(inputImg);
%     
%     % run kmeans
%     [idx, C] = kmeans(inputImg, noOfColors,...
%         'EmptyAction', 'singleton');
% 
%     % round output points to ints, so they can be colors
%     palette = round(C);
%     
%     % find all indexes of the brightest color (will be the white
%     % background)
%     white_idx = find(sum(palette,2) == max(sum(palette,2)));
%     
%     % reduce palette to black and white
%     reduced_palette = palette;
%     for j = 1:size(palette,1)
%         if j == white_idx
%             reduced_palette(j,:) = [255 255 255];
%         else
%             reduced_palette(j,:) = [0 0 0];
%         end
%     end
%     
%     % set up output img matrix
%     outImg_bw = zeros(img_size(1), img_size(2),3);
%     outImg_color = zeros(img_size(1), img_size(2),3);
%     temp = reshape(idx, [img_size(1) img_size(2)]);
%     for j = 1:img_size(1)
%         for k = 1:img_size(2)
%             outImg_color(j,k,:) = palette(temp(j,k),:); % for normal kmeans image
%             outImg_bw(j,k,:) = reduced_palette(temp(j,k),:); % for BW image
%         end
%     end
%     
%     % convert image out
%     outImg_bw = uint8(outImg_bw);
%     outImg_bw = rgb2gray(outImg_bw);
%     
%     % debugging
% %     imshow(outImg_bw, 'InitialMagnification', 400);
%     
%     % Invert for hole filling. Catches most errors.
%     outImg_bw = (255 - outImg_bw);
%     outImg_bw = imfill(outImg_bw, 'holes');
%     outImg_bw = (255 - outImg_bw);
%     
%     % debugging
%     imshow(uint8(outImg_color));
%     pause(0.2)
%     
%     % store
%     Images(i).Silhouette = outImg_bw;
%     palid = sprintf('KMeansPalette_%d', noOfColors);
%     indid = sprintf('KMeansInds_%d', noOfColors);
%     Images(i).(palid) = palette;
%     Images(i).(indid) = idx;
%     fprintf('Done with image %d \n', i)
% end

%% Save silhouettes
% outpath = '/Users/pearlje/Documents/MATLAB/matsumoto/Images/Silhouette_260/';
% 
% for i = 1:length(Images)
%     img = Images(i).Silhouette;
%     [~, name,~] = fileparts(Images(i).FilePath);
%     imwrite(img, strcat(outpath, sprintf('%s.png', name)));
% end

%% Save Images struct
% outpath = '/Users/pearlje/Documents/MATLAB/matsumoto/Images/';
% fname = 'Catdog_260_Images_RgbSilh.mat';
% save(fullfile(outpath, fname), 'Images')

%% Re-load Images struct

clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
figureSavePath = '/Users/jonahpearl/Documents/BJR group/Catdog paper/';
imgPath = 'Images';
fname = 'Catdog_260_Images_RgbSilh.mat';
load(fullfile(EXT_HD, imgPath, fname), 'Images');

%% Scatter cat/dog colors
n = 10;
figure
hold on
grid on
cat_mat = vertcat(Images(1:20).KMeansPalette_6);
dog_mat = vertcat(Images(21:40).KMeansPalette_6);

colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3    ;
    0          0        0    ];
cat_color = colors(1,:);
dog_color = colors(2,:);

zeros_vector = repelem(0, size(cat_mat,1));
scatter3(cat_mat(:,1), cat_mat(:,2), cat_mat(:,3), [], cat_color, 'o')
scatter3(cat_mat(:,1), cat_mat(:,2), zeros_vector, 10, cat_color, 'o', 'filled')
scatter3(cat_mat(:,1), zeros_vector, cat_mat(:,3), 10, cat_color, 'o', 'filled')
% scatter3(zeros_vector, cat_mat(:,2), cat_mat(:,3), 10, 'ro', 'filled')


scatter3(dog_mat(:,1), dog_mat(:,2), dog_mat(:,3), [], dog_color, 'o')
scatter3(dog_mat(:,1), dog_mat(:,2), zeros_vector, 10, dog_color, 'o', 'filled')
scatter3(dog_mat(:,1), zeros_vector, dog_mat(:,3), 10, dog_color, 'o', 'filled')
% scatter3(zeros_vector, dog_mat(:,2), dog_mat(:,3), 10, 'bo', 'filled')


xlabel('Red')
xlim([0 256])
xticks([0 256])
ylabel('Green')
ylim([0 256])
yticks([0 256])
zlabel('Blue')
zlim([0 256])
zticks([0 256])
make_custom_patch_legend({cat_color, dog_color}, {'Cats', 'Dogs'}, 'Location', 'northeastoutside')
set(gca, 'FontSize', 42)
view([-136.8 30.6])

%% Train classifier on stim colors
% accuracy about 55 %
X = vertcat(cat_mat, dog_mat);
Y = vertcat(ones(size(cat_mat,1),1), 2*ones(size(dog_mat,1),1));
loss_vec = zeros(1, n);
for i = 1:n
    model = fitcsvm(X, Y, 'Standardize', true);
    cv = crossval(model);
    loss_vec(i) = kfoldLoss(cv);
end

%% Plot cats and dogs
% cats | dogs
% 10x2 | 10 x 2
% 10x24 | 10 x 24

% figure('Position', [100 100 1500 1500])
% montage({Images(21:260).RGBImage}, 'Size', [15 16])
% % saveas(gcf, 'cats1.png')
% 
% figure('Position', [100 100 1500 1500])
% montage({Images(281:520).RGBImage}, 'Size', [15 16])
% saveas(gcf, 'dogs1.png')

% show random subsets of the 260 sets for conciseness
rnd_cats = randsample(21:260, 20);
rnd_dogs = randsample(281:520, 20);
figure('Position', [100 100 1500 1500])
montage({Images(rnd_cats).RGBImage}, 'Size', [4 5])
figure('Position', [100 100 1500 1500])
montage({Images(rnd_dogs).RGBImage}, 'Size', [4 5])


figure('Position', [100 100 1000 1000])
montage({Images(1:20).RGBImage}, 'Size', [4 5])
% saveas(gcf, 'cats2.png')

figure('Position', [100 100 1000 1000])
montage({Images(261:280).RGBImage}, 'Size', [4 5])
% saveas(gcf, 'dogs2.png')

%% Compare cat/dog silhouettes

% black pixels are 0, white are 255.

num_pixs = zeros(2,length(Images));

for i = 1:length(Images)
    num_pixs(1,i) = Images(i).Category;
    num_pixs(2,i) = sum(sum(Images(i).Silhouette==0));
end

nb_cat = num_pixs(2, num_pixs(1,:)==1); % number of black pixels in cats
nb_dog = num_pixs(2, num_pixs(1,:)==2); % ditto in dogs
mean_diff = mean(nb_cat) - mean(nb_dog);

% histogram
figure
hold on
histogram(nb_cat, 'Normalization', 'probability', 'FaceColor', 'r', 'FaceAlpha', 0.5)
histogram(nb_dog, 'Normalization', 'probability', 'FaceColor', 'b', 'FaceAlpha', 0.5)
legend({'Cats', 'Dogs'})
xlabel('Number black pixels')
ylabel('Probability')
set(gca, 'FontSize', 42)

% ttest
[h, pval, ci] = ttest2(nb_cat, nb_dog);
fprintf('Hyp %d, pval %0.3f, ci %d to %d. True mean diff %d \n', h, pval, ci(1), ci(2), mean_diff)
fprintf('True mean diff is %0.3f %% of total pixels \n', mean_diff / 256^2 * 100)

%% Compare cat/dog colors (method 1)

channels = {'red', 'green', 'blue'};
cat_cols = zeros(260*(210^2), 3); % (num imgs * overestimate num pixels) x num channels (rgb)
dog_cols = zeros(260*(210^2), 3);
iic = 1;
iid = 1;

for i = 1:length(Images)
    numpix = numel(Images(i).Silhouette);
    
    % pre allocate
    rgb = zeros(numpix, 3);
    
    % reshape into 1 x 256^2 vectors
    silh = reshape(Images(i).Silhouette, 1, numpix);
    for c = 1:3
       rgb(:,c) = reshape(Images(i).RGBImage(:,:,c), 1, numpix); 
    end
    
    % find within-animal pixels
    silh_inds = find(silh == 0);
    
    % store color values
    for c = 1:3
        if Images(i).Category == 1    
            cat_cols(iic:(iic+numel(silh_inds)-1), c) = rgb(silh_inds,c);
        elseif Images(i).Category == 2
            dog_cols(iid:(iid+numel(silh_inds)-1), c) = rgb(silh_inds,c);
        end
    end
    
    % iterate
    if Images(i).Category == 1    
        iic = iic + numel(silh_inds);
    elseif Images(i).Category == 2
        iid = iid + numel(silh_inds);
    end
end

% remove excess zeros
cat_cols(iic:end, :) = [];
dog_cols(iid:end, :) = [];


% color histograms
figure
for c = 1:3
    subplot(3,1,c)
    hold on
    histogram(cat_cols(:,c), 'Normalization', 'probability')
    histogram(dog_cols(:,c), 'Normalization', 'probability')
    xlabel(sprintf('%s value', channels{c}))
    ylabel('Probability')
    if c == 1
        legend({'Cats', 'Dogs'})
    end
end

%% Compare cat/dog colors (method 2)
col_names = {'r', 'g', 'b'};
figure
for c = 1:3
    cid = sprintf('RGBHists_%s', col_names{c});
    for i = 1:length(Images)
        numpix = numel(Images(i).Silhouette);
        silh = reshape(Images(i).Silhouette, 1, numpix);
        silh_inds = find(silh == 0);
        color_channel = reshape(Images(i).RGBImage(:,:,c), 1, numpix);
        Images(i).(cid) = imhist(color_channel(silh_inds));
    end
    
    mean_cat = mean([Images(1:260).(cid)],2);
    mean_dog = mean([Images(261:520).(cid)],2);
    
    subplot(3,1,c)
    hold on
    plot(mean_cat, '-')
    plot(mean_dog, '-')
    title(sprintf('Category-mean color (%s) histograms', col_names{c}))
    xlabel('Pixel value')
    ylabel('Number of pixels')
    legend({'Cats', 'Dogs'})
end

%% Compare absolute luminances

% get luminance histogram for each image
for i = 1:length(Images)
    im = rgb2gray(Images(i).RGBImage);
    Images(i).AbsLum = imhist(im);
end

% average across images in each category
mean_cat = mean([Images(1:260).AbsLum],2);
mean_dog = mean([Images(261:520).AbsLum],2);

% plot as lines, so you can clearly see both hists
figure
hold on
plot(mean_cat, '-')
plot(mean_dog, '-')
set(gca, 'YScale', 'log')
title('Category-mean luminance histograms')
xlabel('Absoluate luminance')
ylabel('Number of pixels')
legend({'Cats', 'Dogs'})
set(gca, 'FontSize', 42)

%% Compare cat/dog relative luminance (silh or whole image)

for i = 1:length(Images)
    numpix = numel(Images(i).Silhouette);
    
    % pre allocate
    lums = zeros(numpix, 1);
    
    % get relative luminance values by pixel (0-1)
    img_sRGB = double(Images(i).RGBImage) ./ 255;
    img_bools = img_sRGB < 0.03938; % different conversions for high or low vals
    img_norm_low = img_sRGB ./ 12.93;
    img_norm_low(~img_bools) = 0;
    img_norm_high = ((img_sRGB+0.055) ./ 1.055) .^ 2.4;
    img_norm_high(img_bools) = 0;
    img_RL_3d = img_norm_low + img_norm_high;
    img_RL = img_RL_3d(:,:,1)*0.2126 + img_RL_3d(:,:,2)*0.7152+ img_RL_3d(:,:,3)*0.0722;
    
    % store mean luminance
    Images(i).MeanLum_overall = mean(mean(img_RL));
    
    % store silhouette-only mean luminance
    silh = reshape(Images(i).Silhouette, 1, numpix);
    rl = reshape(img_RL, 1, numpix);
    Images(i).MeanLum_silh = mean(rl(silh == 0));
end


% overall luminance histograms
cat_lums = [Images(1:260).MeanLum_overall];
dog_lums = [Images(261:520).MeanLum_overall];
mean_diff = mean(cat_lums) - mean(dog_lums);
figure
hold on
histogram(cat_lums, 'Normalization', 'probability', 'FaceColor', 'r', 'FaceAlpha', 0.5)
histogram(dog_lums, 'Normalization', 'probability', 'FaceColor', 'b', 'FaceAlpha', 0.5)
[h, pval, ci] = ttest2(cat_lums, dog_lums);
fprintf('Hyp %d, pval %0.3f, ci %d to %d. True mean diff %d \n', h, pval, ci(1), ci(2), mean_diff)
fprintf('True mean diff is %0.3f %% of total pixels \n', mean_diff / 256^2 * 100)
title('Whole-image luminance')
xlabel('Image-mean relative luminance')
ylabel('Probability')
legend({'Cats', 'Dogs'})

% just-animal luminance histograms
% cat_lums = [Images(1:260).MeanLum_silh];
% dog_lums = [Images(261:520).MeanLum_silh];
% figure
% hold on
% histogram(cat_lums, 'Normalization', 'probability')
% histogram(dog_lums, 'Normalization', 'probability')
% % title('Silhouette-only luminance')
% xlabel('Image-mean relative luminance')
% ylabel('Probability')
% legend({'Cats', 'Dogs'})
% set(gca, 'FontSize', 42)

%% Convert to grey

for i = 1:length(Images)
    Images(i).Greyscale = rgb2gray(Images(i).RGBImage);
end

%% Compute FFTs of each image, and average across categories
% Images must be square!

% params
numpix = 256;
numim = length(Images);

% pre allocate
% angs = zeros(numpix, numpix, numim);
% mags = zeros(numpix, numpix, numim);
energies = zeros(numim, ceil(sqrt(2)*numpix/2));
cropped_energies = zeros(numim, floor(numpix/2));

% make meshgrid for r
f1 = -numpix/2:numpix/2-1;
[XX YY] = meshgrid(f1,f1);
[t r] = cart2pol(XX,YY); % polar positions for each pixel in the image
r = round(r); % round so we can use it to index

for i = 1:length(Images)
    
    % shifted fft of image
    fftim1 = fftshift(fft2(Images(i).Greyscale));
    
    % separate phases and magnitudes
%     [angs(:,:,i),mags(:,:,i)] = cart2pol(real(fftim1),imag(fftim1));
    [~, fftim] = cart2pol(real(fftim1),imag(fftim1));
    
    % take sum of energies across all spatial orientations at each frequency.
    % for each value of r (ie, each possible radius in fft space, ie, each (inverse) frequency)
    % add up all the fft magnitudes at that frequency+1.
    % ie, energies(j,i) = sum(fftim(r(:) == j+1));
    energies(i,:) = accumarray(r(:)+1,fftim(:));  % 1 x num freq
    
    % remove frequencies for which sampling is limited (ie, anything larger
    % than the radius of the largest inscribed circle in the image. Eg, the
    % lowest possible frequency is the diagonal of the image (longest
    % radius from center), but there are only four samples for it (the
    % diagonals).
    cropped_energies(i,:) = energies(i, 1:floor(numpix/2)); 
    
end

%% Plot FFT results
cat_color = 'r';
dog_color = 'b';

% compare cats to dogs
figure
hold on
% errorbar(1:size(cropped_energies,2), mean(cropped_energies(1:260, :)),  std(cropped_energies(1:260, :)) / sqrt(260))
% errorbar(1:size(cropped_energies,2), mean(cropped_energies(261:520, :)),  std(cropped_energies(261:520, :)) / sqrt(260))
stdshade_sem(cropped_energies(1:260, :), 0.5, cat_color)
stdshade_sem(cropped_energies(261:520, :), 0.5, dog_color)
set(gca, 'YScale', 'log', 'XScale', 'log')
xlabel('Spatial Frequency (cpi)')
ylabel('Relative Energy')
% title('FFT analysis of image categories')
set(gca, 'FontSize', 24)
