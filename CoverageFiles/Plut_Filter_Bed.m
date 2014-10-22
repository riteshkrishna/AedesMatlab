%% A cleaner version exists in Plut_Analyze_BEDFile.m

% Plutella Data
male_datafile = 'DataFiles/DistanceBased/PLUTELLA/scf/PLUTMALE.sorted.bed';
M = load_data_from_BED(male_datafile);    

female_datafile = 'DataFiles/DistanceBased/PLUTELLA/scf/PLUT_FEMALE_PX8.sorted.bed';
F_1 = load_data_from_BED(female_datafile);
female_datafile = 'DataFiles/DistanceBased/PLUTELLA/scf/PLUT_FEMALE_PX9_1.sorted.bed';
F_2 = load_data_from_BED(female_datafile);
female_datafile = 'DataFiles/DistanceBased/PLUTELLA/scf/PLUT_FEMALE_PX9_2.sorted.bed';
F_3 = load_data_from_BED(female_datafile);
female_datafile = 'DataFiles/DistanceBased/PLUTELLA/scf/PLUT_FEMALE_PX10.sorted.bed';
F_4 = load_data_from_BED(female_datafile);

F_depth = (F_1.depth + F_2.depth + F_3.depth + F_4.depth)/4;
F_baseCovered = (F_1.bases_covered_in_feature + F_2.bases_covered_in_feature + F_3.bases_covered_in_feature + F_4.bases_covered_in_feature) / 4;
F_breadth = (F_1.breadth + F_2.breadth + F_3.breadth + F_4.breadth)/4;
F_contigLength = F_1.contig_length;

FEMALE_DATA = [F_depth F_baseCovered F_contigLength F_breadth];    
MALE_DATA = [M.depth M.bases_covered_in_feature M.contig_length M.breadth];
feature_list = M.features; % Feature list is in the same order for both MALE and FEMALE

% Find features that have zero depth
zero_depth_idx_male = find(MALE_DATA(:,1) == 0);
zero_depth_idx_female = find(FEMALE_DATA(:,1) == 0);

% Find features common between male and female
zero_depth_idx_common = intersect(zero_depth_idx_male, zero_depth_idx_female);

filtered_feature_list = feature_list;
filtered_feature_list(zero_depth_idx_common) = []

Filtered_MALE_DATA = MALE_DATA;
Filtered_MALE_DATA(zero_depth_idx_common,:) = []

Filtered_FEMALE_DATA = FEMALE_DATA;
Filtered_FEMALE_DATA(zero_depth_idx_common,:) = []

MALE_DEPTH = Filtered_MALE_DATA(:,1)
FEMALE_DEPTH = Filtered_FEMALE_DATA(:,1)
FEATURE_LENGTH = Filtered_FEMALE_DATA(:,3) % Same for both sets
FEATURE_LIST = filtered_feature_list;

%% Depth vs Length exp.
figure;
hist([MALE_DEPTH FEMALE_DEPTH])
h = title('Plutella : Depth')
set(h,'FontSize',14)
h = xlabel('Depth (No scaling performed)')
set(h,'FontSize',14)
legend('Male','Female')

% Non-filtered data
figure;
scatter(MALE_DATA(:,3), log10(MALE_DATA(:,1)))
figure;
scatter(FEMALE_DATA(:,3), log10(FEMALE_DATA(:,1)),20,'r')

% Filtered data
figure;
hold on
scatter(FEATURE_LENGTH,log10(MALE_DEPTH))
scatter(FEATURE_LENGTH,log10(FEMALE_DEPTH),20,'r')
hold off


figure;
subplot(1,2,1)
scatter(FEATURE_LENGTH,log10(MALE_DEPTH))
subplot(1,2,2)
scatter(FEATURE_LENGTH,log10(FEMALE_DEPTH),20,'r')
hold off

figure;
hist([log10(FEMALE_DEPTH) log10(MALE_DEPTH)],50)

%% Ratio study
DEPTH_RATIO = MALE_DEPTH ./ FEMALE_DEPTH;
% All 0 values = males zero, all if values = females zero
MALE_ZERO_idx = find(DEPTH_RATIO == 0);
FEMALE_ZERO_idx = find(DEPTH_RATIO == Inf);

% Find all ratios >= 1
RATIO_GT_1_idx = find(DEPTH_RATIO>=1);
% Remove Inf indices from the above
DEPTH_RATIO_1_idx = setdiff(RATIO_GT_1_idx,FEMALE_ZERO_idx);
% Find all ratios >= 2
RATIO_GT_2_idx = find(DEPTH_RATIO>=2);
% Remove Inf indices from the above
DEPTH_RATIO_2_idx = setdiff(RATIO_GT_2_idx,FEMALE_ZERO_idx);

%%
PASSED_RATIOS = DEPTH_RATIO(DEPTH_RATIO_2_idx);
PASSED_RATIOS_features = FEATURE_LIST(DEPTH_RATIO_2_idx)
PASSED_RATIOS_length = Filtered_MALE_DATA(DEPTH_RATIO_2_idx,3)
PASSED_RATIO_MALE_DEPTH = Filtered_MALE_DATA(DEPTH_RATIO_2_idx,1)
PASSED_RATIO_FEMALE_DEPTH = Filtered_FEMALE_DATA(DEPTH_RATIO_2_idx,1)

figure;
scatter(PASSED_RATIO_FEMALE_DEPTH,log10(PASSED_RATIOS_length))
figure;
scatter(PASSED_RATIO_MALE_DEPTH,log10(PASSED_RATIOS_length))

figure;
scatter(MALE_DEPTH(DEPTH_RATIO_2_idx),FEMALE_DEPTH(DEPTH_RATIO_2_idx))
labels = PASSED_RATIOS_features;
text(MALE_DEPTH(DEPTH_RATIO_2_idx),FEMALE_DEPTH(DEPTH_RATIO_2_idx),labels)

%%
%XDATA = log10(MALE_DEPTH);
%YDATA = log10(FEMALE_DEPTH);
%XY_FEATURE_LIST;

% XDATA = log10(MALE_DEPTH(DEPTH_RATIO_2_idx));
% YDATA = log10(FEMALE_DEPTH(DEPTH_RATIO_2_idx));
% XY_FEATURE_LIST = FEATURE_LIST(DEPTH_RATIO_2_idx);

%XDATA = log10(MALE_DEPTH(DEPTH_RATIO_1_idx));
%YDATA = log10(FEMALE_DEPTH(DEPTH_RATIO_1_idx));
XDATA = MALE_DEPTH(DEPTH_RATIO_1_idx);
YDATA = FEMALE_DEPTH(DEPTH_RATIO_1_idx);
XY_FEATURE_LIST = FEATURE_LIST(DEPTH_RATIO_1_idx);

figure;
scatter(XDATA,YDATA);

%% Box method
OUTLIERS = excludedata(XDATA,YDATA,'box',[0 2 0 2])
OUTLIERS_idx = find(OUTLIERS == 0);
figure; scatter(XDATA(OUTLIERS_idx),YDATA(OUTLIERS_idx));
OUTLIERS_FEATURES = XY_FEATURE_LIST(OUTLIERS_idx)

%% Fit method
M_data = XDATA;
F_data = YDATA;
title_label = 'Ratio > 1'
[p,ErrorEst] = polyfit(M_data,F_data,1);
M_hat = polyval(p,M_data,ErrorEst);
resids = F_data - M_hat;

figure;
plot(M_data,F_data,'.')
h = title(title_label);
set(h,'FontSize',14);
h = xlabel('Male');
set(h,'FontSize',14)
h = ylabel('Female');
set(h,'FontSize',14)
    
figure;
plot(M_data, resids,'.')
h = title('Residuals analysis')
set(h,'FontSize',14)
h = xlabel('Male');
set(h,'FontSize',14)
h = ylabel('Residuals');
set(h,'FontSize',14)   

z = 3;
n = size(resids,1);
mu_resid = mean(resids);
std_resid = std(resids);
I = abs(resids - repmat(mu_resid,n,1)) <= z .* repmat(std_resid,n,1);
outlier_idx = find(I == 0);

figure;
scatter(M_data(outlier_idx),F_data(outlier_idx));
h = title('Contigs : Outliers');
set(h,'FontSize',14)
h = xlabel('Male');
set(h,'FontSize',14)
h = ylabel('Female');
set(h,'FontSize',14)
    
sprintf('Total outliers : %d', size(outlier_idx,1))    
outlier_contigs = XY_FEATURE_LIST(outlier_idx);
write_to_file('Plut_outlier_ratio_gt_1.txt',outlier_contigs)
