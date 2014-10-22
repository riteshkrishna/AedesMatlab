%% Cleaned version of Plut_Filter_Bed.m

%% Load Plutella Data
% % Al's assembly
% male_datafile = 'DataFiles/DistanceBased/PLUTELLA/PX_contig/PLUTMALE.sorted.bed';
% M = load_data_from_BED(male_datafile);    
% female_datafile = 'DataFiles/DistanceBased/PLUTELLA/PX_contig/PLUT_FEMALE_PX8.sorted.bed';
% F_1 = load_data_from_BED(female_datafile);
% female_datafile = 'DataFiles/DistanceBased/PLUTELLA/PX_contig/PLUT_FEMALE_PX9_1.sorted.bed';
% F_2 = load_data_from_BED(female_datafile);
% female_datafile = 'DataFiles/DistanceBased/PLUTELLA/PX_contig/PLUT_FEMALE_PX9_2.sorted.bed';
% F_3 = load_data_from_BED(female_datafile);
% female_datafile = 'DataFiles/DistanceBased/PLUTELLA/PX_contig/PLUT_FEMALE_PX10.sorted.bed';
% F_4 = load_data_from_BED(female_datafile);

% US assembly
male_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/PLUTMALE.sorted.bed';
M = load_data_from_BED(male_datafile);    
female_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/PLUT_FEMALE_PX8.sorted.bed';
F_1 = load_data_from_BED(female_datafile);
female_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/PLUT_FEMALE_PX9_1.sorted.bed';
F_2 = load_data_from_BED(female_datafile);
female_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/PLUT_FEMALE_PX9_2.sorted.bed';
F_3 = load_data_from_BED(female_datafile);
female_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/PLUT_FEMALE_PX10.sorted.bed';
F_4 = load_data_from_BED(female_datafile);

%%
% Considering all the female samples
% F_depth = (F_1.depth + F_2.depth + F_3.depth + F_4.depth)/4;
% F_baseCovered = (F_1.bases_covered_in_feature + F_2.bases_covered_in_feature + F_3.bases_covered_in_feature + F_4.bases_covered_in_feature) / 4;
% F_breadth = (F_1.breadth + F_2.breadth + F_3.breadth + F_4.breadth)/4;

% Considering only the good female samples
F_depth = (F_2.depth + F_3.depth)/2;
F_baseCovered = (F_2.bases_covered_in_feature + F_3.bases_covered_in_feature) / 2;
F_breadth = (F_2.breadth + F_3.breadth)/2;

F_contigLength = F_1.contig_length;

FEMALE_DATA = [F_depth F_baseCovered F_contigLength F_breadth];    
MALE_DATA = [M.depth M.bases_covered_in_feature M.contig_length M.breadth];
feature_list = M.features; % Feature list is in the same order for both MALE and FEMALE

%% Clean data
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

%% Ratio study
LOWER_FOLD = 1.5;
UPPER_FOLD = 3;
RATIO_THRESHOLD_LOWER_LIMIT = log2(LOWER_FOLD);
RATIO_THRESHOLD_UPPER_LIMIT = log2(UPPER_FOLD); % Threshold for discarding vary high ratios indices

% DEPTH_RATIO = log2(MALE_DEPTH ./ FEMALE_DEPTH);
% 
% MALE_ZERO_idx = find(DEPTH_RATIO == 0); % All 0 values = males zero 
% FEMALE_ZERO_idx = find(DEPTH_RATIO == Inf); % all inf values = females zero
% 
% RATIO_idx_1 = find(DEPTH_RATIO>=RATIO_THRESHOLD_LOWER_LIMIT); % Find all ratios >= RATIO_THRESHOLD_LOWER_LIMIT
% RATIO_idx_2 = find(DEPTH_RATIO<=RATIO_THRESHOLD_UPPER_LIMIT);
% RATIO_GT_1_idx = intersect(RATIO_idx_1,RATIO_idx_2);
% DEPTH_RATIO_1_idx = setdiff(RATIO_GT_1_idx,FEMALE_ZERO_idx); % Remove Inf indices from the above

DEPTH_RATIO = log2(FEMALE_DEPTH ./ MALE_DEPTH); % ZZ for female > ZW for male
FEMALE_ZERO_idx = find(DEPTH_RATIO == 0); 
MALE_ZERO_idx = find(DEPTH_RATIO == Inf); 
RATIO_idx_1 = find(DEPTH_RATIO>=RATIO_THRESHOLD_LOWER_LIMIT); % Find all ratios >= RATIO_THRESHOLD_LOWER_LIMIT
RATIO_idx_2 = find(DEPTH_RATIO<=RATIO_THRESHOLD_UPPER_LIMIT);
RATIO_GT_1_idx = intersect(RATIO_idx_1,RATIO_idx_2);
DEPTH_RATIO_1_idx = RATIO_GT_1_idx;


%%
%XDATA = log10(MALE_DEPTH);
%YDATA = log10(FEMALE_DEPTH);
%XY_FEATURE_LIST;
% 
XDATA = MALE_DEPTH(DEPTH_RATIO_1_idx);
YDATA = FEMALE_DEPTH(DEPTH_RATIO_1_idx);
XY_LENGTH = FEATURE_LENGTH(DEPTH_RATIO_1_idx);
XY_FEATURE_LIST = FEATURE_LIST(DEPTH_RATIO_1_idx);

figure;
scatter(XDATA,YDATA);

%% Fit method
M_data = XDATA;
F_data = YDATA;
title_label = 'Ratio >= Threshold'
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

z = 1;
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
    
% Sort outliers
[SO,IDX] = sort(resids(outlier_idx))
figure; plot(SO,'.');
outlier_contigs = XY_FEATURE_LIST(outlier_idx(IDX));
sprintf('Total outliers : %d', size(outlier_contigs,1)) 

setdiff(XY_FEATURE_LIST(outlier_idx),outlier_contigs)

out_file = sprintf('Plut_US__outlier_log2_%0.1f_lb_%0.1f_ub.txt',LOWER_FOLD, UPPER_FOLD)
write_to_file(out_file,outlier_contigs)
