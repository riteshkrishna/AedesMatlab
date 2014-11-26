% Script to study the effect of window based coverage for Plutella Data
% US assembly - 10KB window
% male_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/window_10kb/PLUTMALE.sorted.bed';
% M = load_data_from_BED(male_datafile);    
% female_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/window_10kb/PLUT_FEMALE_PX8.sorted.bed';
% F_1 = load_data_from_BED(female_datafile);
% female_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/window_10kb/PLUT_FEMALE_PX9_1.sorted.bed';
% F_2 = load_data_from_BED(female_datafile);
% female_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/window_10kb/PLUT_FEMALE_PX9_2.sorted.bed';
% F_3 = load_data_from_BED(female_datafile);
% female_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/window_10kb/PLUT_FEMALE_PX10.sorted.bed';
% F_4 = load_data_from_BED(female_datafile);

%% US assembly - 1KB window
male_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/window_1kb/PLUTMALE.sorted.bed';
M = load_data_from_BED(male_datafile);    
female_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/window_1kb/PLUT_FEMALE_PX8.sorted.bed';
F_1 = load_data_from_BED(female_datafile);
female_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/window_1kb/PLUT_FEMALE_PX9_1.sorted.bed';
F_2 = load_data_from_BED(female_datafile);
female_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/window_1kb/PLUT_FEMALE_PX9_2.sorted.bed';
F_3 = load_data_from_BED(female_datafile);
female_datafile = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/window_1kb/PLUT_FEMALE_PX10.sorted.bed';
F_4 = load_data_from_BED(female_datafile);


%%
% Considering all the female samples
% F_depth = (F_1.depth + F_2.depth + F_3.depth + F_4.depth)/4;
% F_baseCovered = (F_1.bases_covered_in_feature + F_2.bases_covered_in_feature + F_3.bases_covered_in_feature + F_4.bases_covered_in_feature) / 4;
% F_breadth = (F_1.breadth + F_2.breadth + F_3.breadth + F_4.breadth)/4;

% Considering only the good female samples
F_depth = round((F_2.depth + F_3.depth)/2); % round it as it should be an int
F_baseCovered = round((F_2.bases_covered_in_feature + F_3.bases_covered_in_feature) / 2);% round it as it should be an int
F_breadth = (F_2.breadth + F_3.breadth)/2;

F_contigLength = F_1.contig_length;

FEMALE_DATA = [F_depth F_baseCovered F_contigLength F_breadth];    
MALE_DATA = [M.depth M.bases_covered_in_feature M.contig_length M.breadth];
feature_list = M.features; % Feature list is in the same order for both MALE and FEMALE
col_depth = 1;
col_baseCovered = 2;
col_featureLength = 3;
col_breadth = 4;

%% Find all zero males
zero_depth_idx_male = find(MALE_DATA(:,col_depth) == 0);
% Find all zero females
zero_depth_idx_female = find(FEMALE_DATA(:,col_depth) == 0);
% find zero male and female
zero_depth_idx_common = intersect(zero_depth_idx_male, zero_depth_idx_female);

% Remove features that are zero in both
filtered_feature_list = feature_list;
filtered_feature_list(zero_depth_idx_common) = [];
Filtered_MALE_DATA = MALE_DATA;
Filtered_MALE_DATA(zero_depth_idx_common,:) = [];
Filtered_FEMALE_DATA = FEMALE_DATA;
Filtered_FEMALE_DATA(zero_depth_idx_common,:) = [];

%%
W_male_idx = find(Filtered_MALE_DATA(:,col_depth) == 0); 
W_contig_names = filtered_feature_list(W_male_idx);
W_contig_length = Filtered_MALE_DATA(W_male_idx,col_featureLength);
W_corresponding_Female_depth = Filtered_FEMALE_DATA(W_male_idx,col_depth);

%% Find opposite of above - Present in male, absent in female
opp_W_male_idx = find(Filtered_FEMALE_DATA(:,col_depth) == 0); 
opp_W_contig_names = filtered_feature_list(opp_W_male_idx);
opp_W_contig_length = Filtered_FEMALE_DATA(opp_W_male_idx,col_featureLength);
opp_W_corresponding_male_depth = Filtered_MALE_DATA(opp_W_male_idx,col_depth);
%% Look at the distributions
% Remove zero male and zero females as we know them from above W, opp-W
% analyses
zero_male_female_indices = [W_male_idx; opp_W_male_idx];
Filtered_MALE_DATA(zero_male_female_indices,:) = [];
Filtered_FEMALE_DATA(zero_male_female_indices,:) = [];
filtered_feature_list(zero_male_female_indices,:) = [];

male_depth = log10(Filtered_MALE_DATA(:,col_depth));
female_depth = log10(Filtered_FEMALE_DATA(:,col_depth));

figure;
bin_size = 0.02;
h1 = histogram(male_depth,'Normalization','pdf');
h1.BinWidth = bin_size;
hold on
h2 = histogram(female_depth,'Normalization','pdf');
h2.BinWidth = bin_size;


%% up to here - same as Plut_DistStudy_BEDFiles.m %%
lower_threshold = 1.8;
upper_threshold = 3.5;
f_idx_1 = find(female_depth>=lower_threshold);
f_idx_2 = find(female_depth<=upper_threshold);
f_idx = intersect(f_idx_1,f_idx_2);
figure;
histogram(female_depth(f_idx))
figure;
histogram(male_depth(f_idx))
figure;
bin_size = 0.02;
h1 = histogram(female_depth(f_idx),'Normalization','pdf')
h1.BinWidth = bin_size;
hold on
h2 = histogram(male_depth(f_idx),'Normalization','pdf')
h2.BinWidth = bin_size;

%%
% Fit tri-modal mixture model
[params_est,ressquared] = Plut_fit_Model_optimize(female_depth)

figure;
hold on
histogram(female_depth,'Normalization','probability')
height = 0.07;
plot([params_est(4),params_est(4)],[0,height],'Color','yellow','LineWidth',2)
plot([params_est(4) - params_est(5),params_est(4) - params_est(5)],[0,height],'LineStyle','--','Color','yellow')
plot([params_est(4) + params_est(5),params_est(4) + params_est(5)],[0,height],'LineStyle','--','Color','yellow')

plot([params_est(6),params_est(6)],[0,height],'Color','red','LineWidth',2 )
plot([params_est(6) - params_est(7),params_est(6) - params_est(7)],[0,height],'LineStyle','--' ,'Color','red'  )
plot([params_est(6) + params_est(7),params_est(6) + params_est(7)],[0,height],'LineStyle','--' ,'Color','red'  )

hold off

%% Filter out the features of interest
mean_seek = params_est(4);
sig_seek = params_est(5);
upp_limit = mean_seek + 0.25 * sig_seek; % Note - half of the std-dev
low_limit = mean_seek - 0.25 * sig_seek;
% Find all features that lie between upper and lower limit

selected_idx_within_low_upp_limit = [];
for i=1:numel(female_depth)
    if (female_depth(i) >= low_limit) & (female_depth(i) <=  upp_limit)
        selected_idx_within_low_upp_limit(end+1) = i;
    end
end


%% Plot selected female values from the window along with corresponding male
% values
figure;
histogram(female_depth(selected_idx_within_low_upp_limit));
hold on
histogram(male_depth(selected_idx_within_low_upp_limit));
% Thought line 1 -
% The above figure suggests that we have male values coming from the entire
% range. Need to device a strategy for considering only those windows that
% exhibit the desired pattern,i.e., windows within the selected range for
% female, but should have the corresponding male value exceeding the upper
% limit in female distribution
idx_windowsOutsideFemaleRangeInMale = [];
for i=1:numel(selected_idx_within_low_upp_limit)
    if male_depth(selected_idx_within_low_upp_limit(i)) > upp_limit
        idx_windowsOutsideFemaleRangeInMale(end+1) = selected_idx_within_low_upp_limit(i);
    end
end
figure;
histogram(female_depth(idx_windowsOutsideFemaleRangeInMale));
hold on
histogram(male_depth(idx_windowsOutsideFemaleRangeInMale));

% Thought line 2 - 
% On the other thought, it makes it logical that the male values are coming
% from all spectrum for the selected female values, becuase we expect only
% certain regions in the male regions to be repetitive and the remaining
% ones would be of lower coverage.

% TO-DO something about it....
    

%% Select those candidates that have window size totalling more than 20KB.
idx_to_select = selected_idx_within_low_upp_limit;
%idx_to_select = idx_windowsOutsideFemaleRangeInMale;

selected_features = filtered_feature_list(idx_to_select);
selected_depth = female_depth(idx_to_select);

uniq_selected_features = unique(selected_features);
uniq_selected_features_windowCount = zeros(numel(uniq_selected_features),1);
for i=1:numel(uniq_selected_features)
    uniq_selected_features_windowCount(i) = sum(ismember(selected_features,uniq_selected_features(i)));
end

figure;histogram(log10(uniq_selected_features_windowCount ))

% Each wondow is 1 KB, so we need about several windows 
min_windows_desired = 20;
uniq_selected_features_windowCount_minWin_idx = find(uniq_selected_features_windowCount >= min_windows_desired);

desired_features = uniq_selected_features(uniq_selected_features_windowCount_minWin_idx);
desired_win_counts = uniq_selected_features_windowCount(uniq_selected_features_windowCount_minWin_idx);

figure; histogram(desired_win_counts)

[sorted_desired_win_counts, sorted_desired_win_counts_idx] = sort(desired_win_counts,'descend');

OutputListFile = 'Plut_windowBasedResults.txt'
fileID = fopen(OutputListFile,'w');
fprintf(fileID,'## Plutella window based identifications - 1KB window size for analysis \n');
for i=1:numel(sorted_desired_win_counts)
    fprintf(fileID,'%s \t %d \n',desired_features{sorted_desired_win_counts_idx(i)},sorted_desired_win_counts(i));
    
end
fclose(fileID);

%% Aggregated view = Figures for selected male and female features 
idx = ismember(filtered_feature_list,desired_features);
selected_female_depth = female_depth(idx);
selected_male_depth = male_depth(idx);
figure;
histogram(selected_female_depth,'Normalization','pdf')
hold on
histogram(selected_male_depth,'Normalization','pdf')

%% Create batch file
contig_list = {numel(sorted_desired_win_counts)};
for i=1:numel(sorted_desired_win_counts)
    contig_list{i} = desired_features{sorted_desired_win_counts_idx(i)};
end
Plut_createIGVBatchFile;
