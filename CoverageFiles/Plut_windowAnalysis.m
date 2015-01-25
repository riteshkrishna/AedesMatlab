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


%% Prepare the data-structure
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

%% Clean data by removing zero value entries. Take a note of presence/absence in male/female 
% Find all zero males
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

% Find Present in female, absent in male
W_male_idx = find(Filtered_MALE_DATA(:,col_depth) == 0); 
W_contig_names = filtered_feature_list(W_male_idx);
W_contig_length = Filtered_MALE_DATA(W_male_idx,col_featureLength);
W_corresponding_Female_depth = Filtered_FEMALE_DATA(W_male_idx,col_depth);

% Find opposite of above - Present in male, absent in female
opp_W_male_idx = find(Filtered_FEMALE_DATA(:,col_depth) == 0); 
opp_W_contig_names = filtered_feature_list(opp_W_male_idx);
opp_W_contig_length = Filtered_FEMALE_DATA(opp_W_male_idx,col_featureLength);
opp_W_corresponding_male_depth = Filtered_MALE_DATA(opp_W_male_idx,col_depth);

%% Prepare data for downstream analysis
% Remove zero male and zero females as we know them from above W, opp-W
% analyses
zero_male_female_indices = [W_male_idx; opp_W_male_idx];
Filtered_MALE_DATA(zero_male_female_indices,:) = [];
Filtered_FEMALE_DATA(zero_male_female_indices,:) = [];
filtered_feature_list(zero_male_female_indices,:) = [];

% Finally, the data that needs to be analysed
male_depth = log10(Filtered_MALE_DATA(:,col_depth));
female_depth = log10(Filtered_FEMALE_DATA(:,col_depth));

figure;
bin_size = 0.02;
h1 = histogram(male_depth,'Normalization','pdf');
h1.BinWidth = bin_size;
hold on
h2 = histogram(female_depth,'Normalization','pdf');
h2.BinWidth = bin_size;
legend('Male','Female');
title('Depth computed log scale for 1KB wide windows')


% %% Plot a filtered histogram by setting thresholds
% lower_threshold = 1.8;
% upper_threshold = 3.5;
% f_idx_1 = find(female_depth>=lower_threshold);
% f_idx_2 = find(female_depth<=upper_threshold);
% f_idx = intersect(f_idx_1,f_idx_2);
% %figure; histogram(female_depth(f_idx))
% %figure; histogram(male_depth(f_idx))
% figure;
% bin_size = 0.02;
% h1 = histogram(female_depth(f_idx),'Normalization','pdf')
% h1.BinWidth = bin_size;
% hold on
% h2 = histogram(male_depth(f_idx),'Normalization','pdf')
% h2.BinWidth = bin_size;

% %% Fit a tri-modal mixture model to the data
% [params_est,ressquared] = Plut_fit_Model_optimize(female_depth)
% 
% figure;
% hold on
% histogram(female_depth,'Normalization','probability')
% height = 0.07;
% plot([params_est(4),params_est(4)],[0,height],'Color','yellow','LineWidth',2)
% plot([params_est(4) - params_est(5),params_est(4) - params_est(5)],[0,height],'LineStyle','--','Color','yellow')
% plot([params_est(4) + params_est(5),params_est(4) + params_est(5)],[0,height],'LineStyle','--','Color','yellow')
% 
% plot([params_est(6),params_est(6)],[0,height],'Color','red','LineWidth',2 )
% plot([params_est(6) - params_est(7),params_est(6) - params_est(7)],[0,height],'LineStyle','--' ,'Color','red'  )
% plot([params_est(6) + params_est(7),params_est(6) + params_est(7)],[0,height],'LineStyle','--' ,'Color','red'  )
% 
% hold off
% 
% %% Filter out the features of interest
% sig_level = 1; % Interest boundary - move by this threshold
% mean_seek = params_est(4);
% sig_seek = params_est(5);
% upp_limit = mean_seek + sig_level * sig_seek; 
% low_limit = mean_seek - sig_level * sig_seek;
% 
% % Find all features that lie between upper and lower limit
% selected_idx_within_low_upp_limit = [];
% for i=1:numel(female_depth)
%     if (female_depth(i) >= low_limit) & (female_depth(i) <=  upp_limit)
%         selected_idx_within_low_upp_limit(end+1) = i;
%     end
% end
% 
% 
% % Plot selected female values from the window along with corresponding male values
% figure;
% histogram(female_depth(selected_idx_within_low_upp_limit));
% hold on
% histogram(male_depth(selected_idx_within_low_upp_limit));
% 
% %% Strategies to analyse the selected windows
% %% Thought line 1 -
% % The above figure suggests that we have male values coming from the entire
% % range. Need to device a strategy for considering only those windows that
% % exhibit the desired pattern,i.e., windows within the selected range for
% % female, but should have the corresponding male value exceeding the upper
% % limit in female distribution
% idx_windowsOutsideFemaleRangeInMale = [];
% for i=1:numel(selected_idx_within_low_upp_limit)
%     if male_depth(selected_idx_within_low_upp_limit(i)) > upp_limit
%         idx_windowsOutsideFemaleRangeInMale(end+1) = selected_idx_within_low_upp_limit(i);
%     end
% end
% figure;
% histogram(female_depth(idx_windowsOutsideFemaleRangeInMale));
% hold on
% histogram(male_depth(idx_windowsOutsideFemaleRangeInMale));
% 
% %% Thought line 2 - 
% % On the other thought, it makes it logical that the male values are coming
% % from all spectrum for the selected female values, becuase we expect only
% % certain regions in the male regions to be repetitive and the remaining
% % ones would be of lower coverage.
% 
% % TO-DO something about it....
%     
% %% Thought line 3 - 
% % every window selected by above criteria belong to a contig. Some contigs
% % may have multiple windows in the region of interest. Count how many
% % contigs these windows belong to. Also determine that what is the total 
% % number of windows (in the selected region) for each of those contigs.
% idx_to_select = selected_idx_within_low_upp_limit;
% selected_features = filtered_feature_list(idx_to_select);
% selected_depth = female_depth(idx_to_select);
% selected_lengths = Filtered_FEMALE_DATA(idx_to_select,3);
% 
% uniq_selected_features = unique(selected_features);
% uniq_selected_features_windowCount = zeros(numel(uniq_selected_features),1);
% for i=1:numel(uniq_selected_features)
%     uniq_selected_features_windowCount(i) = sum(ismember(selected_features,uniq_selected_features(i)));
% end
% figure;histogram((uniq_selected_features_windowCount ))
% 
% %% Option 1 - Candidates based on size
% % Select those candidates that have window size totalling more than 20KB.
% % Each wondow is 1 KB, so we need about several windows 
% min_windows_desired = 20;
% uniq_selected_features_windowCount_minWin_idx = find(uniq_selected_features_windowCount >= min_windows_desired);
% desired_features = uniq_selected_features(uniq_selected_features_windowCount_minWin_idx);
% desired_win_counts = uniq_selected_features_windowCount(uniq_selected_features_windowCount_minWin_idx);
% figure; histogram(desired_win_counts)
% 
% % Aggregated view = Figures for selected male and female features 
% idx = ismember(filtered_feature_list,desired_features);
% selected_female_depth = female_depth(idx);
% selected_male_depth = male_depth(idx);
% figure;
% histogram(selected_female_depth,'Normalization','pdf')
% hold on
% histogram(selected_male_depth,'Normalization','pdf')
% 
% data_to_sort = desired_win_counts;
% data_to_sort_names = desired_features;
% OutputListFile = 'TEST_Plut_windowBasedResults.txt';
% file_header = '## Plutella window based identifications - 1KB window size for analysis \n';
% Plut_file_write(OutputListFile,data_to_sort_names,data_to_sort,file_header);
% 
% 
% %% Option - 2
% % Selection of candidates based on size-proportion
% proportion_covered_by_selected_windows = zeros(numel(uniq_selected_features),1);
% total_contig_length  = zeros(numel(uniq_selected_features),1);
% for i=1:numel(uniq_selected_features)
%     idx_for_length = ismember(selected_features,uniq_selected_features(i));
%     length_covered_by_selected_windows = sum(selected_lengths(idx_for_length)); % Since each wondow is not 1KB full, we need to get to raw data
%     idx_in_bed_file = ismember(filtered_feature_list,uniq_selected_features(i));
%     total_contig_length(i) = sum(Filtered_FEMALE_DATA(idx_in_bed_file,3));
%     proportion_covered_by_selected_windows(i) = length_covered_by_selected_windows / total_contig_length(i);
% end
% figure; histogram(proportion_covered_by_selected_windows);
% 
% % Filter on proportion threshold
% proportion_threshold = 0.5;
% idx_selected_proportions = find(proportion_covered_by_selected_windows >= proportion_threshold);
% uniq_selected_features_filtered = uniq_selected_features(idx_selected_proportions);
% proportion_covered_by_selected_windows_filtered = proportion_covered_by_selected_windows(idx_selected_proportions);
% total_contig_length_filtered = total_contig_length(idx_selected_proportions);
% 
% % Filter on total contig length threshold
% contig_length_threshold = 20000;
% idx_filtered_proportion_contiglength = find(total_contig_length_filtered >= contig_length_threshold);
% uniq_selected_features_filtered_proportion_contiglength = uniq_selected_features_filtered(idx_filtered_proportion_contiglength);
% total_contig_length_filtered_proportion_contiglength = total_contig_length_filtered(idx_filtered_proportion_contiglength);
% 
% % Write to a file
% %data_to_sort = proportion_covered_by_selected_windows_filtered;
% %data_to_sort_names = uniq_selected_features_filtered;
% data_to_sort = total_contig_length_filtered_proportion_contiglength;
% data_to_sort_names = uniq_selected_features_filtered_proportion_contiglength;
% 
% OutputListFile = 'TEST_Plut_windowBasedResults_proportion.txt'
% file_header = ['## Plutella - 1KB window,  proportion threshold >= 0.5  , contig length >= 20kb \n' ...
%                 '## Contig names \t proportion of contig in the region of interest \n'];
% Plut_file_write(OutputListFile,data_to_sort_names,data_to_sort,file_header);
% 
% % NOTE - This selects very small profiles mostly - useless like that
% 
% %% Option 3 
% % Look for those windows that have at least double coverage in male
% idx_to_select = selected_idx_within_low_upp_limit;
% selected_features = filtered_feature_list(idx_to_select);
% %selected_depth_female = female_depth(idx_to_select);
% selected_lengths = Filtered_FEMALE_DATA(idx_to_select,3);
% %selected_depth_male = male_depth(idx_to_select);
% 
% selected_depth_male = Filtered_MALE_DATA(idx_to_select,col_depth);
% selected_depth_female = Filtered_FEMALE_DATA(idx_to_select,col_depth);
% 
% ratio_male_female = selected_depth_male ./ selected_depth_female;
% ratio_idx = find (ratio_male_female >= 2);
% ratio_selected_features = selected_features(ratio_idx);
% ratio_selected_lengths = selected_lengths(ratio_idx);
% uniq_ratio_selected_features = unique(ratio_selected_features);
% 
% proportion_covered_by_selected_windows = zeros(numel(uniq_ratio_selected_features),1);
% total_contig_length  = zeros(numel(uniq_ratio_selected_features),1);
% for i=1:numel(uniq_ratio_selected_features)
%     idx_for_length = ismember(ratio_selected_features,uniq_ratio_selected_features(i));
%     length_covered_by_selected_windows = sum(ratio_selected_lengths(idx_for_length)); % Since each wondow is not 1KB full, we need to get to raw data
%     idx_in_bed_file = ismember(filtered_feature_list,uniq_ratio_selected_features(i));
%     total_contig_length(i) = sum(Filtered_FEMALE_DATA(idx_in_bed_file,3));
%     proportion_covered_by_selected_windows(i) = length_covered_by_selected_windows / total_contig_length(i);
% end
% figure; histogram(proportion_covered_by_selected_windows);
% 
% data_to_sort = proportion_covered_by_selected_windows;
% data_to_sort_names = uniq_ratio_selected_features;
% OutputListFile = 'TEST_Plut_windowBasedResults_proportion.txt'
% file_header = ['## Plutella - 1KB window,  Male/Female ratio > 2 windows \n' ...
%                 '## Contig names \t proportion of contig in the region of interest \n'];
% %Plut_file_write(OutputListFile,data_to_sort_names,data_to_sort,file_header);
% [sorted_data, sorted_idx] = sort(data_to_sort,'descend');
% fileID = fopen(OutputListFile,'w');
% fprintf(fileID,file_header);
% for i=1:numel(sorted_data)
%     fprintf(fileID,'%s \t %d \t %d \n',data_to_sort_names{sorted_idx(i)},sorted_data(i), total_contig_length(sorted_idx(i)));
% end
% fclose(fileID);
% % NOTE - Needs further refinements. Need to know the regions for long
% % contigs. Needs to create a global measure for the contigs.
% 
% %% Create batch file for IGV
% 
% [sorted_data, sorted_idx] = sort(data_to_sort,'descend');
% contig_list = {numel(sorted_data)};
% for i=1:numel(sorted_data)
%     contig_list{i} = data_to_sort_names{sorted_idx(i)};
% end
% batchFile = 'igv_Plut_1kbwindow_proportion_27112014.txt';
% Plut_createIGVBatchFile(batchFile,contig_list)
% 
