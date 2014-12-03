%% The analysis part of Plut_windowAnalysis.m moved here..
function Plut_windowAnalysis_part_2()
    
    Plut_windowAnalysis;
    
    [params] = fit_model(female_depth);
    [selected_idx_within_low_upp_limit,upp_limit,low_limit] = filter_region_of_interest(params,female_depth,male_depth);
    
    % Save for use in individual functions to avoid running the whole code
    % again
    save_file_name = 'plut_analysis.mat';
    save(save_file_name);
    
    
    %thought_1(selected_idx_within_low_upp_limit,upp_limit, male_depth,female_depth)
    
    %min_windows_desired = 20;
    %thought_3_option_1(selected_idx_within_low_upp_limit,female_depth,male_depth,...
    %            filtered_feature_list,Filtered_FEMALE_DATA,min_windows_desired);
     
    %thought_3_option_2(selected_idx_within_low_upp_limit,female_depth,male_depth,...
    %            filtered_feature_list,Filtered_FEMALE_DATA)
    
    %thought_3_option_3(selected_idx_within_low_upp_limit,female_depth,male_depth,...
    %            filtered_feature_list,Filtered_FEMALE_DATA,Filtered_MALE_DATA)
    
    
    thought_3_option_4(selected_idx_within_low_upp_limit,filtered_feature_list,Filtered_FEMALE_DATA,Filtered_MALE_DATA)
end

%% Fit a tri-modal mixture model to the data
function [params_est] = fit_model(female_depth)
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
end

%% Filter out the features of interest
function [selected_idx_within_low_upp_limit,upp_limit,low_limit] = filter_region_of_interest(params_est,female_depth,male_depth)
    sig_level = 1; % Interest boundary - move by this threshold
    mean_seek = params_est(4);
    sig_seek = params_est(5);
    upp_limit = mean_seek + sig_level * sig_seek; 
    low_limit = mean_seek - sig_level * sig_seek;

    % Find all features that lie between upper and lower limit
    selected_idx_within_low_upp_limit = [];
    for i=1:numel(female_depth)
        if (female_depth(i) >= low_limit) & (female_depth(i) <=  upp_limit)
            selected_idx_within_low_upp_limit(end+1) = i;
        end
    end
    
    % Plot selected female values from the window along with corresponding male values
    figure;
    histogram(female_depth(selected_idx_within_low_upp_limit));
    hold on
    histogram(male_depth(selected_idx_within_low_upp_limit));

end

%% Strategies to analyse the selected windows
%% Thought line 1 -
% The above figure suggests that we have male values coming from the entire
% range. Need to device a strategy for considering only those windows that
% exhibit the desired pattern,i.e., windows within the selected range for
% female, but should have the corresponding male value exceeding the upper
% limit in female distribution
function [idx_windowsOutsideFemaleRangeInMale] = thought_1(selected_idx_within_low_upp_limit,upp_limit, male_depth,female_depth)
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
end

%% Thought line 2 - 
% On the other thought, it makes it logical that the male values are coming
% from all spectrum for the selected female values, becuase we expect only
% certain regions in the male regions to be repetitive and the remaining
% ones would be of lower coverage.
% TO-DO something about it....
function thought_2()

end
    
%% Thought line 3 - 
% every window selected by above criteria belong to a contig. Some contigs
% may have multiple windows in the region of interest. Count how many
% contigs these windows belong to. Also determine that what is the total 
% number of windows (in the selected region) for each of those contigs.

function thought_3_option_1(selected_idx_within_low_upp_limit,female_depth,male_depth,...
                filtered_feature_list,Filtered_FEMALE_DATA,min_windows_desired)
    
    idx_to_select = selected_idx_within_low_upp_limit;
    selected_features = filtered_feature_list(idx_to_select);
    selected_depth = female_depth(idx_to_select);
    selected_lengths = Filtered_FEMALE_DATA(idx_to_select,3);

    uniq_selected_features = unique(selected_features);
    uniq_selected_features_windowCount = zeros(numel(uniq_selected_features),1);
    for i=1:numel(uniq_selected_features)
        uniq_selected_features_windowCount(i) = sum(ismember(selected_features,uniq_selected_features(i)));
    end
    figure;histogram((uniq_selected_features_windowCount ))
    
    uniq_selected_features_windowCount_minWin_idx = find(uniq_selected_features_windowCount >= min_windows_desired);
    desired_features = uniq_selected_features(uniq_selected_features_windowCount_minWin_idx);
    desired_win_counts = uniq_selected_features_windowCount(uniq_selected_features_windowCount_minWin_idx);
    figure; histogram(desired_win_counts)

    % Aggregated view = Figures for selected male and female features 
    idx = ismember(filtered_feature_list,desired_features);
    selected_female_depth = female_depth(idx);
    selected_male_depth = male_depth(idx);
    figure;
    histogram(selected_female_depth,'Normalization','pdf')
    hold on
    histogram(selected_male_depth,'Normalization','pdf')

    data_to_sort = desired_win_counts;
    data_to_sort_names = desired_features;
    OutputListFile = 'TEST_Plut_windowBasedResults_thought_3_option_1.txt';
    file_header = '## Plutella window based identifications - 1KB window size for analysis \n';
    Plut_file_write(OutputListFile,data_to_sort_names,data_to_sort,file_header);
    
end

%% Option - 2
% Selection of candidates based on size-proportion
function thought_3_option_2(selected_idx_within_low_upp_limit,female_depth,male_depth,...
                filtered_feature_list,Filtered_FEMALE_DATA)
    
    idx_to_select = selected_idx_within_low_upp_limit;
    selected_features = filtered_feature_list(idx_to_select);
    selected_depth = female_depth(idx_to_select);
    selected_lengths = Filtered_FEMALE_DATA(idx_to_select,3);

    uniq_selected_features = unique(selected_features);
    uniq_selected_features_windowCount = zeros(numel(uniq_selected_features),1);
    for i=1:numel(uniq_selected_features)
        uniq_selected_features_windowCount(i) = sum(ismember(selected_features,uniq_selected_features(i)));
    end
    figure;histogram((uniq_selected_features_windowCount ))
    
    proportion_covered_by_selected_windows = zeros(numel(uniq_selected_features),1);
    
    total_contig_length  = zeros(numel(uniq_selected_features),1);
    for i=1:numel(uniq_selected_features)
        idx_for_length = ismember(selected_features,uniq_selected_features(i));
        length_covered_by_selected_windows = sum(selected_lengths(idx_for_length)); % Since each wondow is not 1KB full, we need to get to raw data
        idx_in_bed_file = ismember(filtered_feature_list,uniq_selected_features(i));
        total_contig_length(i) = sum(Filtered_FEMALE_DATA(idx_in_bed_file,3));
        proportion_covered_by_selected_windows(i) = length_covered_by_selected_windows / total_contig_length(i);
    end
    figure; histogram(proportion_covered_by_selected_windows);

    % Filter on proportion threshold
    proportion_threshold = 0.5;
    idx_selected_proportions = find(proportion_covered_by_selected_windows >= proportion_threshold);
    uniq_selected_features_filtered = uniq_selected_features(idx_selected_proportions);
    proportion_covered_by_selected_windows_filtered = proportion_covered_by_selected_windows(idx_selected_proportions);
    total_contig_length_filtered = total_contig_length(idx_selected_proportions);

    % Filter on total contig length threshold
    contig_length_threshold = 20000;
    idx_filtered_proportion_contiglength = find(total_contig_length_filtered >= contig_length_threshold);
    uniq_selected_features_filtered_proportion_contiglength = uniq_selected_features_filtered(idx_filtered_proportion_contiglength);
    total_contig_length_filtered_proportion_contiglength = total_contig_length_filtered(idx_filtered_proportion_contiglength);

    % Write to a file
    %data_to_sort = proportion_covered_by_selected_windows_filtered;
    %data_to_sort_names = uniq_selected_features_filtered;
    data_to_sort = total_contig_length_filtered_proportion_contiglength;
    data_to_sort_names = uniq_selected_features_filtered_proportion_contiglength;

    OutputListFile = 'TEST_Plut_windowBasedResults_thought_3_option_2.txt'
    file_header = ['## Plutella - 1KB window,  proportion threshold >= 0.5  , contig length >= 20kb \n' ...
                '## Contig names \t proportion of contig in the region of interest \n'];
    Plut_file_write(OutputListFile,data_to_sort_names,data_to_sort,file_header);

end

%% Option 3 
% Look for those windows that have at least double coverage in male
function thought_3_option_3(selected_idx_within_low_upp_limit,female_depth,male_depth,...
                filtered_feature_list,Filtered_FEMALE_DATA,Filtered_MALE_DATA)
            
    idx_to_select = selected_idx_within_low_upp_limit;
    selected_features = filtered_feature_list(idx_to_select);
    %selected_depth_female = female_depth(idx_to_select);
    selected_lengths = Filtered_FEMALE_DATA(idx_to_select,3);
    %selected_depth_male = male_depth(idx_to_select);

    selected_depth_male = Filtered_MALE_DATA(idx_to_select,1);
    selected_depth_female = Filtered_FEMALE_DATA(idx_to_select,1);

    ratio_male_female = selected_depth_male ./ selected_depth_female;
    ratio_idx = find (ratio_male_female >= 2);
    ratio_selected_features = selected_features(ratio_idx);
    ratio_selected_lengths = selected_lengths(ratio_idx);
    uniq_ratio_selected_features = unique(ratio_selected_features);

    proportion_covered_by_selected_windows = zeros(numel(uniq_ratio_selected_features),1);
    total_contig_length  = zeros(numel(uniq_ratio_selected_features),1);
    for i=1:numel(uniq_ratio_selected_features)
        idx_for_length = ismember(ratio_selected_features,uniq_ratio_selected_features(i));
        length_covered_by_selected_windows = sum(ratio_selected_lengths(idx_for_length)); % Since each wondow is not 1KB full, we need to get to raw data
        idx_in_bed_file = ismember(filtered_feature_list,uniq_ratio_selected_features(i));
        total_contig_length(i) = sum(Filtered_FEMALE_DATA(idx_in_bed_file,3));
        proportion_covered_by_selected_windows(i) = length_covered_by_selected_windows / total_contig_length(i);
    end
    figure; histogram(proportion_covered_by_selected_windows);

    data_to_sort = proportion_covered_by_selected_windows;
    data_to_sort_names = uniq_ratio_selected_features;
    OutputListFile = 'TEST_Plut_windowBasedResults_thought_3_option_3.txt'
    file_header = ['## Plutella - 1KB window,  Male/Female ratio > 2 windows \n' ...
                '## Contig names \t proportion of contig in the region of interest \n'];
    %Plut_file_write(OutputListFile,data_to_sort_names,data_to_sort,file_header);
    [sorted_data, sorted_idx] = sort(data_to_sort,'descend');
    fileID = fopen(OutputListFile,'w');
    fprintf(fileID,file_header);
    for i=1:numel(sorted_data)
        fprintf(fileID,'%s \t %d \t %d \n',data_to_sort_names{sorted_idx(i)},sorted_data(i), total_contig_length(sorted_idx(i)));
    end
    fclose(fileID);
    % NOTE - Needs further refinements. Need to know the regions for long
    % contigs. Needs to create a global measure for the contigs.
end

%% Option 4 - Approach based on counting neighbouring windows
function thought_3_option_4(selected_idx_within_low_upp_limit,filtered_feature_list,Filtered_FEMALE_DATA,Filtered_MALE_DATA)
    
    save_file_name = 'plut_analysis.mat';
    load(save_file_name);
    
    col_depth = 1;
    %col_baseCovered = 2;
    %col_featureLength = 3;
    %col_breadth = 4;

    idx_to_select = selected_idx_within_low_upp_limit;
    selected_features = filtered_feature_list(idx_to_select);
    uniq_selected_features = unique(selected_features);
    
    % Concentrate just on window counting rather than proportion etc.
    % Create a data-structure from female data we are processing
    cell_feature_win = cell(numel(uniq_selected_features),1);
    for i=1:numel(uniq_selected_features)
        idx_in_window = idx_to_select(ismember(selected_features,uniq_selected_features(i))); 
        idx_bed_file = find(ismember(filtered_feature_list,uniq_selected_features(i))==1);
        cell_feature_win{i,1} = uniq_selected_features(i); % name of contig
        cell_feature_win{i,2} = idx_in_window; % index in selected region
        cell_feature_win{i,3} = idx_bed_file; % index in bed file
        cell_feature_win{i,4} = ismember( cell_feature_win{i,3},cell_feature_win{i,2}); % get a binary matrix showing proximity of windows on genome
        cell_feature_win{i,5} = numel(cell_feature_win{i,2}); % total number of windows in region     
        cell_feature_win{i,6} = sum(cell_feature_win{i,4}) / numel(cell_feature_win{i,4}); % proportion of windows in the region wrt total windows
    end
    
    % Get the above information from male data as well
    male_cell_feature_win = cell(size(cell_feature_win,1),1);
    for i=1:size(cell_feature_win,1)
        male_cell_feature_win{i,1} = cell_feature_win{i,1}; % name
        male_cell_feature_win{i,2} = cell_feature_win{i,3}; % idx in bed file
        
        female_depth_for_bed_idx = Filtered_FEMALE_DATA([male_cell_feature_win{i,2}],col_depth); % non-log female depth
        male_depth_for_bed_idx = Filtered_MALE_DATA([male_cell_feature_win{i,2}],col_depth);     % non-log male depth   
        ratio_male_female_depth = male_depth_for_bed_idx ./ female_depth_for_bed_idx;
        
        male_cell_feature_win{i,3} = ratio_male_female_depth; % male-female depth ratio for those bed indices
        
        % Late addition
        male_cell_feature_win{i,4} = male_depth_for_bed_idx; % male depth
        male_cell_feature_win{i,5} = female_depth_for_bed_idx; % female depth
    end
    
    % thoughts -
    % cell_feature_win{i,4} gives binary rep. of bed indices in that are in selected region
    % cell_feature_win{i,6} gives proportion of windows in selected region
    % male_cell_feature_win{i,3} gives male-female ratio of bed indices
    
    % Select on size
    % min numner of windows that must fall within the region
    win_threshold = 5;
    idx_gt_win_threshold = find (([cell_feature_win{:,5} ] > win_threshold) == 1);
    
    % then - select on male-female ratio
    ratio_threshold = 2;
    idx_gt_win_ratio_threshold = [];
    for i=1:numel(idx_gt_win_threshold)
        idx = idx_gt_win_threshold(i);
        mean_ratio = mean(male_cell_feature_win{idx,3});
        std_ratio = std(male_cell_feature_win{idx,3});
        if (mean_ratio + std_ratio) >= ratio_threshold
            idx_gt_win_ratio_threshold(end+1,1) = idx;
        end
    end
    
    cells_filtered_on_win_ratio_threshold = cell_feature_win(idx_gt_win_ratio_threshold,:);
    male_cells_filtered_on_win_ratio_threshold = male_cell_feature_win(idx_gt_win_ratio_threshold,:);
    
    % check how they look
    wins = zeros(size(cells_filtered_on_win_ratio_threshold,1),1);
    for i=1:size(cells_filtered_on_win_ratio_threshold,1)
        wins(i,1) = numel([cells_filtered_on_win_ratio_threshold{i,3}]);
    end
    scatter(wins,[cells_filtered_on_win_ratio_threshold{:,6}])
    figure; histogram(wins)
    figure; histogram([cells_filtered_on_win_ratio_threshold{:,6}])
    
    % Full data
    male_depth = log10(Filtered_MALE_DATA(:,col_depth));
    female_depth = log10(Filtered_FEMALE_DATA(:,col_depth));
    % Selected data
    male_selected_depth = [];
    female_selected_depth = [];
    for i=1:size(male_cells_filtered_on_win_ratio_threshold,1)
        male_selected_depth = [male_selected_depth;male_cells_filtered_on_win_ratio_threshold{i,4}];
        female_selected_depth = [female_selected_depth;male_cells_filtered_on_win_ratio_threshold{i,5}];
    end
    
    % plot selected data
    figure; 
    histogram(log10(male_selected_depth),'Normalization','pdf');
    hold on
    histogram(log10(female_selected_depth),'Normalization','pdf');
    hold off
    % plot full and selected data together
    figure;
    bin_size = 0.02;
    h1 = histogram(male_depth,'Normalization','pdf');
    h1.BinWidth = bin_size;
    hold on
    h2 = histogram(female_depth,'Normalization','pdf');
    h2.BinWidth = bin_size;
    h3 = histogram(log10(male_selected_depth),'Normalization','pdf');
    h3.BinWidth = bin_size;
    h4 = histogram(log10(female_selected_depth),'Normalization','pdf');
    h4.BinWidth = bin_size;
    % %%%%% Finished plotting %%%%%%
    
    
    % Determine total contig-lengths
    contig_lengths = zeros(size(cells_filtered_on_win_ratio_threshold,1),1);
    for i=1:size(cells_filtered_on_win_ratio_threshold,1)
        idx_in_bed_file = [cells_filtered_on_win_ratio_threshold{i,3}];
        contig_lengths(i) = sum(Filtered_FEMALE_DATA(idx_in_bed_file,3));
    end
    
    % Write to file - ordered by contig length
    data_to_sort = contig_lengths;
    [sorted_data, sorted_idx] = sort(data_to_sort,'descend');
    
    OutputListFile = 'TEST_Plut_windowBasedResults_thought_3_option_4.txt'
    file_header = ['## Plutella - 1KB window,  Male/Female ratio > 2 windows \n' ...
                '## Contig-name \t Contig-length \t male-female mean ratio \n'];
    fileID = fopen(OutputListFile,'w');
    fprintf(fileID,file_header);
    for i=1:numel(sorted_data)
        idx = sorted_idx(i);
        contig_name = char(cells_filtered_on_win_ratio_threshold{idx,1});
        mean_ratio = mean(male_cells_filtered_on_win_ratio_threshold{idx,3});
        fprintf(fileID,'%s \t %d \t %f \n',contig_name,sorted_data(i),mean_ratio);
    end
    fclose(fileID);
    
    % Write to IGV
    contig_list = {numel(sorted_data)};
    for i=1:numel(sorted_data)
        idx = sorted_idx(i);
        contig_list{i} = char(cells_filtered_on_win_ratio_threshold{idx,1});
    end
    batchFile = 'igv_TEST_Plut_windowBasedResults_thought_3_option_4.txt';
    Plut_createIGVBatchFile(batchFile,contig_list)
    
end
%% Create batch file for IGV
function create_igv_file(data_to_sort_names,data_to_sort)
    [sorted_data, sorted_idx] = sort(data_to_sort,'descend');
    contig_list = {numel(sorted_data)};
    for i=1:numel(sorted_data)
        contig_list{i} = data_to_sort_names{sorted_idx(i)};
    end
    batchFile = 'igv_Plut_1kbwindow_proportion_27112014.txt';
    Plut_createIGVBatchFile(batchFile,contig_list)
end

