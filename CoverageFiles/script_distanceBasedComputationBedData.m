%script_distanceBasedComputationBedData
num_features_to_select = 1000 

%% For CGR data -
male_datafile = 'DataFiles/contigs/CGR_Data/ALL_MALE.bed'
M = load_data_from_BED(male_datafile)    
female_datafile = 'DataFiles/contigs/CGR_Data/ALL_FEMALE.bed'
F = load_data_from_BED(female_datafile)

feature_list = M.features; % Feature list is in the same order for both MALE and FEMALE

MALE_DATA = [M.depth M.bases_covered_in_feature M.contig_length M.breadth];
FEMALE_DATA = [F.depth F.bases_covered_in_feature F.contig_length F.breadth];    
MALE_DATA_SCALED = [scale_between_0_1(M.depth) scale_between_0_1(M.bases_covered_in_feature) scale_between_0_1(M.contig_length) M.breadth];
FEMALE_DATA_SCALED = [scale_between_0_1(F.depth) scale_between_0_1(F.bases_covered_in_feature) scale_between_0_1(F.contig_length) F.breadth];    

D = sqrt((MALE_DATA_SCALED(:,1) - FEMALE_DATA_SCALED(:,1)).^2 + (MALE_DATA_SCALED(:,4) - FEMALE_DATA_SCALED(:,4)).^2);

[sort_D,sort_idx] = sort(D);
% Select the num_features_to_select farthest ones
CGR_contigs = feature_list(sort_idx(end-num_features_to_select:end))

%% For VT data
male_datafile = 'DataFiles/contigs/VT_Data/MALE_SRA.sorted.bed'
M = load_data_from_BED(male_datafile)    
female_datafile = 'DataFiles/contigs/VT_Data/FEMALE_SRA.sorted.bed'
F = load_data_from_BED(female_datafile)

feature_list = M.features; % Feature list is in the same order for both MALE and FEMALE

MALE_DATA = [M.depth M.bases_covered_in_feature M.contig_length M.breadth];
FEMALE_DATA = [F.depth F.bases_covered_in_feature F.contig_length F.breadth];    
MALE_DATA_SCALED = [scale_between_0_1(M.depth) scale_between_0_1(M.bases_covered_in_feature) scale_between_0_1(M.contig_length) M.breadth];
FEMALE_DATA_SCALED = [scale_between_0_1(F.depth) scale_between_0_1(F.bases_covered_in_feature) scale_between_0_1(F.contig_length) F.breadth];    

D = sqrt((MALE_DATA_SCALED(:,1) - FEMALE_DATA_SCALED(:,1)).^2 + (MALE_DATA_SCALED(:,4) - FEMALE_DATA_SCALED(:,4)).^2);

[sort_D,sort_idx] = sort(D);
% Select the num_features_to_select farthest ones
VT_contigs = feature_list(sort_idx(end-num_features_to_select:end))

%% Compare contigs
CGR_interesct_VT = intersect(CGR_contigs,VT_contigs)

%% Write results
OutputListFile = strcat('CGR_VT_OutOfTop_',num2str(num_features_to_select));
OutputListFile = strcat(OutputListFile,'.txt')

fileID = fopen(OutputListFile,'w');
fprintf(fileID,'## CGR-VT intersection - Out of %d\n',num_features_to_select)
for i=1:size(CGR_interesct_VT)
    fprintf(fileID,'%s\n',CGR_interesct_VT{i})
end

fclose(fileID);


