%script_distanceBasedComputationBedData - for Plutella
num_features_to_select = 50 

%% For CGR data -
male_datafile = 'DataFiles/DistanceBased/PLUTELLA/PLUTMALE.sorted.bed'
M = load_data_from_BED(male_datafile)    
female_datafile = 'DataFiles/DistanceBased/PLUTELLA/PLUT_FEMALE_PX9_2.sorted.bed'
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

%% Write results
OutputListFile = strcat('Plut_OutOfTop_',num2str(num_features_to_select));
OutputListFile = strcat(OutputListFile,'.txt')

fileID = fopen(OutputListFile,'w');
fprintf(fileID,'## Plutella data - Out of %d\n',num_features_to_select)
for i=1:size(CGR_contigs)
    fprintf(fileID,'%s\n',CGR_contigs{i})
end

fclose(fileID);

% For Plut_createIGVBatchFile.m
contig_list = CGR_contigs;
