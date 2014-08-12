%filename = 'DataFiles/scaffolds/merged_bed_lane1.txt'
filename = 'DataFiles/contigs/merged_bed_lane1_contigs.txt'

delimiterIn = ',';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

% Apriori information based on the sample and filename relationship. The
% first column in the imported file indicates the samples. Samples 1 to 12
% belong to Female and samples 13-24 are for males. The following
% information indicates that which "line" in the file belong to male/female
% samples

Fem_data_rows = sort([11 17:1:24 1 2 3])
Male_data_rows = sort([4:1:10 12:1:16])

Lane_1_Fem = A.data(Fem_data_rows,:)
Lane_1_Male = A.data(Male_data_rows,:)

% Import Lane 2
%filename = 'DataFiles/merged_bed_lane2.txt'
filename = 'DataFiles/contigs/merged_bed_lane2_contigs.txt'
delimiterIn = ',';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

% Apriori information based on the sample and filename relationship. The
% first column in the imported file indicates the samples. Samples 1 to 12
% belong to Female and samples 13-24 are for males. The following
% information indicates that which "line" in the file belong to male/female
% samples

Fem_data_rows = sort([11 17:1:24 1 2 3])
Male_data_rows = sort([4:1:10 12:1:16])

Lane_2_Fem = A.data(Fem_data_rows,:)
Lane_2_Male = A.data(Male_data_rows,:)

feature_names = {A.textdata{1,2:end}}

%%
MU_lane_1_male = mean(Lane_1_Male);
SIG_lane_1_male = std(Lane_1_Male);
MU_lane_1_female = mean(Lane_1_Fem);
SIG_lane_1_female = std(Lane_1_Fem);

MU_lane_2_male = mean(Lane_2_Male);
SIG_lane_2_male = std(Lane_2_Male);
MU_lane_2_female = mean(Lane_2_Fem);
SIG_lane_2_female = std(Lane_2_Fem);

figure;
hold on
scatter(MU_lane_1_female,MU_lane_1_male)
scatter(MU_lane_2_female,MU_lane_2_male,'r')
hold off

figure;
hold on
plot(MU_lane_1_female,'*-')
plot(MU_lane_1_male,'--rs')
scrollplot

figure;
diff_mu = MU_lane_1_female - MU_lane_1_male;
plot(diff_mu)
scrollplot
figure;
plot(sort(abs(diff_mu)),'*-')

%%
% zero_val_idx_male_lane1 = find(MU_lane_1_male == 0)
% zero_val_idx_male_lane2 = find(MU_lane_2_male == 0)
% zero_val_idx_female_lane1 = find(MU_lane_1_female == 0)
% zero_val_idx_female_lane2 = find(MU_lane_2_female == 0)
% 
% zero_female_idx = union(zero_val_idx_female_lane1,zero_val_idx_female_lane2)
% zero_male_idx = union(zero_val_idx_male_lane1,zero_val_idx_male_lane2)
% zero_malefemale_intersect_idx = intersect(zero_female_idx,zero_male_idx)
% zero_F_minus_M_idx = setdiff(zero_female_idx,zero_male_idx)
% zero_M_minus_F_idx = setdiff(zero_male_idx,zero_female_idx)
% 
% zero_F_minus_M_names=feature_names(zero_F_minus_M_idx)
% zero_M_minus_F_names=feature_names(zero_M_minus_F_idx)

%%
MU_allLane_male = mean([MU_lane_1_male;MU_lane_2_male])
MU_allLane_female = mean([MU_lane_1_female;MU_lane_2_female])

zero_val_idx_male_allLane = find(MU_allLane_male == 0)
zero_val_idx_female_allLane = find(MU_allLane_female == 0)

allLane_zero_malefemale_intersect_idx = intersect(zero_val_idx_female_allLane,zero_val_idx_male_allLane)
allLane_zero_F_minus_M_idx = setdiff(zero_val_idx_female_allLane,zero_val_idx_male_allLane)
allLane_zero_M_minus_F_idx = setdiff(zero_val_idx_male_allLane,zero_val_idx_female_allLane)

allLane_zero_F_minus_M_names=feature_names(allLane_zero_F_minus_M_idx)
allLane_zero_M_minus_F_names=feature_names(allLane_zero_M_minus_F_idx)



