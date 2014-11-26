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

%% Z and W analysis
% W candidates = male features with zero depth and non-zero depth in female

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

%% Find male-W now - Simpler way
W_male_idx = find(Filtered_MALE_DATA(:,col_depth) == 0); 

h = figure;
plot(Filtered_MALE_DATA(W_male_idx,col_depth),'*-')
hold on
plot(Filtered_FEMALE_DATA(W_male_idx,col_depth),'ro-')
title('W candidates from zero males and non-zero females')
legend('male coverage','female coverage')
xlabel('contigs');
ylabel('Depth');
%saveas(h,'w_zeromale_nonzerofemale','png');
h = figure;
plot(Filtered_MALE_DATA(W_male_idx,col_featureLength),'bs-'); % W candidate contig_lengths
title('W candidates from zero males and non-zero females')
xlabel('contigs');
ylabel('size');
%saveas(h,'w_zeromale_contig_lengths','png');

W_contig_names = filtered_feature_list(W_male_idx);
W_contig_length = Filtered_MALE_DATA(W_male_idx,col_featureLength);
W_corresponding_Female_depth = Filtered_FEMALE_DATA(W_male_idx,col_depth);

% out_file = 'Plut_W_candidates.txt'
% fileID = fopen(out_file,'w');
% fprintf(fileID,'##Selection of W candidates that have zero coverage in Male samples. \n## Names \t lengths \t coverage in female samples \n');
% for i=1:size(W_contig_names,1)
%     fprintf(fileID,'%s \t %d \t %d \n',W_contig_names{i}, W_contig_length(i),round(W_corresponding_Female_depth(i)));
% end
% fclose(fileID);

%% Find opposite of above - Present in male, absent in female
opp_W_male_idx = find(Filtered_FEMALE_DATA(:,col_depth) == 0); 

h = figure;
plot(Filtered_FEMALE_DATA(opp_W_male_idx,col_depth),'*-')
hold on
plot(Filtered_MALE_DATA(opp_W_male_idx,col_depth),'ro-')
title('Candidates from zero females and non-zero males')
legend('female coverage','male coverage')
xlabel('contigs');
ylabel('Depth');
%saveas(h,'zerofemale_nonzeromale','png');
h = figure;
plot(Filtered_FEMALE_DATA(opp_W_male_idx,col_featureLength),'bs-'); % W candidate contig_lengths
title('Candidates from zero females and non-zero males')
xlabel('contigs');
ylabel('size');
%saveas(h,'zerofemale_contig_lengths','png');

opp_W_contig_names = filtered_feature_list(opp_W_male_idx);
opp_W_contig_length = Filtered_FEMALE_DATA(opp_W_male_idx,col_featureLength);
opp_W_corresponding_male_depth = Filtered_MALE_DATA(opp_W_male_idx,col_depth);

% out_file = 'Plut_oppositeOf_W_candidates.txt'
% fileID = fopen(out_file,'w');
% fprintf(fileID,'##Selection of candidates that have zero coverage in female samples. \n## Names \t lengths \t coverage in male samples \n');
% for i=1:size(opp_W_contig_names,1)
%     fprintf(fileID,'%s \t %d \t %d \n',opp_W_contig_names{i}, opp_W_contig_length(i),round(opp_W_corresponding_male_depth(i)));
% end
% fclose(fileID);

%% Look at the distributions
% Remove zero male and zero females as we know them from above W, opp-W
% analyses
zero_male_female_indices = [W_male_idx; opp_W_male_idx];
Filtered_MALE_DATA(zero_male_female_indices,:) = [];
Filtered_FEMALE_DATA(zero_male_female_indices,:) = [];
filtered_feature_list(zero_male_female_indices,:) = [];

male_depth = log10(Filtered_MALE_DATA(:,col_depth));
female_depth = log10(Filtered_FEMALE_DATA(:,col_depth));

%%
figure;
h1 = histogram(male_depth,'Normalization','pdf');
hold on
y = 0:0.1:7;
mu = mean(male_depth);
sigma = std(male_depth);
f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(y,f,'LineWidth',1.5)

figure;
h2 = histogram(female_depth,'Normalization','pdf');
hold on
y = 0:0.1:7;
mu = mean(female_depth);
sigma = std(female_depth);
f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(y,f,'LineWidth',1.5)

figure;
h1 = histogram(male_depth,'Normalization','pdf');
hold on
h2 = histogram(female_depth,'Normalization','pdf');
y = 0:0.1:7;
mu = mean(female_depth) + mean(male_depth);
mu = mu / 2;
sigma = std(female_depth) + std(male_depth);
sigma = sigma / 2;
f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(y,f,'LineWidth',1.5)
