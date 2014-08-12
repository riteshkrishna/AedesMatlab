%Analysis of coverage_allBed.txt, coverage here means "breadth only" as
%reported by coverageBED. 

%% % Settings For CGR data
%filename = 'DataFiles/contigs/coverage_allBed.txt'
% all_Fem_data_row = 1
% all_Male_data_row = 2
% Fem_data_rows = [3:1:14]
% Male_data_rows = [15:1:26]


%% % Settings For Virginia Tech SRA data
filename = 'DataFiles/contigs/VT_Data/merged_VT_SRA.txt'
all_Fem_data_row = 1
all_Male_data_row = 2
Fem_data_rows = [3:1:4]
Male_data_rows = [5:1:6]

%% Generic code
delimiterIn = ',';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

feature_names = {A.textdata{1,2:end}}

Fem_samples = A.data(Fem_data_rows,:)
Male_samples = A.data(Male_data_rows,:)
all_Fem = A.data(all_Fem_data_row,:)
all_Male = A.data(all_Male_data_row,:)

% Plot individual samples - Female
figure;
for i=1:size(Fem_samples,1)
    subplot(4,3,i)
    hist(Fem_samples(i,:),100);
    sample_name = A.textdata(Fem_data_rows(i) + 1)
    % Search-replace '_' with a '-' to make the figure titles looks better
    expression = '_';
    replace = '-';
    sample_name = regexprep(sample_name,expression,replace)
    title(sample_name)
end

% Plot individual samples - Male
figure;
for i=1:size(Male_samples,1)
    subplot(4,3,i)
    hist(Male_samples(i,:),100);
    sample_name = A.textdata(Male_data_rows(i) +  1)
    
    expression = '_'; 
    replace = '-';
    sample_name = regexprep(sample_name,expression,replace)
    title(sample_name)
end

% Plot global samples
figure;
hist(all_Fem,100)
h = title('Contigs: All female Sample combined together')
set(h,'FontSize',14)
h = xlabel('Fraction of breadth covered')
set(h,'FontSize',14)

figure;
hist(all_Male,100)
h = title('Contigs: All male Sample combined together')
set(h,'FontSize',14)
h = xlabel('Fraction of breadth covered')
set(h,'FontSize',14)

figure;
scatter(all_Male,all_Fem)
h = xlabel('All Male');
set(h,'FontSize',14)
h = ylabel('All Female');
set(h,'FontSize',14)
