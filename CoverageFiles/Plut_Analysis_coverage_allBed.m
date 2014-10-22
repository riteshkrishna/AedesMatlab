%Analysis of coverage_allBed.txt, coverage here means "breadth only" as
%reported by coverageBED. 

%% % Settings For plut data
filename = 'DataFiles/DistanceBased/PLUTELLA/scf_contigs/coverage_allBed.txt'

% For Al's assembly
%filename = 'DataFiles/DistanceBased/PLUTELLA/PX_contig/coverage_PX_contig_allBed.txt'
%filename = 'DataFiles/DistanceBased/PLUTELLA/PX_scaffold/coverage_PX_scaffold_allBed.txt'

% For US assembly
%filename = 'DataFiles/DistanceBased/PLUTELLA/US_assembly/coverage_polished_assembly_allBed.txt'

Fem_data_rows = [1:1:4]
Male_data_rows = [5]

%% Generic code
delimiterIn = ',';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

feature_names = {A.textdata{1,2:end}}

Fem_samples = A.data(Fem_data_rows,:)
Male_samples = A.data(Male_data_rows,:)

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


