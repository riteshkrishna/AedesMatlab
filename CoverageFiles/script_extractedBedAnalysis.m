contig_file1 = 'DataFiles/contigs/CGR_data/Outlier_contigs_coverage__CGR.txt'
contig_file2 = 'DataFiles/contigs/VT_data/Outlier_contigs_coverage_VT_SRA.txt'
set_operation = 'intersect'
output_dir = '/Users/ritesh/Dropbox/RAGtag_oxitec/Writeups/Analysis/VT_CGR';
outFile = strcat(output_dir,'/set_output.txt');

[contig_list] = compare_two_contig_lists(contig_file1,contig_file2,set_operation, outFile)

bedColumn_breadth = 6
bedColumn_depth = 3
bedColumn_length = 5

male_bedFile = 'DataFiles/contigs/VT_Data/MALE_SRA.sorted.bed'
female_bedFile = 'DataFiles/contigs/VT_Data/FEMALE_SRA.sorted.bed'

 
[male_breadth,male_contig_idx_breadth] = extractDataFromBedForContigList(male_bedFile,contig_list,bedColumn_breadth);
[male_depth,male_contig_idx_depth] = extractDataFromBedForContigList(male_bedFile,contig_list,bedColumn_depth);
[male_length,male_contig_idx_depth] = extractDataFromBedForContigList(male_bedFile,contig_list,bedColumn_length);
male_coverage = male_depth .*  male_breadth;
MALE_DATA = [male_breadth male_depth male_coverage male_length]

[female_breadth,female_contig_idx_breadth] = extractDataFromBedForContigList(female_bedFile,contig_list,bedColumn_breadth);
[female_depth,female_contig_idx_depth] = extractDataFromBedForContigList(female_bedFile,contig_list,bedColumn_depth);
[female_length,female_contig_idx_depth] = extractDataFromBedForContigList(male_bedFile,contig_list,bedColumn_length);
female_coverage = female_depth .* female_breadth;
FEMALE_DATA = [female_breadth female_depth female_coverage female_length]

figure;
plot(male_depth,'*-')
hold on
plot(female_depth,'ro-')
title('Depth')
legend('Male','Female')
figure_file = strcat(output_dir,'/depth.png');
print('-dtiff','-r300',figure_file)
        
figure;
plot(male_breadth,'*-')
hold on
plot(female_breadth,'ro-')
title('Breadth')
legend('Male','Female')
figure_file = strcat(output_dir,'/breadth.png');
print('-dtiff','-r300',figure_file)


figure;
plot(male_coverage,'*-')
hold on
plot(female_coverage,'ro-')
title('Coverage')
legend('Male','Female')
figure_file = strcat(output_dir,'/coverage.png');
print('-dtiff','-r300',figure_file)

figure;
plot(sort(male_depth),'*-')
hold on
plot(sort(female_depth),'ro-')
title('Coverage')
legend('Male','Female')



%%
male_bedFile = 'DataFiles/contigs/CGR_Data/ALL_MALE.bed'
female_bedFile = 'DataFiles/contigs/CGR_Data/ALL_FEMALE.bed'

[male_breadth,male_contig_idx_breadth_2] = extractDataFromBedForContigList(male_bedFile,contig_list,bedColumn_breadth);
[male_depth,male_contig_idx_depth_2] = extractDataFromBedForContigList(male_bedFile,contig_list,bedColumn_depth);
[female_breadth,female_contig_idx_breadth_2] = extractDataFromBedForContigList(female_bedFile,contig_list,bedColumn_breadth);
[female_depth,female_contig_idx_depth_2] = extractDataFromBedForContigList(female_bedFile,contig_list,bedColumn_depth);

figure;
plot(male_depth,'*-')
hold on
plot(female_depth,'ro-')
title('Depth')
legend('Male','Female')
figure_file = strcat(output_dir,'/depth.png');
print('-dtiff','-r300',figure_file)
        
figure;
plot(male_breadth,'*-')
hold on
plot(female_breadth,'ro-')
title('Breadth')
legend('Male','Female')
figure_file = strcat(output_dir,'/breadth.png');
print('-dtiff','-r300',figure_file)
male_coverage = male_depth .*  male_breadth;
female_coverage = female_depth .* female_breadth;
figure;
plot(male_coverage,'*-')
hold on
plot(female_coverage,'ro-')
title('Coverage')
legend('Male','Female')
figure_file = strcat(output_dir,'/coverage.png');
print('-dtiff','-r300',figure_file)
