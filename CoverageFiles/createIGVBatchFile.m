%% Create an IGV batch file for the list of contigs we are interested in looking at

genome = 'AedesAegypti'
bam_file = '/pub9/ritesh/Aedes_Aegypti/BowtieMappings_single/SRA_data/FEMALE_SRA.sorted.bam,/pub9/ritesh/Aedes_Aegypti/BowtieMappings_single/SRA_data/MALE_SRA.sorted.bam '
output_dir = '/pub9/ritesh/Aedes_Aegypti/BowtieMappings_single/SRA_data/igv_snapshots'

batchFile = 'igv_batch_35_CGR_VT.txt'

fileID = fopen(batchFile,'w');
fprintf(fileID,'new \n');
fprintf(fileID,'load %s \n',bam_file);
fprintf(fileID,'snapshotDirectory %s \n',output_dir);
fprintf(fileID,'genome %s \n',genome);

for i=1:size(contig_list,1)
     fprintf(fileID,'goto %s\n',contig_list{i});
     fprintf(fileID,'sort base \n');
     fprintf(fileID,'collapse \n');
     fprintf(fileID,'snapshot \n');
end
fclose(fileID);
    
