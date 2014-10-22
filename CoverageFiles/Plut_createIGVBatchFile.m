%% Create an IGV batch file for the list of PLUTELLA contigs we are interested in looking at

genome = 'Plutella_RB_assemb'
bam_file = '/pub9/ritesh/Plutella/BowtieMappings/scaffold_contig/PLUT_FEMALE_PX9_2.sorted.bam,/pub9/ritesh/Plutella/BowtieMappings/scaffold_contig/PLUTMALE.sorted.bam '
output_dir = '/pub9/ritesh/Plutella/Analysis/scaffold_contig/igv_snapshots'

batchFile = 'igv_Plut_OutOfTop_200_19092014.txt'

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
    
