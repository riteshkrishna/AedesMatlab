%% Create an IGV batch file for the list of contigs we are interested in looking at

genome = 'Plutella_polishedAssembly'
bam_file = '/pub9/ritesh/Plutella/BowtieMappings/polished_assembly/PLUT_FEMALE_PX9_2.sorted.bam,/pub9/ritesh/Plutella/BowtieMappings/polished_assembly/PLUTMALE.sorted.bam '
output_dir = '/pub9/ritesh/Plutella/Analysis/polished_assembly'

batchFile = 'igv_plut_w_candidates_03112014.txt'

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
    
