function Plut_file_write(OutputListFile,data_to_sort_names,data_to_sort, file_header)

    [sorted_data, sorted_idx] = sort(data_to_sort,'descend');
    %OutputListFile = 'Plut_windowBasedResults_proportion.txt'
    fileID = fopen(OutputListFile,'w');
    fprintf(fileID,file_header);
    %fprintf(fileID,'## Plutella window based identifications - 1KB window size for analysis \n');
    %fprintf(fileID,'## Contig names \t proportion of contig in the region of interest \n');
    for i=1:numel(sorted_data)
        fprintf(fileID,'%s \t %d \n',data_to_sort_names{sorted_idx(i)},sorted_data(i));
    end
    fclose(fileID);
end
