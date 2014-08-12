function write_to_file(out_file,dataToWrite)

    fileID = fopen(out_file,'w');
        for i=1:size(dataToWrite,1)
            fprintf(fileID,'%s\n',dataToWrite{i});
        end
    fclose(fileID);
end
