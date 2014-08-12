%% Return the content of a BED file in a struct. 
function [M] = load_data_from_BED(male_datafile)
    M = struct
    delimiterIn = '\t';
    A = importdata(male_datafile,delimiterIn);

    M.features = A.textdata;
    M.depth = A.data(:,3);
    M.bases_covered_in_feature = A.data(:,4);
    M.contig_length = A.data(:,5);
    M.breadth = A.data(:,6);
end

