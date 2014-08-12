%% Coverage
% contig_file1 = 'DataFiles/contigs/CGR_data/Outlier_contigs_coverage__CGR.txt'
% contig_file2 = 'DataFiles/contigs/VT_data/Outlier_contigs_coverage_VT_SRA.txt'
% set_operation = 'intersect'/'A_minus_B'/'B_minus_A'/'union'
% outFile = 'set_output.txt'
%% Breadth
% contig_file1 = 'DataFiles/contigs/CGR_data/Outlier_contigs_breadth_CGR.txt'
% contig_file2 = 'DataFiles/contigs/VT_data/Outlier_contigs_breadth_VT_SRA.txt'
% set_operation = 'intersect'/'A_minus_B'/'B_minus_A'/'union'
% outFile = 'set_output.txt'
%%

function [contig_list] = compare_two_contig_lists(contig_file1,contig_file2,set_operation, outFile)
    
    [contig_list] = compare_contig_files(contig_file1, contig_file2,set_operation);
    
    sprintf('\n Count : %s \t %d',set_operation, size(contig_list,1))
    
    write_to_file(outFile, contig_list);
end

function [res] = compare_contig_files(contig_file1, contig_file2,operation)
    
    A = importdata(contig_file1);
    B = importdata(contig_file2);
    
    switch(operation)
        case 'union'
            res = union(A,B)
        case 'intersect'
            res = intersect(A,B)
        case 'A_minus_B'
            res = setdiff(A,B)
        case 'B_minus_A'
            res = setdiff(B,A)
        otherwise
                sprintf('ERROR')
    end
    
end
