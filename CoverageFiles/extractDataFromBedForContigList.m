%% Given a contig list, extract the corresponding rows from the bed file for the specified column
% bedFile = 'DataFiles/contigs/VT_Data/MALE_SRA.sorted.bed'
% contigList = A cell array
% bedColumn = 6/3 etc (for breadth/depth resp.)
%
% Returns both scaled and non-scaled data
% [bedColumn_data_scaled,,contig_idx,bedColumn_data] = extractDataFromBedForContigList(bedFile,contigList,bedColumn)
function [bedColumn_data_scaled,contig_idx,bedColumn_data] = extractDataFromBedForContigList(bedFile,contigList,bedColumn)
    
    delimiterIn = '\t';
    Num_data = importdata(bedFile,delimiterIn);

    contig_idx = zeros(size(contigList,1),1);
    for i=1:size(contigList)
        contig_idx(i) = find(ismember(Num_data.rowheaders,contigList{i})==1);
    end
    
    %contig_idx = sort(contig_idx)
    
    bedColumn_data = Num_data.data(contig_idx,bedColumn);
    bedColumn_data_scaled = scale_between_0_1(bedColumn_data);
end