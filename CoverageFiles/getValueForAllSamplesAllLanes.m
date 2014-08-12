% For a given index, check what's the coverage like for all samples across
% all lanes
function [filtered_idx] = getValueForAllSamplesAllLanes(allLane_zero_idx,Lane_1_Male,...,
                            Lane_2_Male,Lane_1_Fem,Lane_2_Fem,gender)

    t = size(allLane_zero_idx,2)
    filtered_idx = [];
    
    SAMPLE_THRESHOLD = 3
    
    for n=1:t
        idx = allLane_zero_idx(n);
        M_samples = [Lane_1_Male(:,idx) Lane_2_Male(:,idx)];
        F_samples = [Lane_1_Fem(:,idx) Lane_2_Fem(:,idx)];
        
        non_zero_M_samples = sum(sum(M_samples > 0));
        non_zero_F_samples = sum(sum(F_samples > 0));
        
        if gender == 'M'
            if non_zero_F_samples >= SAMPLE_THRESHOLD
                filtered_idx = [filtered_idx; idx];
            end
        elseif gender == 'F'
            if non_zero_M_samples >= SAMPLE_THRESHOLD
                filtered_idx = [filtered_idx; idx];
            end
        end
        
    end
    
end