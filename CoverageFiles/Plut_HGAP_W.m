%% New W data available and mapped against HGAP_NOV_2014 assembly.
function Plut_HGAP_W()
    
    prepare_all_datasets_and_result_files()
    
end

function analyse_results()
    
%     load('DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_result_pxPP6.mat')
%     for i=1:numel(left_right_ratio_pass_candidates_with_enoughWindows_enoughLongC)
%         plot(M.depth(feature_depth{left_right_ratio_pass_candidates_with_enoughWindows_enoughLongC(end-i),2}))
%         pause
%     end

%     % Each one contains - contig_name contig_length no_of_windows_occupied
%     contigs_GP9 = importdata('DataFiles/DistanceBased/PLUTELLA/HGAP/output_pxGP9.txt');
%     contigs_PP6 = importdata('DataFiles/DistanceBased/PLUTELLA/HGAP/output_pxPP6.txt');
%     contigs_PP4 = importdata('DataFiles/DistanceBased/PLUTELLA/HGAP/output_pxPP4.txt');
%     contigs_PP7 = importdata('DataFiles/DistanceBased/PLUTELLA/HGAP/output_pxPP7.txt');
%     % Find all common contig names
%     common_contigs = intersect(intersect(intersect(contigs_GP9.rowheaders, contigs_PP6.rowheaders), contigs_PP4.rowheaders), contigs_PP7.rowheaders);
%     
    % Alternatively, load the result mat files and plot from there..
    pp6_results = load('DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_result_pxPP6.mat');
    pp4_results = load('DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_result_pxPP4.mat');
    gp9_results = load('DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_result_pxGP9.mat');
    pp7_results = load('DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_result_pxPP7.mat');
    
    pp6_contig_names = pp6_results.feature_depth(pp6_results.left_right_ratio_pass_candidates_with_enoughWindows_enoughLongC,1);
    pp4_contig_names = pp4_results.feature_depth(pp4_results.left_right_ratio_pass_candidates_with_enoughWindows_enoughLongC,1);
    gp9_contig_names = gp9_results.feature_depth(gp9_results.left_right_ratio_pass_candidates_with_enoughWindows_enoughLongC,1);
    pp7_contig_names = pp7_results.feature_depth(pp7_results.left_right_ratio_pass_candidates_with_enoughWindows_enoughLongC,1);

    [common_contigs_all_4dsets] = intersect(intersect(intersect(gp9_contig_names, pp6_contig_names), pp4_contig_names),pp7_contig_names);
    
    for i=1:numel(common_contigs_all_4dsets)
        contig_name = common_contigs_all_4dsets{i};
        
        pp4_idx = find(strcmp(contig_name,pp4_results.feature_depth(:,1)))
        pp4_win_idx =  pp4_results.feature_depth(pp4_idx,2);
        pp4_data = pp4_results.M.depth(pp4_win_idx{:});
        
        pp6_idx = find(strcmp(contig_name,pp6_results.feature_depth(:,1)))
        pp6_win_idx =  pp6_results.feature_depth(pp6_idx,2);
        pp6_data = pp6_results.M.depth(pp6_win_idx{:});
        
        pp7_idx = find(strcmp(contig_name,pp7_results.feature_depth(:,1)))
        pp7_win_idx =  pp7_results.feature_depth(pp7_idx,2);
        pp7_data = pp7_results.M.depth(pp7_win_idx{:});
        
        gp9_idx = find(strcmp(contig_name,gp9_results.feature_depth(:,1)))
        gp9_win_idx =  gp9_results.feature_depth(gp9_idx,2);
        gp9_data = gp9_results.M.depth(gp9_win_idx{:});
        
        subplot(2,2,1);
        stem(pp4_data);
        title('PP4')
        subplot(2,2,2);
        stem(pp6_data);
        title('PP6')
        subplot(2,2,3);
        stem(pp7_data);
        title('PP7')
        subplot(2,2,4);
        stem(gp9_data);
        title('GP9')
        pause
    end
    
end
function prepare_all_datasets_and_result_files()

    [matfile_name,result_matfile_name] = prepare_pxPP6();
    analyze_data(matfile_name, 'DataFiles/DistanceBased/PLUTELLA/HGAP/output_verbose_pxPP6.txt',...
                               'DataFiles/DistanceBased/PLUTELLA/HGAP/output_pxPP6.txt',result_matfile_name);
    
    [matfile_name,result_matfile_name] = prepare_pxPP4();
    analyze_data(matfile_name, 'DataFiles/DistanceBased/PLUTELLA/HGAP/output_verbose_pxPP4.txt',...
                                'DataFiles/DistanceBased/PLUTELLA/HGAP/output_pxPP4.txt',result_matfile_name);
   
    [matfile_name,result_matfile_name] = prepare_pxGP9()
    analyze_data(matfile_name, 'DataFiles/DistanceBased/PLUTELLA/HGAP/output_verbose_pxGP9.txt',...
                                'DataFiles/DistanceBased/PLUTELLA/HGAP/output_pxGP9.txt',result_matfile_name);
    
    [matfile_name,result_matfile_name] = prepare_pxPP7()
    analyze_data(matfile_name, 'DataFiles/DistanceBased/PLUTELLA/HGAP/output_verbose_pxPP7.txt',...
                                'DataFiles/DistanceBased/PLUTELLA/HGAP/output_pxPP7.txt',result_matfile_name);
   
end

function [matfile_name,result_matfile_name] = prepare_pxPP6()
    input_bed = 'DataFiles/DistanceBased/PLUTELLA/HGAP/PLUTHGAP.sorted.window_1kb.bed';
    matfile_name = 'DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_data_pxPP6.mat';
    result_matfile_name = 'DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_result_pxPP6.mat';
    prepare_pxdata(input_bed, matfile_name);
end

function [matfile_name,result_matfile_name] = prepare_pxPP4()
    input_bed = 'DataFiles/DistanceBased/PLUTELLA/HGAP/PLUTHGAP.pxPP4.sorted.window_1kb.bed';
    matfile_name = 'DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_data_pxPP4.mat';
    result_matfile_name = 'DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_result_pxPP4.mat';
    prepare_pxdata(input_bed, matfile_name);
end

function [matfile_name,result_matfile_name] = prepare_pxGP9()
    input_bed = 'DataFiles/DistanceBased/PLUTELLA/HGAP/PLUTHGAP.pxGP9.sorted.window_1kb.bed';
    matfile_name = 'DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_data_pxGP9.mat';
    result_matfile_name = 'DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_result_pxGP9.mat';
    prepare_pxdata(input_bed, matfile_name);
end

function [matfile_name,result_matfile_name] = prepare_pxPP7()
    input_bed = 'DataFiles/DistanceBased/PLUTELLA/HGAP/PLUTHGAP.pxPP7.sorted.window_1kb.bed';
    matfile_name = 'DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_data_pxPP7.mat';
    result_matfile_name = 'DataFiles/DistanceBased/PLUTELLA/HGAP/workspace_result_pxPP7.mat';
    prepare_pxdata(input_bed, matfile_name);
end


function prepare_pxdata(input_bed, matfile_name)

    M = struct
    delimiterIn = '\t';
    A = importdata(input_bed,delimiterIn);
    M.features = A.textdata;
    M.depth = A.data(:,3);
    M.bases_covered_in_feature = A.data(:,4);
    M.contig_length = A.data(:,5);
    M.breadth = A.data(:,6);

    % Get names of features that are not zero-covered
    non_zerodepth_idx = find(M.depth ~= 0);
    non_zerodepth_feature = unique(M.features(non_zerodepth_idx));

    feature_depth = cell(numel(non_zerodepth_feature),2);

    for i = 1:numel(non_zerodepth_feature)
        %idx_for_this_feature = find(~cellfun(@isempty,strfind(M.features, non_zerodepth_feature{i})));
        idx_for_this_feature = find(strcmp(non_zerodepth_feature{i},M.features));
        feature_depth{i,1} = non_zerodepth_feature{i};
        feature_depth{i,2} = idx_for_this_feature;
    end
    
    save(matfile_name)
end


function analyze_data(mat_workspace, output_file_verbose, output_file_result,result_matfile_name)
    
    % We expect to see at least as many as these windows with non-zero
    % values
    MIN_FILLED_WINDOW_THRESHOLD = 20;
    MIN_CONTIGLENGTH_THRESHOLD = 20000;
    
    
    load(mat_workspace);
    
    fout_verbose = fopen(output_file_verbose,'w')
    fout_results = fopen(output_file_result,'w')
    
    left_right_ratio_pass_candidates = [];
    left_right_ratio_pass_candidates_with_enoughWindows = [];
    left_right_ratio_pass_candidates_with_enoughWindows_enoughLongContig = [];
    
    for i=1:size(feature_depth,1)
        feature_name = feature_depth{i,1};
        window_list = M.depth(feature_depth{i,2});
        contig_length = sum(M.contig_length(feature_depth{i,2}));
        [pass_status_occ,out_st1] = check_occurence_pattern(i, feature_name, window_list);
        [pass_status_spread, out_st2, flags] = check_spread_pattern(i,feature_name, window_list);
        
        if (pass_status_occ == 1) & (pass_status_spread == 1)
            left_right_ratio_pass_candidates(end+1,1) = i;
            
            total_occupied_windows = sum(sum(flags));
            if total_occupied_windows >= MIN_FILLED_WINDOW_THRESHOLD
                left_right_ratio_pass_candidates_with_enoughWindows(end+1,1) = i;
                
                if contig_length >= MIN_CONTIGLENGTH_THRESHOLD
                    left_right_ratio_pass_candidates_with_enoughWindows_enoughLongContig(end+1,1) = i;
                    
                    out_str = sprintf('%s \t %d \t %d \n', feature_name, contig_length, total_occupied_windows);
                    fprintf(fout_results,'%s', out_str);
                end
                
            end
            
        end
        
        fprintf(fout_verbose,'%s', out_st1);
        fprintf(fout_verbose,'%s', out_st2);
    end
    
    fclose(fout_verbose);
    fclose(fout_results);
    
    save(result_matfile_name);
    
%     for i=1:441
%         plot(M.depth(feature_depth{left_right_ratio_pass_candidates(end-i),2}))
%         pause
%     end

    
end

function [pass_status,out_string] = check_occurence_pattern(counter, feature_name, window_list)
    mid_point = round(numel(window_list)/2);
    win_left = numel(find(window_list(1:mid_point) ~= 0));
    win_right = numel(find(window_list(mid_point+1:end) ~= 0));
    
    ratio_LeftRight = win_left / win_right;
    
    threshold = 0.25; % less than 25% difference means skewd on left or right
    if (ratio_LeftRight < threshold) | (ratio_LeftRight > 1/threshold)
        pass_status = 0;
    else
        pass_status = 1;
    end
    
    out_string = sprintf('%d :occ %d => \t %s \t %f \n',counter, pass_status,feature_name, ratio_LeftRight);
    %fprintf(1,'%d :occ %d => \t %s \t %f \n',counter, pass_status,feature_name, ratio_LeftRight)
    fprintf(1,'%s',out_string);
end

function  [pass_status,out_string, flags] = check_spread_pattern(counter,feature_name, window_list)
    
    mid_point = round(numel(window_list) * 0.5);
    left_part = window_list(1:mid_point);
    right_part = window_list(mid_point:end);
    
    win = left_part;
    quarter_point = round(numel(win) * 0.25);
    mid_point = round(numel(win) * 0.5);
    threequarter_point = round(numel(win) * 0.75);
    win_quarter = numel(find(win(1:quarter_point) ~= 0));
    win_quarter_mid = numel(find(win(quarter_point+1:mid_point) ~= 0));
    win_mid_threequarter = numel(find(win(mid_point+1:threequarter_point) ~= 0));
    win_threequarter_end = numel(find(win(threequarter_point+1:end) ~= 0));
    left_flag = [(win_quarter > 0) (win_quarter_mid > 0) (win_mid_threequarter > 0) (win_threequarter_end > 0) ];
    flags =  [win_quarter win_quarter_mid win_mid_threequarter win_threequarter_end];
    
    win = right_part;
    quarter_point = round(numel(win) * 0.25);
    mid_point = round(numel(win) * 0.5);
    threequarter_point = round(numel(win) * 0.75);
    win_quarter = numel(find(win(1:quarter_point) ~= 0));
    win_quarter_mid = numel(find(win(quarter_point+1:mid_point) ~= 0));
    win_mid_threequarter = numel(find(win(mid_point+1:threequarter_point) ~= 0));
    win_threequarter_end = numel(find(win(threequarter_point+1:end) ~= 0));
    right_flag = [(win_quarter > 0) (win_quarter_mid > 0) (win_mid_threequarter > 0) (win_threequarter_end > 0)];
    flags = [flags; win_quarter win_quarter_mid win_mid_threequarter win_threequarter_end];
    
    
    flag_and = left_flag & right_flag;
    
    if sum(flag_and) >= 2
        pass_status = 1;
    else
        pass_status = 0;
    end
    %fprintf(1,'%d :spr %d => \t %s \t %s\n',counter, pass_status,feature_name, mat2str(flag_and))
    out_string = sprintf('%d :spr %d => \t %s \t %s\n',counter, pass_status,feature_name, mat2str(flag_and));
    fprintf(1,'%s',out_string);
end

