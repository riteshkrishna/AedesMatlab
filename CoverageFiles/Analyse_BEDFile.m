%% Code to analyze BED files using linear regression
%% For VT data -
% male_datafile = 'DataFiles/contigs/VT_Data/MALE_SRA.sorted.bed'
% female_datafile = 'DataFiles/contigs/VT_Data/FEMALE_SRA.sorted.bed'
% operation = 'coverage'/'breadth'/'depth'
% out_file = 'Outlier_contigs.txt';
%% For CGR data -
% male_datafile = 'DataFiles/contigs/CGR_Data/ALL_MALE.bed'
% female_datafile = 'DataFiles/contigs/CGR_Data/ALL_FEMALE.bed'
% operation = 'coverage'/'breadth'/'depth'
% out_file = 'Outlier_contigs.txt';
%%
function [M_operation,F_operation] = Analyse_BEDFile(male_datafile,female_datafile,operation,out_file)
    
    [M,F] = load_data(male_datafile,female_datafile);
    [Data_filtered] = cleanData(M,F);
    [Data_scaled] = scaleData(Data_filtered);
    
    M_operation = 0
    F_operation = 0
    switch (operation)
        case 'coverage'
            % Defining coverage as product of depth and breadth
            M_coverage = Data_scaled.male_depth_filtered .* Data_scaled.male_breadth_filtered;
            F_coverage = Data_scaled.female_depth_filtered .* Data_scaled.female_breadth_filtered;
            M_operation = M_coverage
            F_operation = F_coverage
        case 'breadth'
            M_operation = Data_scaled.male_breadth_filtered;
            F_operation = Data_scaled.female_breadth_filtered;
        case 'depth'
            M_operation = Data_scaled.male_depth_filtered;
            F_operation = Data_scaled.female_depth_filtered;
        otherwise
            sprintf('ERROR')
    end
    
    resids = fit_data(M_operation, F_operation,operation);
    z = 3;
    outlier_idx = filter_residuals(resids,z)
    analyse_outliers(M_operation,F_operation,outlier_idx,Data_scaled,out_file)
            
    %print ' '
end

function [M,F] = load_data(male_datafile,female_datafile)
    M = struct
    delimiterIn = '\t';
    A = importdata(male_datafile,delimiterIn);

    M.male_features = A.textdata;
    M.male_depth = A.data(:,3);
    M.male_bases_covered_in_feature = A.data(:,4);
    M.male_contig_length = A.data(:,5);
    M.male_breadth = A.data(:,6);
    
    %% FEMALE data
    F = struct
    delimiterIn = '\t';
    B = importdata(female_datafile,delimiterIn);

    F.female_features = B.textdata;
    F.female_depth = B.data(:,3);
    F.female_bases_covered_in_feature = B.data(:,4);
    F.female_contig_length = B.data(:,5);
    F.female_breadth = B.data(:,6);
    
    plot_figures(M,F,1)
    
end

%%
function plot_figures(M,F,type)
   
    switch type
        case 1
            %%
            figure;
            hist(M.male_breadth,100)
            figure;
            hist(F.female_breadth,100)
        case 2
            
    end

end

%%
% Overall Coverage over a contig should be defined as the product of depth
% and breadth. The deapth shall first be normalized before a product is
% taken.
% Input from - [Data_filtered] = cleanData(M,F)
function [Data_scaled] = scaleData(Data_filtered)
    Data_scaled = Data_filtered;
    
    F_data = Data_scaled.female_depth_filtered;
    M_data = Data_scaled.male_depth_filtered;
    
    [F_scaled] = scale_between_0_1(F_data)
    [M_scaled] = scale_between_0_1(M_data)
    
    Data_scaled.female_depth_filtered = F_scaled
    Data_scaled.male_depth_filtered = M_scaled;
end

function [Data_filtered] = cleanData(M,F)
    [temp,I_female] = fun_removeOutliers(F.female_depth);
    [temp,I_male] = fun_removeOutliers(M.male_depth);

    outlier_female_idx = find(I_female == 0);
    outlier_male_idx = find(I_male == 0);
    all_outlier_idx = union(outlier_male_idx,outlier_female_idx);
    common_outlier_idx = intersect(outlier_male_idx,outlier_female_idx);
    outlier_female_specific = setdiff(outlier_female_idx,outlier_male_idx);
    outlier_male_specific = setdiff(outlier_male_idx,outlier_female_idx);
    
    plot_outlier(M,F,all_outlier_idx)
    
    % Both the profiles for both the conditions look same, so we can remove
    % them
    Data_filtered = struct
    Data_filtered.contigs_to_keep = all((I_female & I_male),2);
    idx_contigs_to_keep = find(Data_filtered.contigs_to_keep==1);
    
    Data_filtered.female_depth_filtered = F.female_depth(idx_contigs_to_keep,:);
    Data_filtered.male_depth_filtered   = M.male_depth(idx_contigs_to_keep,:);
    
    Data_filtered.female_breadth_filtered = F.female_breadth(idx_contigs_to_keep,:);
    Data_filtered.male_breadth_filtered   = M.male_breadth(idx_contigs_to_keep,:);
    
    Data_filtered.male_contigs_filtered = M.male_features(idx_contigs_to_keep,:)
    Data_filtered.female_contigs_filtered = F.female_features(idx_contigs_to_keep,:)
end

%%
% Plot male and female outliers wrt each other
function plot_outlier(M,F,all_outlier_idx)
    figure
    plot(log(F.female_depth(all_outlier_idx)),'o-')
    hold on
    plot(log(M.male_depth(all_outlier_idx)),'rs-')
    figure
    plot(log(F.female_breadth(all_outlier_idx)),'o-')
    hold on
    plot(log(M.male_breadth(all_outlier_idx)),'rs-')
end

%%
function [resid] = fit_data(M_data,F_data,title_label)
    
    [p,ErrorEst] = polyfit(M_data,F_data,1);
    M_hat = polyval(p,M_data,ErrorEst);
    resid = F_data - M_hat;
    
    plot_fit(M_data,F_data,resid,title_label)
end
%%
function plot_fit(M_data,F_data,resid,title_label)
    figure;
    plot(M_data,F_data,'.')
    %grid on
    %h = title('Contigs : Coverage as a product of deapth and breadth');
    h = title(title_label);
    set(h,'FontSize',14);
    h = xlabel('Male');
    set(h,'FontSize',14)
    h = ylabel('Female');
    set(h,'FontSize',14)
    
    figure;
    plot(M_data, resid,'.')
    h = title('Residuals analysis')
    set(h,'FontSize',14)
    h = xlabel('Male');
    set(h,'FontSize',14)
    h = ylabel('Residuals');
    set(h,'FontSize',14)   

end
%% Find the outliers among residuals, the function returns indices of
% residuals. z is 1,2 or 3 corresponding to signam, 2-sigma, 3-sigma etc
function [outlier_idx] = filter_residuals(resids,z)
    n = size(resids,1);
    mu_resid = mean(resids);
    std_resid = std(resids);
    I = abs(resids - repmat(mu_resid,n,1)) <= z .* repmat(std_resid,n,1);
    outlier_idx = find(I == 0);
end

%%
function analyse_outliers(M_coverage,F_coverage,outlier_idx,Data_scaled,out_file)
    
    figure;
    scatter(M_coverage(outlier_idx),F_coverage(outlier_idx));
    h = title('Contigs : Outliers');
    set(h,'FontSize',14)
    h = xlabel('Male');
    set(h,'FontSize',14)
    h = ylabel('Female');
    set(h,'FontSize',14)
    
    sprintf('Total outliers : %d', size(outlier_idx,1))
    
    outlier_contigs = Data_scaled.male_contigs_filtered(outlier_idx);
    
    write_to_file(out_file,outlier_contigs)
    
end

function temp()
%%
%http://www.mathworks.co.uk/help/bioinfo/ug/exploring-genome-wide-differences-in-dna-methylation-profiles.html

counts_1 = female_depth;
counts_2 = male_depth;

nbp = nbinfit(counts_1);
figure
hold on
emphist = histc(counts_1,0:100); % calculate the empirical distribution
bar(0:100,emphist./sum(emphist),'c','grouped') % plot histogram
plot(0:100,nbinpdf(0:100,nbp(1),nbp(2)),'b','linewidth',2); % plot fitted model
%axis([0 50 0 .001])
legend('Empirical Distribution','Negative Binomial Fit')
ylabel('Frequency')
xlabel('Counts')
title('Frequency of counts for 100 bp windows (HCT116-1)')

rtnbinpdf = @(x,p1,p2,t) nbinpdf(x,p1,p2) ./ nbincdf(t-1,p1,p2);
rtnbinfit = @(x,t) mle(x,'pdf',@(x,p1,p2) rtnbinpdf(x,p1,p2,t),'start',nbinfit(x),'lower',[0 0]);
nbp = [0.5 0.2];              % Known coefficients
x = nbinrnd(nbp(1),nbp(2),10000,1); % Random sample
trun = 6;                     % Set a truncation threshold

nbphat1 = nbinfit(x);         % Fit non-truncated model to all data
nbphat2 = nbinfit(x(x<trun)); % Fit non-truncated model to truncated data (wrong)
nbphat3 = rtnbinfit(x(x<trun),trun); % Fit truncated model to truncated data

figure
hold on
emphist = histc(x,0:100);     % Calculate the empirical distribution
bar(0:100,emphist./sum(emphist),'c','grouped') % plot histogram
h1 = plot(0:100,nbinpdf(0:100,nbphat1(1),nbphat1(2)),'b-o','linewidth',2);
h2 = plot(0:100,nbinpdf(0:100,nbphat2(1),nbphat2(2)),'r','linewidth',2);
h3 = plot(0:100,nbinpdf(0:100,nbphat3(1),nbphat3(2)),'g','linewidth',2);
axis([0 25 0 .2])
legend([h1 h2 h3],'Neg-binomial fitted to all data',...
                  'Neg-binomial fitted to truncated data',...
                  'Truncated neg-binomial fitted to truncated data')
ylabel('Frequency')
xlabel('Counts')

trun = 4;  % Set a truncation threshold (as in [1])
pn1 = rtnbinfit(counts_1(counts_1<trun),trun); % Fit to HCT116-1 counts
pn2 = rtnbinfit(counts_2(counts_2<trun),trun); % Fit to HCT116-2 counts

pval1 = 1 - nbincdf(counts_1,pn1(1),pn1(2));
pval2 = 1 - nbincdf(counts_2,pn2(1),pn2(2));

fdr1 = mafdr(pval1,'bhfdr',true);
fdr2 = mafdr(pval2,'bhfdr',true);

w1 = fdr1<.01; % logical vector indicating significant windows in HCT116-1
w2 = fdr2<.01; % logical vector indicating significant windows in HCT116-2
w12 = w1 & w2; % logical vector indicating significant windows in both replicates

Number_of_sig_windows_HCT116_1 = sum(w1)
Number_of_sig_windows_HCT116_2 = sum(w2)
Number_of_sig_windows_HCT116 = sum(w12)
end
