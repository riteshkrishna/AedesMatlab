%% Script to fit mixture model to Plutella Female data.
function [params_est,ressquared] = Plut_fit_Model_optimize(data_for_trimodal_fit)

x = data_for_trimodal_fit;

total_datapoints = numel(x);

% Bin the data
n = numel(unique(x));
xx = sort(unique(x));

databin = zeros(n,1);
for i=1:n
    databin(i) = numel(find(xx(i) == x));
end
databin = databin/total_datapoints;

figure; plot(xx,databin,'.-')

% Fit three independent distributions based on the observations from data

% Exponential dist
f1 = exppdf(xx,1/2.35) % Funny, mu in R, should be 1/mu in Matlab
f1 = (f1 * 0.4) / sum(f1) 
%figure; plot(xx,f1)

% Normal dist - 1
f2 = normpdf(xx,2.1,0.18)
f2 = (f2 * 0.3) / sum(f2)
%figure; plot(xx,f2)

% Normal dist -2 
f3 = normpdf(xx,2.42,0.07)
f3 = (f3 * 0.3) / sum(f3)
%figure; plot(xx,f3)

% All the above together
mf = f1 + f2 + f3;
% figure;
% plot(xx,databin);
% hold on
% plot(xx,mf,'rs-')
% hold off

% Using an optimizer to fit the curve
% bounds = [p1 exp_mu p2 norm_mu1 norm_sig1 (1 - p1 + p2) norm_mu2 norm_sig2]
lower_bound = [0.35 2.2 0.25 2.05 0.16 2.40 0.05]
upper_bound = [0.45 2.4 0.35 2.15 0.20 2.45 0.10]
params = (lower_bound + upper_bound) / 2
start = params;

% Define the function as a mixture of exp + normal + normal distribution
F_mixture = @(params, xin) params(1) * exppdf(xin,1/params(2))/sum(exppdf(xin,1/params(2))) ...
                         + params(3) * normpdf(xin,params(4),params(5))/sum(normpdf(xin,params(4),params(5)))...
                         + (1 - params(1) - params(3)) * normpdf(xin,params(6),params(7))/sum(normpdf(xin,params(6),params(7)));

F_sumsquares = @(params) sum((databin - F_mixture(params, xx)).^2)

opts = optimoptions('fminunc','Algorithm','quasi-newton','TolFun',1e-08,'Display','iter');
[params_est,ressquared,eflag,outputu] = fminunc(F_sumsquares,start,opts)

% Check the data with the estimations obtained
f1_est = exppdf(xx,1/params_est(2)) 
f1_est = (f1_est * params_est(1)) / sum(f1_est) 
% Normal dist - 1
f2_est = normpdf(xx,params_est(4),params_est(5))
f2_est = (f2_est * params_est(3)) / sum(f2_est)
% Normal dist -2 
f3_est = normpdf(xx,params_est(6),params_est(7))
f3_est = (f3_est * (1 - params_est(1) - params_est(3))) / sum(f3_est)

mf_est = f1_est + f2_est + f3_est;
figure;
plot(xx,databin,'*-');
hold on
plot(xx,mf_est,'rs-')
hold off

