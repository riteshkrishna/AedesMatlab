% data = female_depth_filtered;
% 
% % E-step
% mu = [2.1 2.4];
% sig = [0.2 0.4];
% pr = @(x) exp(-0.5*((x-mu)./sig).^2)./sig;
% prn = @(x) pr(x)./sum(pr(x));
% 
% prns = zeros([numel(data),2]);
% for j=1:numel(data) 
%     prns(j,:)=prn(data(j)); 
% end;
% prns(100:110,:)
% 
% % M step 
% mu = sum(prns.*repmat(data,[1,2]), 1) ./ sum(prns,1) 
% xmmu = repmat(data,[1,2]) - repmat(mu,[numel(data),1]); 
% sig = sqrt(sum(prns .* xmmu.^2, 1) ./ sum(prns,1))
% pop = sum(prns,1)/numel(data)
% 
% % Test
% %mu = [randsample(data,1) randsample(data,1)]; 
% %sig = [.3 .3];

data = female_depth_filtered;
mu = [2.1 2.4];
sig = [0.2 0.4];
for jj=1:50,
    pr = @(x) exp(-0.5*((x-mu)./sig).^2)./(2.506*sig);
    prn = @(x) pr(x)./sum(pr(x));
    for j=1:numel(data); 
        prns(j,:)=prn(data(j)); 
    end;
    
    mu = sum(prns.*repmat(data,[1,2]), 1) ./ sum(prns,1); 
    xmmu = repmat(data,[1,2]) - repmat(mu,[numel(data),1]); 
    sig = sqrt(sum(prns .* xmmu.^2, 1) ./ sum(prns,1)); 
    pop = sum(prns,1)/numel(data);
    
    thefunc = @(x) sum(pop.*pr(x),2);
    
    x = 1:.01:4;
    f = arrayfun(thefunc,x);
    plot(x,f,'b');
    hold on; 
end;

mu
sig

[f x] = ksdensity(data);
plot(x,f,'r')
hold off;