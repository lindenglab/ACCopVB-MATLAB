function [ x ] = cnormrnd( mu, sig, llim, ulim)
% MAIN FUNCTION:
% cnormrnd: generate from a constrained normal distribution
% 
% INPUT:
% mu: mean list (1 x n)
% sig: standard deviation list (1 x n)
% llim: lower bound
% ulim: upper bound
% 
% OUTPUT:
% x: random sample list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Ole Maneesoonthorn || updated by Lin Deng with vector input 

if(ulim<llim)
 error('Fatal Error are cigrnd: upper limits <= lower limits');
end

if(isinf(llim))
    fa=0; %lower limit set to -inf
else
    fa=normcdf(llim,mu,sig); %lower limit set to F(llim), F(.) cdf of N(mu,sig)
end

if(isinf(ulim))
    fb=1; %upper limit set to +inf
else
    fb=normcdf(ulim,mu,sig);  %upper limit set to F(ulim), F(.) cdf of N(mu,sig)
end

n = size(mu, 2);
% uni=rand.*(fb-fa)+fa;

randList = rand(size(mu));
uni=randList.*(fb-fa)+fa;

x=norminv(uni,mu,sig);

% if sum(isinf(x)) > 0
%     inf_list = isinf(x);
%     x(inf_list) = norminv_ext(uni(inf_list), mu(inf_list), sig(inf_list));
% end
end

