function [dep] = qtdep_emp( CopulaData, u )
% EMPQTDEP Calculates the quantile dependence of empirical copula function
% 
% Inputs
%  U - Copula data (uniform) which should be a n x D input, where n is the
%      number of samples and D is the dimensionality
%  q - quantile value
% Outputs
%  LL - Lower Left dependence
%  UR - Upper Right dependence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(u);
[LL, UR, LR, UL] = deal(zeros(n,1));
for i = 1:n
    idL = CopulaData(:,1) < u(i);
    idR = CopulaData(:,1) > 1-u(i);
    

    LL(i) = mean( CopulaData(idL,2) < u(i)  ); % LL: P( U2<q   | U1<q )
    UR(i) = mean( CopulaData(idR,2) > 1-u(i)); % UR: P( U2>1-q | U1>1-q )
    LR(i) = mean( CopulaData(idR,2) < u(i)  ); % LR: P( U2<q   | U1>1-q )
    UL(i) = mean( CopulaData(idL,2) > 1-u(i)); % UL: P( U2>1-q | U1<q )
end
dep.LL = LL;
dep.UR = UR;
dep.LR = LR;
dep.UL = UL;
end