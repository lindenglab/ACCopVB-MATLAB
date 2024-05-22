function [LL, UR, LR, UL] = empqtdep( U, q )
% Calculates the quantile dependence of empirical copula function
% Inputs
%  U - Copula data (uniform) which should be a n x D input, where n is the
%      number of samples and D is the dimensionality
%  q - quantile value
% Outputs
%  LL - Lower Left dependence
%  UR - Upper Right dependence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_q = length(q);

[LL, UR, LR, UL] = deal(zeros(n_q,1));
for i = 1:n_q
    idL = U(:,1) < q(i);
    idR = U(:,1) > 1-q(i);
    
    LL(i) = mean( U(idL,2) < q(i)  ); % LL: P( U2<q   | U1<q )
    UR(i) = mean( U(idR,2) > 1-q(i)); % UR: P( U2>1-q | U1>1-q )
    LR(i) = mean( U(idR,2) < q(i)  ); % LR: P( U2<q   | U1>1-q )
    UL(i) = mean( U(idL,2) > 1-q(i)); % UL: P( U2>1-q | U1<q )
end
end