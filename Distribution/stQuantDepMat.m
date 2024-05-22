function [LL_mat, UR_mat, LR_mat, UL_mat] = stQuantDepMat(u, delta, OmegaBar, nu)
d = size(delta,1);
[LL_mat, UR_mat, LR_mat, UL_mat] = deal(zeros(size(OmegaBar)));
for a = 1:d
    for b = 1:d
%         if a ~= b
        if b < a    
            deltaA = delta([a,b]);
            rho = OmegaBar(a,b);
            OmegaBarA = [1, rho; rho, 1];
            [LL_mat(a,b), UR_mat(a,b), LR_mat(a,b), UL_mat(a,b)] = stquantaildep_AC(u, deltaA, OmegaBarA, nu);
        end
    end
end
end