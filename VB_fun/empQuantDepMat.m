function [LL_mat, UR_mat, LR_mat, UL_mat] = empQuantDepMat(U, q)
d = size(U,2);
[LL_mat, UR_mat, LR_mat, UL_mat] = deal(zeros(d,d));
for a = 1:d
    for b = 1:d
        %         if a ~= b
        if b < a
            [LL_mat(a,b), UR_mat(a,b), LR_mat(a,b), UL_mat(a,b)] = empqtdep(U(:,[a,b]), q );
        end
    end
end
end