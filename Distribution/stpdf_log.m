function [pdf] = stpdf_log(X, alpha, nu)

Q = alpha .* X' .* sqrt((nu+1)./(nu+X'.^2));

pdf = log(tpdf(X', nu)) + log(tcdf(Q, nu+1));


end