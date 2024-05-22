function [invOmegaBar] = inv_OmegaBar(G)
    [P, K] = size(G);
    %     sizeG = size(G);
    D = eye(size(G,1));
    V = G*G'+D;
    V1 = diag(diag(V).^(-1/2));
    G1 = V1*G;
    D1 = V1 * D * V1;
    invD1 = diag(1./diag(D1));
    invOmegaBar = invD1 - invD1*G1*(eye(K) + G1'*invD1*G1)^(-1)*G1'*invD1;
end