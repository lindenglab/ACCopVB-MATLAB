function [output] = gradTG(C, V, V1, G)
    [P, K] = size(G);
    T1 = V*V1*C;
    T2 = C*V1*V;
    T3 = V1*C*V1;
    V2 = diag(V).^(-3/2).*eye(P);
    output = -1/2*(V2.*(T1 + T1' + T2 + T2'))*G + (T3 + T3')*G;
%     output = -V2.*(T1 + T1')*G + 2*T3*G;
end