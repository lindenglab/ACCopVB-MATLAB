clc;clear;
addpath('Distribution/')


n = 1040;
G = [1, -1, -10;
    0, 2, 10]';
K = size(G,2);
OmegaBar_true = G2OmegaBar(G);
alpha_true = [-2,-2,3]';
delta_true = alpha2delta(alpha_true, OmegaBar_true);
nu_true = 4;

X = mvstrnd(alpha_true, OmegaBar_true, nu_true, n);

d = size(X,2);
U = zeros(n,d);
for j = 1:d
    U(:,j) = stcdf(X(:,j),delta_true(j),nu_true);
end

filename = sprintf("Data/Simulation_Data_3D.mat");
parsave(filename, X, U, G, OmegaBar_true, alpha_true, delta_true, nu_true);
