function [alpha] = delta2alpha(deltavec, OmegaBar)
    invOmegaBar = OmegaBar^(-1);
    alpha = (1 - deltavec' * invOmegaBar * deltavec)^(-1/2) * invOmegaBar * deltavec;
end