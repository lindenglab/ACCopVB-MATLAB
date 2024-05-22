function [delta] = alpha2delta(alphavec, OmegaBar)
    delta = (1 + alphavec' * OmegaBar * alphavec)^(-1/2) * (OmegaBar * alphavec);
end

