function [output] = gradTAlpha(Ci, C1, alpha, OmegaBar)
    output = C1^(-2)*OmegaBar*alpha*alpha'*OmegaBar*(Ci+Ci')*OmegaBar*alpha...
        - C1^(-1)*OmegaBar*(Ci+Ci')*OmegaBar*alpha;
end