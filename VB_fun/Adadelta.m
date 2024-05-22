function [newValue, state] = Adadelta(gradients, state)

if nargin == 1
    state = struct;
end
if ~isfield(state, 'epsilon')
    state.epsilon = 1e-6;
end
if ~isfield(state, 'rho')
    state.rho = 0.95;
end
if ~isfield(state, 'gradient')
    state.accGrad = zeros(size(gradients));
end
if ~isfield(state, 'value')
    state.accVal = zeros(size(gradients));
end

% Accumulate Gradients
state.accGrad = state.rho * state.accGrad + (1 - state.rho) * gradients.^2;
% Update Values
newValue = gradients .* ( sqrt(state.accVal + state.epsilon) ./ sqrt(state.accGrad + state.epsilon) );
% Accumulate Update Values
state.accVal = state.rho * state.accVal + (1 - state.rho) * newValue.^2;

end
