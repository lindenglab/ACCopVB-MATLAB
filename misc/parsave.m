function parsave(varargin)
% PARSAVE: save mat file within parfor loop

savefile = varargin{1}; % first input argument
for i = 2:nargin
    savevar.(inputname(i)) = varargin{i}; % other input arguments
end
save(savefile,'-struct','savevar')

end