function [x] = GetStaircaseVariable(staircase)

if nargin < 1
    error('missing input argument.');
end

x = staircase.x(staircase.i);

end