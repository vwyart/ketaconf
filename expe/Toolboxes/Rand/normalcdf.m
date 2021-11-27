function [p] = normalcdf(x,mu,sigma)

if nargin < 3
    sigma = 1;
end
if nargin < 2
    mu = 0;
end
if nargin < 1
    error('Wrong input argument list.');
end

p = 0.5*erfc(-(x-mu)./(sqrt(2).*sigma));

end