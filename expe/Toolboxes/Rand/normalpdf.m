function [p] = normalpdf(x,mu,sigma,ispdf)

if nargin < 4
    ispdf = true;
end
if nargin < 3
    sigma = 1;
end
if nargin < 2
    mu = 0;
end
if nargin < 1
    error('Wrong input argument list.');
end

p = exp(-0.5*((x-mu)./sigma).^2);
if ispdf
    p = p./(sqrt(2*pi).*sigma);
end

end