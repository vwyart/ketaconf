function [d,k] = vmpar(x,dim,n)
%  VMPAR  Get von Mises parameters from array of angles
%
%  Usage: [d,k] = VMPAR(x,dim,n)
%
%  where x   - array of angles (pi-periodic)
%        dim - dimension of interest
%        n   - number of approximation steps (exact solution if < 0)
%        d   - von Mises direction
%        k   - von Mises concentration
%
%  Valentin Wyart <valentin.wyart@ens.fr>

% check input arguments
if nargin < 3
    n = [];
end
if nargin < 2
    dim = 1;
end
if nargin < 1
    error('missing input array of angles!');
end

% reshape input array if necessary
if isscalar(x)
    error('input array must not be scalar!');
elseif isvector(x)
    x = x(:);
end

% define function handle
modcc = @(x)mod(x+pi/2,pi)-pi/2; % pi-periodic modulus

% average evidence on complex plane
z = mean(exp(1i*2*x),dim);

d = modcc(angle(z)/2); % evidence direction
r = abs(z); % evidence norm
k = vmr2k(r,n); % evidence concentration

end