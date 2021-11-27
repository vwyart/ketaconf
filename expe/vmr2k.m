function [k] = vmr2k(r,n)
%  VMR2K  Get von Mises concentratrion from coherence
%
%  Usage: [k] = VMR2K(r,n)
%
%  where r - von Mises coherence (vector norm)
%        n - number of approximation steps (exact solution if < 0)
%        k - von Mises concentration
%
%  Valentin Wyart <valentin.wyart@ens.fr>

% check input arguments
if nargin < 2 || isempty(n)
    n = 0;
end
if nargin < 1
    error('missing input vector norm!');
end

% define function handle 
A = @(x)besseli(1,x)./besseli(0,x); % Bessel ratio

if n < 0
    % get exact solution (slower)
    k = nan(size(r));
    for i = 1:numel(r)
        k(i) = fzero(@(x)A(x)-r(i),0);
    end
else
    % get approximate solution (faster)
    % 1. get initial estimate
    k = r.*(2-r.^2)./(1-r.^2);
    % 2. iterate Newton approximation
    for i = 1:n
        k = k-(A(k)-r)./(1-A(k).^2-A(k)./k);
    end
end