function x = tnmean(a,b,m,s)
%  TNMEAN  Mean of truncated normal distribution
%
%  Usage: x = tnmean(a,b) for (m = 0, s = 1)
%         x = tnmean(a,b,m,s)
%  where a <= b (lower and upper truncation bounds, respectively)
%
%  Valentin Wyart <valentin.wyart@ens.fr>
%  adapted from github.com/cossio/TruncatedNormal.jl
if ~ismember(nargin,[2,4])
    error('invalid argument list!');
end
if nargin == 2
    if a < b
        x = sqrt(2/pi)*F1(a/sqrt(2),b/sqrt(2));
    elseif a == b
        x = a;
    else
        error('a should be less or equal to b!');
    end
elseif nargin == 4
    x = m+tnmean((a-m)/s,(b-m)/s)*s;
end
end

function [z] = F1(x,y,thresh)
% returns (exp(-x^2)-exp(-y^2))/(erf(y)-erf(x)) without catastrophic cancellation
if nargin < 3
    thresh = 1e-7;
end
if ~isscalar(x) || ~isscalar(y)
    error('x and y should be scalar!');
end
if abs(x) > abs(y)
    z = F1(y,x,thresh);
    return
elseif isinf(y)
    z = sign(y)/(erfcx(sign(y)*x));
    return
elseif abs(x-y) <= thresh
    e = y-x;
    p = sqrt(pi);
    z = p*x+(p/2+(-p*x/6+(-p/12+x*(p/90+(p*x^2)/90)*e)*e)*e)*e;
    return
end
d = exp(x^2-y^2);
if max(x,y) < 0
    z = (1-d)/(d*erfcx(-y)-erfcx(-x));
elseif min(x,y) > 0 || isinf(y)
    z = (1-d)/(erfcx(x)-d*erfcx(y));
else
    z = exp(-x^2)*(1-d)/(erf(y)-erf(x));
end
end