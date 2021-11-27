function [r] = vmk2r(k)
%  VMK2R  Get von Mises coherence from concentration
%
%  Usage: [r] = VMK2R(k)
%
%  where k - von Mises concentration
%        r - von Mises coherence (vector norm)
%
%  Valentin Wyart <valentin.wyart@ens.fr>

if nargin < 1
    error('missing input concentration!');
end

r = besseli(1,k)./besseli(0,k);

end