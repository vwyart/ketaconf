function [noccurences] = HasDistantValues(x,ndistant,xexcluded)
%  [noccurences] = HasDistantValues(x,[ndistant],[xexcluded])

if nargin < 3, xexcluded = []; end
if nargin < 2, ndistant = 2; end
if nargin < 1, error('Not enough input arguments'); end

if isempty(x)
    noccurences = 0;
    return
end

xval = setdiff(unique(x),xexcluded);

noccurences = 0;
for i = 1:length(xval)
    noccurences = noccurences+HasConsecutiveValues(x ~= xval(i),ndistant);
end

end