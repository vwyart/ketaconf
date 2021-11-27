function [results] = GetStaircaseThreshold(results,iout)

if nargin < 2
    error('missing input argument(s).');
end

iout = min(iout,length(results.istp)-1);
i1st = results.istp(iout)+1;

xthres = 10^mean(log10(results.x(results.istp(iout+1:end))));
pthres = [];

rthres = results.r(i1st:end);
if all(ismember(rthres,[0 1]))
    pthres = mean(rthres);
else
    error('unable to get staircase threshold!');
end

results.xthres = xthres;
results.pthres = pthres;

end