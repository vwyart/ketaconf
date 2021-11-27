function [proc,out] = get_proc(lab,qcat,roctype)
%  GET_PROC  Get multinomial ROC sensitivity to categorical labels
%
%  Usage:
%    >> [proc,out] = GET_PROC(lab,qcat,roctype)
%
%  Input arguments:
%    - lab  = category labels
%    - qcat = multivariate data
%    - roctype = multinomial sensitivity type ('onevsall', 'onevsone')
%
%  Valentin Wyart <valentin.wyart@inserm.fr>

% check input
if nargin < 3
    roctype = 'onevsall';
end
if nargin < 2
    error('missing input!');
end

% preprocess category labels
lab = lab(:); % columnize category labels
[lab_lst,~,lab] = unique(lab); % re-label using positive integers
ncat = length(lab_lst); % number of categories

% check data format
if size(qcat,1) ~= length(lab) || size(qcat,2) ~= ncat
    error('invalid data format!');
end

% compute multinomial sensitivity
switch roctype
    case 'onevsall'
        lab_cmb = (1:ncat)';
        proc_cmb = nan(ncat,1);
        for icat = 1:ncat
            pcat = exp(qcat(:,icat))./sum(exp(qcat),2);
            proc_cmb(icat) = get_aroc(lab == icat,pcat);
        end
    case 'onevsone'
        lab_cmb = combnk(1:ncat,2);
        ncmb = size(lab_cmb,1);
        proc_cmb = nan(ncmb,1);
        for icmb = 1:ncmb
            isub = ismember(lab,lab_cmb(icmb,:));
            icat = lab_cmb(icmb,1);
            ialt = lab_cmb(icmb,2);
            pcat_cmb = exp(qcat(isub,icat))./sum(exp(qcat(isub,[icat,ialt])),2);
            proc_cmb(icmb) = get_aroc(lab(isub) == icat,pcat_cmb);
        end
    otherwise
        error('invalid input argument!');
end
proc = mean(proc_cmb);

if nargout > 1
    % create output structure
    out          = [];
    out.roctype  = roctype; % ROC type
    out.lab_lst  = lab_lst; % lookup table for category labels
    out.lab_cmb  = lab_cmb; % category label pairs
    out.proc     = proc; % multinomial sensitivity
    out.proc_cmb = proc_cmb; % pair-wise sensitivities
end

end