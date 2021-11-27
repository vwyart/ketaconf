function [out] = decode_pinv(C,D,nfold,sdtrn,Winp)
%  DECODE_PINV  Decode continuous variable from multivariate data using
%               pseudoinverse (i.e., an inverted encoding model)
%
%  Usage: [out] = DECODE_PINV(C,D,nfold,sdtrn,Winp)
%
%  where C is the continuous variable to be decoded
%        D is the multivariate data array
%        nfold defines the number of interleaved cross-validation folds
%        sdtrn defines whether to standardize the data on each training fold
%        Winp allows optionally to provide pre-defined decoding weights
%
%  The function outputs out, which contains Chat = predicted (cross-validated)
%  values for the continuous variable, and W = decoding weights.
%
%  Valentin Wyart <valentin.wyart@inserm.fr>

% check input arguments
if nargin < 5
    Winp = {};
end
if nargin < 4
    sdtrn = true;
end
if nargin < 3
    error('missing input arguments!');
end
if size(D,1) ~= size(C,1)
    error('mismatching input sizes!');
end
if ~isempty(Winp) && numel(Winp) ~= nfold
    error('mismatching decoding weights!');
end

n = size(C,1); % number of observations

% define training/test sets
if iscell(nfold)
    itrn = nfold(:,1);
    itst = nfold(:,2);
    nfold = size(nfold,1);
else
    itrn = cell(nfold,1);
    itst = cell(nfold,1);
    for ifold = 1:nfold
        itst{ifold} = find(mod((1:n)',nfold) == ifold-1);
        itrn{ifold} = find(mod((1:n)',nfold) ~= ifold-1);
    end
end
% ignore trials where the regressor is undefined
for ifold = 1:nfold
    itst{ifold} = intersect(itst{ifold},find(all(~isnan(C),2)));
    itrn{ifold} = intersect(itrn{ifold},find(all(~isnan(C),2)));
end

% add intercept term to design matrix
C = cat(2,C,ones(n,1));

% compute predictions
Chat = nan(size(C));
Wout = cell(nfold,1);
for ifold = 1:nfold
    if sdtrn
        % standardize data w.r.t. training set
        D = bsxfun(@minus,D,mean(D(itrn{ifold},:),1));
        D = bsxfun(@rdivide,D,std(D(itrn{ifold},:),[],1));
    end
    % compute weights on training set
    Dtrn = D(itrn{ifold},:)';
    Ctrn = C(itrn{ifold},:)';
    if isempty(Winp)
        Wenc = Dtrn*Ctrn'*pinv(Ctrn*Ctrn'); % encoding weights
        Wdec = pinv(Wenc'*Wenc)*Wenc'; % decoding weights
        Wout{ifold} = Wdec;
    else
        Wdec = Winp{ifold};
    end
    % compute predictions on test set
    Dtst = D(itst{ifold},:)';
    Ctst = Wdec*Dtst;
    Chat(itst{ifold},:) = Ctst';
end

% create output structure
out      = [];
out.Ctru = C(:,1:end-1);
out.Chat = Chat(:,1:end-1);
if isempty(Winp)
    out.W = Wout;
end

end