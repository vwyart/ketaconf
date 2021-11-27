function [out] = decode_pinv_tgen(C,D,nfold,sdtrn,Dgen)
%  DECODE_PINV_TGEN  Decode continuous variable from multivariate data using
%                    pseudoinverse (i.e., an inverted encoding model) with
%                    temporal generalization
%
%  Usage: [out] = DECODE_PINV_TGEN(C,D,nfold,sdtrn,Dgen)
%
%  where C is the continuous variable to be decoded
%        D is the multivariate data array
%        nfold defines the number of interleaved cross-validation folds
%        sdtrn defines whether to standardize the data on each training fold
%        Dgen is the generalization data array
%
%  The function outputs out, which contains Chat = predicted (cross-validated)
%  values for the continuous variable and Chat_gen = predicted values for the
%  continuous variable based on the generalization data array.
%
%  Valentin Wyart <valentin.wyart@inserm.fr>

% check input arguments
if nargin < 4
    sdtrn = true;
end
if nargin < 3
    error('missing input arguments!');
end
if size(D,1) ~= size(C,1)
    error('mismatching input sizes!');
end

n     = size(C,1); % number of observations
nregr = size(C,2); % number of regressors
ntgen = size(Dgen,3);

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
Chat = nan(n,nregr+1);
Chat_gen = nan(n,nregr+1,ntgen);
for ifold = 1:nfold
    if sdtrn
        % standardize data w.r.t. training set
        D = bsxfun(@minus,D,mean(D(itrn{ifold},:),1));
        D = bsxfun(@rdivide,D,std(D(itrn{ifold},:),[],1));
        Dgen = bsxfun(@minus,Dgen,mean(Dgen(itrn{ifold},:,:),1));
        Dgen = bsxfun(@rdivide,Dgen,std(Dgen(itrn{ifold},:,:),[],1));
    end
    % compute weights on training set
    Dtrn = D(itrn{ifold},:)';
    Ctrn = C(itrn{ifold},:)';
    Wenc = Dtrn*Ctrn'*pinv(Ctrn*Ctrn'); % encoding weights
    Wdec = pinv(Wenc'*Wenc)*Wenc'; % decoding weights
    % compute predictions on test set
    Chat(itst{ifold},:) = (Wdec*D(itst{ifold},:)')';
    for itgen = 1:ntgen
        Chat_gen(itst{ifold},:,itgen) = (Wdec*Dgen(itst{ifold},:,itgen)')';
    end
end

% create output structure
out      = [];
out.Ctru = C(:,1:end-1);
out.Chat = Chat(:,1:end-1);
out.Chat_gen = squeeze(Chat_gen(:,1:end-1,:));

end