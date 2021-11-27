function [out] = classify_proj(lab,dat,nfolds,roctype,mfeats,dat_gen)
%  CLASSIFY_PROJ  Decode discrete labels using simple linear classifier based on
%                 multivariate normal model.
%
%  Usage: [out] = CLASSIFY_PROJ(lab,dat,nfolds,roctype,mfeats,dat_gen)
%
%  where lab is the set of discrete labels to be decoded
%        dat is the multivariate data array
%        nfolds is the number of cross-validated training folds
%        roctype is the ROC type flag ('onevsall' or 'onevsone')
%        mfeats is the number of discriminant data features to use (optional)
%        dat_gen is the generalization data array (optional)
%
%  Set mfeats = [] to use all data features for classification (default).
%
%  Valentin Wyart <valentin.wyart@inserm.fr>

% check input
if nargin < 6
    dat_gen = [];
end
if nargin < 5
    mfeats = [];
end
if nargin < 4
    roctype = 'onevsall';
end
if nargin < 3
    error('missing input!');
end

% preprocess category labels
lab = lab(:); % columnize vector of category labels
[lab_lst,~,lab] = unique(lab); % re-label using positive integers
ncat = length(lab_lst); % number of categories

% check data format
nevnts = size(dat,1); % number of epochs
nfeats = size(dat,2); % number of data features

% number of selected data features
if isempty(mfeats)
    mfeats = nfeats;
end

% build lookup table
ilab = sub2ind([nevnts,ncat],(1:nevnts)',lab);

if mfeats < nfeats
    kfeat = nan(nfolds,mfeats);
end

if ~isempty(dat_gen)
    ngen = size(dat_gen,3);
    xgen = nan(nevnts,ngen);
end

% compute category-wise Q-values
qcat = nan(nevnts,ncat);
xprj = nan(nevnts,1);

for ifold = 1:nfolds
    
    % split into training/test sets
    itst = mod((1:nevnts)',nfolds) == mod(ifold,nfolds);
    itrn = ~itst; 
    
    % center and scale data w.r.t. training set
    dat = bsxfun(@minus,dat,mean(dat(itrn,:),1));
    dat = bsxfun(@rdivide,dat,std(dat(itrn,:),[],1));
    
    if ~isempty(dat_gen)
        dat_gen = bsxfun(@minus,dat_gen,mean(dat_gen(itrn,:,:),1));
        dat_gen = bsxfun(@rdivide,dat_gen,std(dat_gen(itrn,:,:),[],1));
    end
        
    if nfeats > 1 && mfeats > 1 % multivariate data
        % select most discriminant data features w.r.t. training set
        if mfeats < nfeats
            pfeat = [];
            for ifeat = 1:nfeats
                [~,pfeat(ifeat)] = get_aroc(lab(itrn),dat(itrn,ifeat));
            end
            [~,ifeat] = sort(pfeat);
            ifeat = ifeat(1:mfeats);
            kfeat(ifold,:) = ifeat;
        else
            ifeat = 1:nfeats;
        end
        % compute category averages
        xdir_trn = nan(mfeats,ncat-1);
        xdir_tst = nan(mfeats,ncat-1);
        for icat = 1:ncat-1
            xdir_trn(:,icat) = mean(dat(itrn & lab == icat,ifeat),1);
            xdir_tst(:,icat) = mean(dat(itst & lab == icat,ifeat),1);
        end
        % normalize category averages to unit Euclidean norm
        xdir_trn = bsxfun(@rdivide,xdir_trn,sqrt(sum(xdir_trn.^2,1)));
        xdir_tst = bsxfun(@rdivide,xdir_tst,sqrt(sum(xdir_tst.^2,1)));
        % project data onto category averages
        xcat = dat(:,ifeat)*xdir_trn;
        
        if ~isempty(dat_gen)
            for igen = 1:ngen
                xgen(itst,igen) = dat_gen(itst,ifeat,igen)*xdir_trn;
            end
        end
        
    else % univariate data
        xcat = dat;
    end
    
    xprj(itst) = xcat(itst);
    
    % compute category-wise Q-values
    if ncat == 2 | nfeats == 1 % univariate decision space
        for icat = 1:ncat
            xavg = mean(xcat(itrn & lab == icat));
            xstd = std(xcat(itrn & lab == icat));
            qcat(itst,icat) = normllh(xcat(itst),xavg,xstd);
        end
    else % multivariate decision space
        for icat = 1:ncat
            xavg = mean(xcat(itrn & lab == icat,:));
            xcov = cov(xcat(itrn & lab == icat,:));
            qcat(itst,icat) = mvnllh(xcat(itst,:),xavg,xcov);
        end
    end
    
end

% convert Q-values to probabilities
qcat = bsxfun(@minus,qcat,qcat(ilab));
pcat = bsxfun(@rdivide,exp(qcat),sum(exp(qcat),2));
ptru = pcat(ilab);

% compute classifier sensitivity
proc = get_proc(lab,qcat,roctype);

% compute classifier accuracy
[~,hat] = max(pcat,[],2);
acc = hat == lab;
pacc = mean(acc);

% build output structure
out         = [];
out.nfolds  = nfolds; % number of cross-validation folds
out.roctype = roctype; % ROC type
out.nfeats  = nfeats; % total number of features
out.mfeats  = mfeats; % number of selected features
out.lab     = lab; % true category label
out.hat     = hat; % decoded category label
out.acc     = acc; % correct decoding?
out.xprj    = xprj; % projected discriminant component
out.qcat    = qcat; % category-wise Q-values
out.pcat    = pcat; % category-wise probabilities
out.ptru    = ptru; % probability of true category
out.proc    = proc; % sensitivity
out.pacc    = pacc; % accuracy

if mfeats < nfeats
    pfeat = nan(1,nfeats);
    for ifeat = 1:nfeats
        pfeat(ifeat) = nnz(kfeat(:) == ifeat)/nfolds;
    end
    out.pfeat = pfeat;
end

if ~isempty(dat_gen)
    out.xgen = xgen;
end

end