function [out] = fit_model_resp(cfg)
%  FIT_MODEL_RESP  Fit inference model to binary responses
%
%  Usage: [out] = FIT_MODEL_RESP(cfg)
%
%  where cfg is the configuration structure
%        out is the output structure
%
%  The configuration structure should contain the following fields:
%    * stype  = session type => 1:placebo or 2:ketamine
%    * seqang = sequence-wise angles (rad)
%    * catdir = categorization axis (rad)
%    * resp   = response => 1 or 2
%    * parmod = list of parameters modulated by session type
%
%  The model parameters are the following:
%    * ec_inf = inference compressive non-linearity (0=unbiased)
%    * el_inf = inference leaky accumulation (0=unbiased)
%    * sd_inf = inference noise s.d.
%    * criter = decision criterion
%
%  Any combination of parameters can be provided as fields of the configuration
%  as fixed parameters that will not be fitted.
%
%  The function uses the interior-point algorithm of the fmincon function to
%  maximize the log-likelihood of the model.
%
%  Valentin Wyart <valentin.wyart@inserm.fr>

% check configuration parameters
if ~all(isfield(cfg,{'stype','seqang','catdir','resp'}))
    error('Incomplete configuration structure!');
end
if ~isfield(cfg,'parmod')
    cfg.parmod = {};
end
if ~isfield(cfg,'only_incr')
    cfg.only_incr = false;
end

% get configuration parameters
stype  = cfg.stype(:); % session type => 1:placebo or 2:ketamine
seqang = cfg.seqang(:); % sequence-wise angles (rad)
catdir = cfg.catdir(:); % categorization axis (rad)
resp   = cfg.resp(:); % response => 1 or 2
parmod = cfg.parmod; % list of parameters modulated by session type

% set useful function handles
col = @(x)x(:);
row = @(x)x(:)';
sub = @(x,k)x(min(k,numel(x)));
logit = @(x)log(x./(1-x));

% get sequence lengths
seqlen = cellfun(@length,seqang);
[seqlen_lst,~,seqlen_ind] = unique(seqlen);

nseq = length(seqlen); % number of sequences
maxlen = max(seqlen); % maximum sequence length

% create matrix of sequence angles
seqang = cellfun(@(x,n)cat(2,fliplr(row(x(1:n))),nan(1,maxlen-n)), ...
    seqang,num2cell(seqlen),'UniformOutput',false);
seqang = cat(1,seqang{:});
% express sequence angles as tilt from categorization axis
seqang = mod(bsxfun(@minus,seqang,catdir),pi);

% set response lookup indices
iobs = find(ismember(resp,[1,2])); % exclude trials with missing response
if isfield(cfg,'ifilt')
    iobs = intersect(iobs,cfg.ifilt);
end
nobs = length(iobs);
ipost = sub2ind([nseq,2],iobs,resp(iobs));

% define model parameters
npar = 0;
pnam = {}; % parameter name
psiz = []; % parameter size
pini = []; % parameter initialization value
pmin = []; % parameter minimum value
pmax = []; % parameter maximum value
% inference compressive non-linearity
npar = npar+1;
pnam{npar,1} = 'ec_inf';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini(npar,1) = 0;
pmin(npar,1) = -2;
pmax(npar,1) = +2;
if ismember('ec_inf',parmod) && isfield(cfg,'ec_inf') && cfg.only_incr
    % restrict to ketamine increase when placebo is provided
    pini(npar) = cfg.sd_inf(1);
    pmin(npar) = cfg.sd_inf(1);
end
% inference leaky accumulation
npar = npar+1;
pnam{npar,1} = 'el_inf';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini(npar,1) = 0;
pmin(npar,1) = -1;
pmax(npar,1) = 10;
if ismember('el_inf',parmod) && isfield(cfg,'el_inf') && cfg.only_incr
    % restrict to ketamine increase when placebo is provided
    pini(npar) = cfg.sd_inf(1);
    pmin(npar) = cfg.sd_inf(1);
end
% inference noise
npar = npar+1;
pnam{npar,1} = 'sd_inf';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini(npar,1) = 1;
pmin(npar,1) = 0;
pmax(npar,1) = 10;
if ismember('sd_inf',parmod) && isfield(cfg,'sd_inf') && cfg.only_incr
    % restrict to ketamine increase when placebo is provided
    pini(npar) = cfg.sd_inf(1);
    pmin(npar) = cfg.sd_inf(1);
end
% decision criterion
npar = npar+1;
pnam{npar,1} = 'criter';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini(npar,1) = 0;
pmin(npar,1) = -5;
pmax(npar,1) = +5;

% define fixed parameters
pfix = cell(npar,1);
for i = 1:npar
    if isfield(cfg,pnam{i})
        pfix{i} = reshape(cfg.(pnam{i}),[psiz(i),1]);
    end
end

% define free parameters to fit
ifit = cell(npar,1);
pfit_nam = {};
pfit_ini = [];
pfit_min = [];
pfit_max = [];
n = 1;
for i = 1:npar
    if isempty(pfix{i})
        ifit{i} = n+[1:psiz(i)]-1;
        pfit_nam = cat(1,pfit_nam,pnam(i*ones(psiz(i),1)));
        pfit_ini = cat(1,pfit_ini,pini(i*ones(psiz(i),1)));
        pfit_min = cat(1,pfit_min,pmin(i*ones(psiz(i),1)));
        pfit_max = cat(1,pfit_max,pmax(i*ones(psiz(i),1)));
        n = n+psiz(i);
    elseif numel(pfix{i}) ~= psiz(i)
        error('wrong size for fixed parameter %s!',pnam{i});
    elseif any(isnan(pfix{i})) % partly fixed parameter
        nfit = nnz(isnan(pfix{i}));
        ifit{i} = n+[1:nfit]-1;
        pfit_nam = cat(1,pfit_nam,pnam(i*ones(nfit,1)));
        pfit_ini = cat(1,pfit_ini,pini(i*ones(nfit,1)));
        pfit_min = cat(1,pfit_min,pmin(i*ones(nfit,1)));
        pfit_max = cat(1,pfit_max,pmax(i*ones(nfit,1)));
        n = n+nfit;
    end
end
nfit = length(pfit_ini);

% fit free parameters
if nfit > 0
    A = [];
    b = [];
    if cfg.only_incr
        % constrain to placebo < ketamine
        tnam = {'ec_inf','el_inf','sd_inf'};
        for t = 1:3
            if ismember(tnam{t},parmod) && ~isfield(cfg,tnam{t})
                i = ismember(pfit_nam,tnam{t});
                At = zeros(1,nfit);
                At(i) = [+1,-1];
                A = [A;At];
                b = [b;0];
            end
        end
    end
    pval = fmincon(@fmin,pfit_ini,A,b,[],[],pfit_min,pfit_max,[], ...
        optimset('Display','notify','FunValCheck','on','Algorithm','interior-point','TolX',1e-20,'MaxFunEvals',1e6));
    [~,phat] = fmin(pval);
else
    phat = pfix;
end

% get best-fitting parameters
out = cell2struct(phat,pnam);

% save configuration structure
fnames = fieldnames(out);
fnames = cat(1,{'cfg'},fnames);
out.cfg = cfg;
out = orderfields(out,fnames);

% get quality of fit
out.nfit = nfit; % number of fitted parameters
out.nobs = nobs; % number of observations
out.iobs = iobs; % indices of observations
out.llh  = get_llh(phat{:}); % model log-likelihood
out.aic  = -2*out.llh+2*nfit+2*nfit*(nfit+1)./(nobs-nfit+1); % model AICc
out.bic  = -2*out.llh+nfit*log(nobs); % model BIC

% get posterior probabilities
out.ppost = get_ppost(phat{:});

% get accumulated evidence
[~,out.xevi,out.sevi,out.tevi] = get_ppost(phat{:});

    function [f,pfit] = fmin(p)
        % create list of model parameters
        pfit = cell(npar,1);
        for i = 1:npar
            if isempty(pfix{i}) % free parameter
                pfit{i} = p(ifit{i});
            elseif any(isnan(pfix{i})) % partly fixed parameter
                pfit{i} = pfix{i};
                pfit{i}(isnan(pfit{i})) = p(ifit{i});
            else % fully fixed parameter
                pfit{i} = pfix{i};
            end
        end
        % return negative model log-likelihoods
        f = -get_llh(pfit{:});
    end

    function [llh,ppost,llh_stype] = get_llh(varargin)
        % get model posterior probabilities
        ppost = get_ppost(varargin{:});
        % get model log-likelihood
        llh = sum(log(max(ppost(ipost),eps)));
        if nargout > 2
            llh_stype = [ ...
                sum(log(max(ppost(ipost(stype == 1)),eps))), ...
                sum(log(max(ppost(ipost(stype == 2)),eps)))];
        end 
    end

    function [ppost,xevi,sevi,tevi] = get_ppost(ec_inf,el_inf,sd_inf,criter)
        ppost_all = nan(nseq,2); % posterior probabilities
        xevi_all  = nan(nseq,1); % sequence evidence
        sevi_all  = nan(nseq,1); % internal noise
        tevi_all  = nan(nseq,1); % accumulated evidence conditioned on responses
        for s = 1:2
            % filter session-wise trials
            ifilt = find(stype == s);
            nfilt = length(ifilt);
            % get parameters
            ec_inf_sub = sub(ec_inf,s);
            el_inf_sub = sub(el_inf,s);
            sd_inf_sub = sub(sd_inf,s);
            criter_sub = sub(criter,s);
            % compute evidence compression
            ndev = 10; % number of developments of infinite series
            xang = col(seqang(ifilt,:));
            xang = xang+0.5/besseli(0,ec_inf_sub)* ...
                sum(bsxfun(@times,sin(bsxfun(@times,1:ndev,xang*4)),besseli(1:ndev,ec_inf_sub)./(1:ndev)),2);
            xang = reshape(xang,[nfilt,maxlen]);
            xevi = sin(2*xang);
            xevi(isnan(xevi)) = 0;
            % compute leaky evidence accumulation
            xevi = nansum(bsxfun(@times,xevi,exp(-el_inf_sub).^(0:maxlen-1)),2);
            sevi = sqrt(arrayfun(@(n)sum(exp(-el_inf_sub).^(2*(0:n-1))),seqlen_lst(:))*sd_inf_sub^2);
            sevi = sevi(seqlen_ind(ifilt));
            xevi_all(ifilt) = xevi;
            sevi_all(ifilt) = sevi;
            % truncate accumulated evidence wrt response
            tevi = nan(nfilt,1);
            for u = 1:nfilt
                if resp(ifilt(u)) == 1
                    tevi(u) = +tnmean(0,+inf,xevi(u)-criter_sub,sevi(u));
                else
                    tevi(u) = -tnmean(-inf,0,xevi(u)-criter_sub,sevi(u));
                end
            end
            fevi = logit(0.75)/mean(tevi); % normalization factor
            tevi = tevi*fevi; % normalize truncated evidence
            tevi_all(ifilt) = tevi;
            % compute response probabilities
            ppost = normcdf(xevi,criter_sub,sevi);
            ppost = [ppost,1-ppost];
            % append session-wise predictions
            ppost_all(ifilt,:) = ppost;
        end
        ppost = ppost_all;
        xevi  = xevi_all;
        sevi  = sevi_all;
        tevi  = tevi_all;
    end

end