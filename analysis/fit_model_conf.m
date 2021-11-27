function [out] = fit_model_conf(cfg)
%  FIT_MODEL_CONF  Fit inference and confidence model to binary responses and
%                  opt-out judgments
%
%  Usage: [out] = FIT_MODEL_CONF(cfg)
%
%  where cfg is the configuration structure
%        out is the output structure
%
%  The configuration structure should contain the following fields:
%    * stype  = session type => 1:placebo or 2:ketamine
%    * seqang = sequence-wise angles (rad)
%    * seqcat = sequence category => 1:orange or 2:blue
%    * catdir = categorization axis (rad)
%    * lotval = lottery probability of success
%    * resp   = response => 1:orange or 2:blue
%    * conf   = confirmation response => 1:confirm 2:opt-out
%    * parmod = list of parameters modulated by session type
%
%  The model parameters are the following:
%    * ec_inf = inference compressive non-linearity (0=unbiased)
%    * el_inf = inference leaky accumulation (0=unbiased)
%    * sd_inf = inference noise s.d.
%    * criter = decision criterion
%    * eg_cnf = confirmation gain
%    * eb_cnf = confirmation bias
%    * sd_cnf = confirmation noise (metacognitive noise)
%
%  Any combination of parameters can be provided as fields of the configuration
%  as fixed parameters that will not be fitted.
%
%  The function uses the interior-point algorithm of the fmincon function to
%  maximize the log-likelihood of the model.
%
%  Valentin Wyart <valentin.wyart@inserm.fr>


% check configuration parameters
if ~all(isfield(cfg,{'stype','seqang','catdir','lotval','resp','conf'}))
    error('Incomplete configuration structure!');
end
if ~isfield(cfg,'parmod')
    % fix all parameters across sessions
    cfg.parmod = {};
end
if ~isfield(cfg,'norm_ref')
    % fixed normalization reference at 75% correct
    cfg.norm_ref = 'fixd';
end
if ~isfield(cfg,'norm_typ')
    % normalization wrt to probability correct
    cfg.norm_typ = 'p';
end
if ~isfield(cfg,'get_fwhm')
    % do not compute full width at half maximum for key parameters
    cfg.get_fwhm = false;
end

% check normalization parameters
if ~ismember(cfg.norm_ref,{'fixd','accu'}) % normalization reference
    error('Undefined normalization reference!');
end
if ~ismember(cfg.norm_typ,{'p','l'}) % normalization type
    error('Undefined normalization type!');
end
norm_ref = cfg.norm_ref;
norm_typ = cfg.norm_typ;

% get full width at half maximum for key parameters?
get_fwhm = cfg.get_fwhm;

% get configuration parameters
parmod = cfg.parmod; % list of parameters modulated by session type
stype  = cfg.stype(:); % session type => 1:placebo 2:ketamine
seqcat = cfg.seqcat(:); % sequence category => 1:orange 2:blue
seqang = cfg.seqang(:); % sequence angles (rad)
catdir = cfg.catdir(:); % categorization direction (rad)
lotval = cfg.lotval(:); % lottery probability of success
resp   = cfg.resp(:); % response => 1:orange 2:blue
conf   = cfg.conf(:); % confirmation response => 1:confirm 2:opt-out

% set normalization reference
switch norm_ref
    case 'fixd' % instructed = 75%
        p0 = mean(unique(lotval))*[1,1];
    case 'accu' % accuracy
        accu = resp == seqcat;
        p0 = [mean(accu(stype == 1)),mean(accu(stype == 2))];
    otherwise
        error('Undefined reference probability correct!');
end

% set useful function handles
col = @(x)x(:);
row = @(x)x(:)';
sub = @(x,k)x(min(k,numel(x)));
logit = @(p)log(p./(1-p));

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
nobs = length(iobs);
ipost = sub2ind([nseq,4],iobs,2*(resp(iobs)-1)+conf(iobs));

% define model parameters
npar = 0;
pnam = {}; % parameter name
psiz = []; % parameter size
pini = {}; % parameter initialization value
pmin = []; % parameter minimum value
pmax = []; % parameter maximum value
% inference compression
npar = npar+1;
pnam{npar,1} = 'ec_inf';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini{npar,1} = 0;
pmin(npar,1) = -2;
pmax(npar,1) = +2;
if isfield(cfg,'ini') && isfield(cfg.ini,pnam{end})
    % custom initialization
    pini{end} = cfg.ini.(pnam{end})(:);
end
% inference leak
npar = npar+1;
pnam{npar,1} = 'el_inf';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini{npar,1} = 0;
pmin(npar,1) = -1;
pmax(npar,1) = 10;
if isfield(cfg,'ini') && isfield(cfg.ini,pnam{end})
    % custom initialization
    pini{end} = cfg.ini.(pnam{end})(:);
end
% inference noise
npar = npar+1;
pnam{npar,1} = 'sd_inf';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini{npar,1} = 1;
pmin(npar,1) = 0;
pmax(npar,1) = 10;
if isfield(cfg,'ini') && isfield(cfg.ini,pnam{end})
    % custom initialization
    pini{end} = cfg.ini.(pnam{end})(:);
end
% decision criterion
npar = npar+1;
pnam{npar,1} = 'criter';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini{npar,1} = 0;
pmin(npar,1) = -5;
pmax(npar,1) = +5;
if isfield(cfg,'ini') && isfield(cfg.ini,pnam{end})
    % custom initialization
    pini{end} = cfg.ini.(pnam{end})(:);
end
% confirmation gain
npar = npar+1;
pnam{npar,1} = 'eg_cnf';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini{npar,1} = 0;
pmin(npar,1) = -3;
pmax(npar,1) = +3;
if isfield(cfg,'ini') && isfield(cfg.ini,pnam{end})
    % custom initialization
    pini{end} = cfg.ini.(pnam{end})(:);
end
% confirmation bias
npar = npar+1;
pnam{npar,1} = 'eb_cnf';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini{npar,1} = 0;
pmin(npar,1) = -3;
pmax(npar,1) = +3;
if isfield(cfg,'ini') && isfield(cfg.ini,pnam{end})
    % custom initialization
    pini{end} = cfg.ini.(pnam{end})(:);
end
% confirmation noise
npar = npar+1;
pnam{npar,1} = 'sd_cnf';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini{npar,1} = 0;
pmin(npar,1) = 0;
pmax(npar,1) = 100;
if isfield(cfg,'ini') && isfield(cfg.ini,pnam{end})
    % custom initialization
    pini{end} = cfg.ini.(pnam{end})(:);
end

% define fixed parameters
pfix = cell(npar,1);
for i = 1:npar
    if isfield(cfg,pnam{i})
        pfix{i} = reshape(cfg.(pnam{i}),[psiz(i),1]);
    end
end

% define free parameters to fit
ifit = cell(npar,1);
pfit_ini = [];
pfit_min = [];
pfit_max = [];
n = 1;
for i = 1:npar
    if isempty(pfix{i}) % fully fixed parameter
        ifit{i} = n+[1:psiz(i)]-1;
        if numel(pini{i}) == 1
            pfit_ini = cat(1,pfit_ini,pini{i}*ones(psiz(i),1));
        elseif numel(pini{i}) == psiz(i)
            pfit_ini = cat(1,pfit_ini,pini{i}(:));
        else
            error('wrong size for initialization value!');
        end  
        pfit_min = cat(1,pfit_min,pmin(i)*ones(psiz(i),1));
        pfit_max = cat(1,pfit_max,pmax(i)*ones(psiz(i),1));
        n = n+psiz(i);
    elseif any(isnan(pfix{i})) % partly fixed parameter
        if numel(pfix{i}) ~= psiz(i)
            error('wrong size for fixed parameter!');
        end
        nfit = nnz(isnan(pfix{i}));
        ifit{i} = n+[1:nfit]-1;
        if numel(pini{i}) == 1
            pfit_ini = cat(1,pfit_ini,pini{i}*ones(nfit,1));
        elseif numel(pini{i}) == psiz(i) && all(isnan(pini{i}) ~= isnan(pfix{i}))
            pfit_ini = cat(1,pfit_ini,pini{i}(~isnan(pini{i})));
        else
            error('wrong size for initialization value!');
        end  
        pfit_min = cat(1,pfit_min,pmin(i)*ones(nfit,1));
        pfit_max = cat(1,pfit_max,pmax(i)*ones(nfit,1));
        n = n+nfit;
    end
end
nfit = length(pfit_ini);

% fit free parameters
if nfit > 0
    pval = fmincon(@fmin,pfit_ini,[],[],[],[],pfit_min,pfit_max,[], ...
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
out.aic  = -2*out.llh+2*nfit+2*nfit*(nfit+1)/(nobs-nfit+1); % model AIC
out.bic  = -2*out.llh+nfit*log(nobs); % model BIC

% get best-fitting posterior probabilities
out.ppost = get_ppost(phat{:});

% get confirmation probability conditioned on observed response
out.pconf = nan(nseq,2);
for iseq = 1:nseq
    out.pconf(iseq,:) = out.ppost(iseq,(resp(iseq)-1)*2+[1,2]);
end
out.pconf = bsxfun(@rdivide,out.pconf,sum(out.pconf,2));

% get accumulated evidence
[~,out.xevi,out.sevi,out.tevi] = get_ppost(phat{:});

% get full width at half maximum for key parameters
if get_fwhm
    % get full width at half maximum for sd_inf
    if numel(out.sd_inf) == 2
        [~,~,ls] = get_llh(phat{:});
        sd_inf_fwhm = nan(2,2);
        for t = 1:2
            sd_inf_fwhm(t,1) = fzero(@(p)ls(t)-get_llh_stype(t,out.ec_inf,out.el_inf,p,out.criter, ...
                out.eg_cnf,out.eb_cnf,out.sd_cnf)+log(0.5),[0.1,out.sd_inf(t)]);
            sd_inf_fwhm(t,2) = fzero(@(p)ls(t)-get_llh_stype(t,out.ec_inf,out.el_inf,p,out.criter, ...
                out.eg_cnf,out.eb_cnf,out.sd_cnf)+log(0.5),[out.sd_inf(t),2.0]);
        end
        out.sd_inf_fwhm = sd_inf_fwhm;
    end
    % get full width at half maximum for eb_cnf
    if numel(out.eb_cnf) == 2
        [~,~,ls] = get_llh(phat{:});
        eb_cnf_fwhm = nan(2,2);
        for t = 1:2
            eb_cnf_fwhm(t,1) = fzero(@(p)ls(t)-get_llh_stype(t,out.ec_inf,out.el_inf,out.sd_inf,out.criter, ...
                out.eg_cnf,p,out.sd_cnf)+log(0.5),[-2,out.eb_cnf(t)]);
            eb_cnf_fwhm(t,2) = fzero(@(p)ls(t)-get_llh_stype(t,out.ec_inf,out.el_inf,out.sd_inf,out.criter, ...
                out.eg_cnf,p,out.sd_cnf)+log(0.5),[out.eb_cnf(t),+2]);
        end
        out.eb_cnf_fwhm = eb_cnf_fwhm;
    end
    % get full width at half maximum for sd_cnf
    if numel(out.sd_cnf) == 2
        [~,~,ls] = get_llh(phat{:});
        sd_cnf_fwhm = nan(2,2);
        for t = 1:2
            try
                sd_cnf_fwhm(t,1) = fzero(@(p)ls(t)-get_llh_stype(t,out.ec_inf,out.el_inf,out.sd_inf,out.criter, ...
                    out.eg_cnf,out.eb_cnf,p)+log(0.5),[0,out.sd_cnf(t)]);
            catch
                sd_cnf_fwhm(t,1) = 0;
            end
            try
                sd_cnf_fwhm(t,2) = fzero(@(p)ls(t)-get_llh_stype(t,out.ec_inf,out.el_inf,out.sd_inf,out.criter, ...
                    out.eg_cnf,out.eb_cnf,p)+log(0.5),[out.sd_cnf(t),10]);
            catch
                sd_cnf_fwhm(t,2) = 10;
            end
        end
        out.sd_cnf_fwhm = sd_cnf_fwhm;
    end 
end

    function [f,pfit] = fmin(p)
    % OBJECTIVE FUNCTION = NEGATIVE MODEL LOG-LIKELIHOOD
    % create list of model parameters
    pfit = cell(npar,1);
    for i = 1:npar
        if isempty(pfix{i}) % free parameter
            pfit{i} = p(ifit{i});
        else % fixed parameter
            if any(isnan(pfix{i})) % partly fixed parameter
                pfit{i} = pfix{i};
                pfit{i}(isnan(pfit{i})) = p(ifit{i});
            else % fully fixed parameter
                pfit{i} = pfix{i};
            end
        end
    end
    % return negative model log-likelihood
    f = -get_llh(pfit{:});
    end

    function [llh,ppost,llh_stype] = get_llh(varargin)
    % GET MODEL LOG-LIKELIHOOD
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

    function [llh] = get_llh_stype(s,varargin)
    % GET MODEL LOG-LIKELIHOOD FOR EITHER SESSION
    [~,~,llh_stype] = get_llh(varargin{:});
    llh = llh_stype(s);
    end

    function [ppost,xevi,sevi,tevi] = get_ppost(ec_inf,el_inf,sd_inf,criter,eg_cnf,eb_cnf,sd_cnf)
    % GET MODEL POSTERIOR PROBABILITIES FOR EACH RESPONSE COMBINATION
    % column#1 -> resp=1 + confirm
    % column#2 -> resp=1 + opt-out
    % column#3 -> resp=2 + confirm
    % column#4 -> resp=2 + opt-out
    ppost_all = nan(nseq,4); % posterior probabilities
    xevi_all  = nan(nseq,1); % sequence evidence
    sevi_all  = nan(nseq,1); % internal noise
    tevi_all  = nan(nseq,1); % accumulated evidence conditioned on responses
    for s = 1:2
        % filter w.r.t. session type
        ifilt = find(stype == s);
        nfilt = length(ifilt);
        % get parameters
        ec_inf_sub = sub(ec_inf,s);
        el_inf_sub = sub(el_inf,s);
        sd_inf_sub = sub(sd_inf,s);
        criter_sub = sub(criter,s);
        eg_cnf_sub = sub(eg_cnf,s);
        eb_cnf_sub = sub(eb_cnf,s);
        sd_cnf_sub = sub(sd_cnf,s);
        % compute evidence compression
        ndev = 10; % development of infinite series
        xang = col(seqang(ifilt,:));
        xang = xang+0.5/besseli(0,ec_inf_sub)* ...
            sum(bsxfun(@times,sin(bsxfun(@times,1:ndev,xang*4)),besseli(1:ndev,ec_inf_sub)./(1:ndev)),2);
        xang = reshape(xang,[nfilt,maxlen]);
        xevi = sin(2*xang);
        xevi(isnan(xevi)) = 0;
        % compute leaky evidence accumulation
        xevi = nansum(bsxfun(@times,xevi,exp(-el_inf_sub).^(0:maxlen-1)),2);
        sevi = sqrt(arrayfun(@(n)sum(exp(-el_inf_sub).^(2*(0:n-1))),seqlen_lst(:))*sd_inf_sub^2);
        xevi_all(ifilt) = xevi;
        sevi_all(ifilt) = sevi(seqlen_ind(ifilt));
        % compute response probabilities
        presp = normcdf(xevi,criter_sub,sevi(seqlen_ind(ifilt)));
        presp = [presp,1-presp];
        % compute evidence and lottery gains
        gevi = (eg_cnf_sub <= 0)*exp(+eg_cnf_sub)+(eg_cnf_sub > 0); % for evidence
        glot = (eg_cnf_sub <= 0)+(eg_cnf_sub > 0)*exp(-eg_cnf_sub); % for lottery
        % truncate accumulated evidence wrt response
        tevi = nan(nfilt,1); % truncated evidence
        for u = 1:nfilt
            if resp(ifilt(u)) == 1
                tevi(u) = +tnmean(0,+inf,xevi(u)-criter_sub,sevi(seqlen_ind(ifilt(u))));
            else
                tevi(u) = -tnmean(-inf,0,xevi(u)-criter_sub,sevi(seqlen_ind(ifilt(u))));
            end
        end
        % define normalization factor
        switch norm_typ
            case 'p'
                fevi = fzero(@(f)mean(1./(1+exp(-f*tevi)))-p0(s),[1e-6,1e+6]);
            case 'l'
                fevi = logit(p0(s))/mean(tevi);
            otherwise
                error('Undefined normalization type!');
        end
        tevi = tevi*fevi; % normalize truncated evidence
        mevi = mean(tevi); % mean truncated evidence
        tevi_all(ifilt) = tevi;
        % compute posterior probabilities
        ppost = [];
        for r = 1:2
            pconf = nan(nfilt,1);
            for i = 1:length(seqlen_lst)
                j = find(seqlen(ifilt) == seqlen_lst(i));
                mx = xevi(j)*(3-2*r)*fevi*gevi;
                bx = criter_sub*(3-2*r)*fevi*gevi;
                sx = sevi(i)*fevi*gevi;
                by = logit(lotval(ifilt(j)))*glot+mevi*(gevi-glot)-eb_cnf_sub;
                sy = sd_cnf_sub;
                if sy == 0
                    pconf(j) = (normcdf((bx+by-mx)/sx)-normcdf((bx-mx)/sx))./normcdf((mx-bx)/sx);
                else
                    pconf(j) = mvncdf([mx-bx,bx+by-mx],[0,0],[sx^2,-sx^2;-sx^2,sx^2+sy^2])./normcdf((mx-bx)/sx);
                end
            end
            pconf = max(pconf,0);
            pconf = [1-pconf,pconf];
            ppost = cat(2,ppost,bsxfun(@times,pconf,presp(:,r)));
        end
        % append session
        ppost_all(ifilt,:) = ppost;
    end
    ppost = ppost_all;
    xevi  = xevi_all;
    sevi  = sevi_all;
    tevi  = tevi_all;
    end

end