function [out] = fit_model_noisesource(cfg)
%  FIT_MODEL_NOISESOURCE  Fit inference model with different noise sources
%
%  Usage: [out] = FIT_MODEL_NOISESOURCE(cfg)
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
%    * sd_ang = sensory noise s.d.
%    * sd_inf = inference noise s.d.
%    * sd_sel = selection noise s.d.
%    * criter = decision criterion
%    * plapse = response lapse rate
%    * kappa  = generative sequence coherence
%
%  Any combination of parameters can be provided as fields of the configuration
%  as fixed parameters that will not be fitted.
%
%  The function uses the interior-point algorithm of the fmincon function to
%  maximize the log-likelihood of the model.
%
%  It is wise to fix kappa to an arbitrary value, constant across all fits, so
%  that other parameter values (e.g., sd_inf) can be compared across fits.
%  Indeed, kappa is highly colinear with sd_inf (and sd_sel) and therefore
%  values of sd_inf and sd_sel can only be compared across fits for a same kappa
%  value.
%
%  Set cfg.only_incr = true for allowing only increases in noise sources under
%  ketamine when noise source levels are provided as fixed values for the
%  placebo condition.
%
%  Valentin Wyart <valentin.wyart@inserm.fr>

% check configuration parameters
if ~all(isfield(cfg,{'stype','seqang','catdir','resp'}))
    error('Incomplete configuration structure!');
end

% set default parameters
if ~isfield(cfg,'parmod')
    cfg.parmod = {};
end
if ~isfield(cfg,'only_incr')
    cfg.only_incr = false;
end

% get configuration parameters
stype  = cfg.stype(:); % session type => 1:placebo or 2:ketamine
seqang = cfg.seqang(:); % sequence angles (rad)
catdir = cfg.catdir(:); % categorization axis (rad)
resp   = cfg.resp(:); % response => 1 or 2
parmod = cfg.parmod; % list of parameters modulated by session type

% set useful function handles
col = @(x)x(:);
row = @(x)x(:)';
sub = @(x,k)x(min(k,numel(x)));

% get sequence lengths
seqlen = cellfun(@length,seqang);
[seqlen_lst,~,seqlen_ind] = unique(seqlen);

nseq = length(seqlen); % number of sequences
maxlen = max(seqlen); % maximum sequence length

% create matrix of angles
seqang = cellfun(@(x,n)cat(2,row(x(1:n)),nan(1,maxlen-n)), ...
    seqang,num2cell(seqlen),'UniformOutput',false);
seqang = cat(1,seqang{:});
% express angles as tilts from 1st category mean
seqtlt = mod(bsxfun(@minus,seqang,catdir+pi/4)+pi/2,pi)-pi/2;

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
% sensory noise
npar = npar+1;
pnam{npar,1} = 'sd_ang';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini(npar,1) = pi/180*10;
pmin(npar,1) = 0;
pmax(npar,1) = pi/180*90;
if ismember('sd_ang',parmod) && isfield(cfg,'sd_ang') && cfg.only_incr
    % restrict to ketamine increase when placebo is provided
    pini(npar) = cfg.sd_ang(1);
    pmin(npar) = cfg.sd_ang(1);
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
% selection noise
npar = npar+1;
pnam{npar,1} = 'sd_sel';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini(npar,1) = 10;
pmin(npar,1) = 0;
pmax(npar,1) = 100;
if ismember('sd_sel',parmod) && isfield(cfg,'sd_sel') && cfg.only_incr
    % restrict to ketamine increase when placebo is provided
    pini(npar) = cfg.sd_sel(1);
    pmin(npar) = cfg.sd_sel(1);
end
% decision criterion
npar = npar+1;
pnam{npar,1} = 'criter';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini(npar,1) = 0;
pmin(npar,1) = -5;
pmax(npar,1) = +5;
% lapse probability
npar = npar+1;
pnam{npar,1} = 'plapse';
psiz(npar,1) = 1+ismember(pnam{end},parmod);
pini(npar,1) = 0.01;
pmin(npar,1) = 0;
pmax(npar,1) = 0.99;
% generative coherence
npar = npar+1;
pnam{npar,1} = 'kappa';
psiz(npar,1) = 1;
pini(npar,1) = 0.5;
pmin(npar,1) = 0.01;
pmax(npar,1) = 10;

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
        % constrain to sd[placebo] < sd[ketamine]
        tnam = {'sd_ang','sd_inf','sd_sel'};
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

% get model log-likelihood
out.nfit = nfit; % number of fitted parameters
out.nobs = nobs; % number of observations
out.iobs = iobs; % indices of observations
out.llh  = get_llh(phat{:}); % model log-likelihood

% compute additional quality-of-fit metrics
[~,~,llh_s] = get_llh(phat{:});
nfit_s = nfit-numel(parmod);
nobs_s = [nnz(stype == 1),nnz(stype == 2)];
out.aic_stype = -2*llh_s+2*nfit_s+2*nfit_s*(nfit_s+1)./(nobs_s-nfit_s+1);
out.bic_stype = -2*llh_s+nfit_s*log(nobs_s);
out.aic = sum(out.aic_stype);
out.bic = sum(out.bic_stype);

% get best-fitting posterior probabilities
out.ppost = get_ppost(phat{:});

% compute statistics for decision noise *variance* (optional)
if isfield(cfg,'s2_dec_vec')
    svec = cfg.s2_dec_vec;
    lvec = nan(size(svec));
    for k = 1:numel(svec)
        lvec(k) = get_llh(0,0,sqrt(svec(k)),0,0,0.5);
    end
    pvec = exp(lvec-max(lvec));
    pvec = pvec/trapz(svec,pvec);
    out.s2_dec_vec = svec;
    out.s2_dec_llh = lvec;
    out.s2_dec_pdf = pvec;
    out.s2_dec_avg = trapz(svec,pvec.*svec);
    out.s2_dec_std = sqrt(trapz(svec,pvec.*(svec-out.s2_dec_avg).^2));
end

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
        % return negative model log-likelihood
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

    function [ppost,xpost,spost] = get_ppost(sd_ang,sd_inf,sd_sel,criter,plapse,kappa)
        xpost = nan(nseq,1); % decision signal
        spost = nan(nseq,1); % decision noise standard deviation
        ppost = nan(nseq,1); % posterior probability
        for s = 1:2
            % filter trials
            ifilt = find(stype == s);
            if numel(ifilt) > 0
                % get parameters
                sd_ang_sub = sub(sd_ang,s);
                sd_inf_sub = sub(sd_inf,s);
                sd_sel_sub = sub(sd_sel,s);
                criter_sub = sub(criter,s);
                plapse_sub = sub(plapse,s);
                % compute effect of sensory noise
                mu_llr = 2*kappa*exp(-2*sd_ang_sub^2)*cos(2*(seqtlt(ifilt,:)));
                s2_llr = 2*kappa^2*(1+exp(-8*sd_ang_sub^2)*cos(4*(seqtlt(ifilt,:))));
                sd_llr = sqrt(max(s2_llr-mu_llr.^2,0));
                % set missing values to zero
                mu_llr(isnan(mu_llr)) = 0;
                sd_llr(isnan(sd_llr)) = 0;
                % compute decision signal
                xpost(ifilt) = sum(mu_llr,2);
                % compute decision noise s.d.
                spost(ifilt) = sqrt( ...
                    sum(sd_llr.^2,2)+ ... % sensory noise
                    (seqlen(ifilt)-1)*sd_inf_sub^2+ ... % inference noise
                    sd_sel_sub^2); % selection noise
                % compute posterior probability
                ppost(ifilt) = normcdf(xpost(ifilt),criter_sub,spost(ifilt));
                ppost(ifilt) = ppost(ifilt)*(1-plapse_sub)+0.5*plapse_sub;
            end
        end
        ppost = [ppost,1-ppost];
    end

end