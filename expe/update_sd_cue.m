function [out] = update_sd_cue(cfg)

if ~all(isfield(cfg,{'mu_cue','seqcat','seqlen','resp'}))
    error('Incomplete configuration structure!');
end
if ~all(isfield(cfg,{'sd_val','sd_llf'}))
    
end

mu_cue = cfg.mu_cue;
seqcat = cfg.seqcat;
seqlen = cfg.seqlen;
resp   = cfg.resp;
sd_val = cfg.sd_val;
sd_llf = cfg.sd_llf;

nseq = length(resp);
iseq = find(ismember(resp,[1,2]));

iresp = sub2ind([nseq,2],iseq,resp(iseq));

out        = [];
out.sd_val = sd_val;
out.sd_llf = sd_llf;
for i = 1:length(sd_val)
    out.sd_llf(i) = out.sd_llf(i)+get_llh(sd_val(i));
end
out.sd_llf = out.sd_llf-max(out.sd_llf);
out.sd_cue = get_xmax(sd_val,out.sd_llf);

    function [llh] = get_llh(p)
        presp = get_presp(p);
        llh = sum(log(max(presp(iresp),eps)));
    end

    function [presp] = get_presp(sd_cue)
        mu_inf = mu_cue.*seqlen.*(3-2*seqcat);
        sd_inf = sd_cue.*sqrt(seqlen);
        presp = normcdf(mu_inf./sd_inf);
        presp = cat(2,presp,1-presp);
    end

end

function [xmax] = get_xmax(x,y)
fmin = @(xi)-interp1(x,y,xi,'spline');
xmax = fminbnd(fmin,min(x),max(x));
end