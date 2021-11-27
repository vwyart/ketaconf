function [out] = gen_seqang(cfg)
%  GEN_SEQANG  Generate sequence angles with titrated evidence
%
%  Usage: [out] = gen_seqang(cfg)
%
%  where cfg - configuration structure
%        out - output structure
%
%  The configuration structure needs to contain the following fields:
%    * sd_cue     - current estimate of evidence s.d. per cue
%    * pcor       - target p(correct)
%    * seqsig_lst - list of signal strengths
%    * seqlen_lst - list of sequence lengths
%    * seqcat     - sequence categories (n=36)
%    * seqsig     - sequence signal strengths (n=36)
%    * seqlen     - sequence lengths (n=36)
%
%  The output structure contains the following fields:
%    * mu_cue - mean evidence per cue (nat)
%    * seqang - sequence angles (rad)
%    * catdir - categorization directions (rad)
%
%  Valentin Wyart <valentin.wyart@ens.fr>

% get configuration parameters
sd_cue     = cfg.sd_cue; % current estimate of evidence s.d. per cue
pcor       = cfg.pcor; % target p(correct)
seqsig_lst = cfg.seqsig_lst; % list of signal strengths
seqlen_lst = cfg.seqlen_lst; % list of sequence lengths
seqcat     = cfg.seqcat; % sequence categories (n=36)
seqsig     = cfg.seqsig; % sequence signal strengths (n=36)
seqlen     = cfg.seqlen; % sequence lengths (n=36)

% define discretized orientations
xv = linspace(-pi/2,+pi/2,1001);

% define function handles
col     = @(x)x(:); % columnize array
vmpdf   = @(x,t,k)exp(k.*cos(2*(x-t)))./(pi*besseli(0,k)); % von Mises pdf
vmrnd   = @(n,t,k)itrnd(n,xv,vmpdf(xv,t,k),'interp'); % von Mises random generator
modcc   = @(x)mod(x+pi/2,pi)-pi/2; % pi-periodic modulus
shuffle = @(x)x(randperm(numel(x))); % shuffle array

% define scaling of signal strengths w.r.t. target p(correct)
pcor_fun = @(s)mean(col(normcdf(s*bsxfun(@times,seqsig_lst',sqrt(seqlen_lst)))));
s0 = fzero(@(s)pcor_fun(s)-pcor,0);

% compute evidence mean and generative concentration per signal strength
mu_cue_lst = sd_cue*s0*seqsig_lst;
kappa_lst = nan(1,3);
xi = pi/2*(-1:0.001:+1);
for isig = 1:3
    mu_cue_fun = @(k)trapz(xi,2*k*sin(2*xi).*vmpdf(xi,+pi/4,k));
    kappa_lst(isig) = fzero(@(k)mu_cue_fun(k)-mu_cue_lst(isig),[0,10]);
end

% generate categorization directions
catdir = [];
for i = 1:3
    catdir_tmp = pi*((1:12)-0.5)/12;
    while true
        catdir_tmp = shuffle(catdir_tmp);
        if min(abs(modcc(diff(catdir_tmp)))) > pi/12
            break
        end
    end
    catdir = cat(2,catdir,catdir_tmp);
end

% generate sequences angles
seqang = cell(1,36);
for isig = 1:3
    for ilen = 1:3
        % locate sequences
        iseq = find(seqsig == seqsig_lst(isig) & seqlen == seqlen_lst(ilen));
        nseq = length(iseq);
        % generate 10^5 sequences
        x = vmrnd([1e5,seqlen_lst(ilen)],+pi/4,kappa_lst(isig));
        % compute sequence-wise statistics
        y = mean(2*kappa_lst(isig)*sin(2*x),2); % mean evidence per cue
        [d,k] = vmpar(x,2,5); % mean direction & concentration
        % compute error ranks
        [~,~,yerr] = unique(abs(y-mu_cue_lst(isig)));
        [~,~,derr] = unique(abs(modcc(d-pi/4)));
        [~,~,kerr] = unique(abs(k-kappa_lst(isig)));
        % identify sequences with smallest sums of error ranks
        [~,irnk] = sort(yerr+derr+kerr);
        x = x(shuffle(irnk(1:nseq)),:);
        x = bsxfun(@times,x,col(3-2*seqcat(iseq)));
        x = bsxfun(@plus,x,col(catdir(iseq)));
        seqang(iseq) = mat2cell(mod(x,pi),ones(1,nseq));
    end
end

% create output structure
out        = [];
out.mu_cue = mu_cue_lst(seqsig); % mean evidence per cue
out.seqang = seqang; % sequence angles (rad)
out.catdir = catdir; % categorization directions (rad)

end