%% Analysis script for response preparation power analysis
%
%  This script affords to define a frequency band of interest for mu activity.
%  This script also affords to define a preferred pruning of discriminant features
%  for decoding responses (and, subsequently, latent decision signals).
%
%  The script computes two decoding metrics:
%    * a parametric metric based on the z-scored discriminant feature
%    * a non-parametric metric based on the area under the ROC curve
%    * a parametric metric for the encoding of the decision variable (DV)
%
%  Valentin Wyart <valentin.wyart@inserm.fr>

% clear MATLAB workspace
clear all
close all
clc

% add Fieldtrip toolbox to MATLAB path
% addpath('~/Dropbox/MATLAB/Toolboxes/fieldtrip-20180404');
% ft_defaults;

% get list of subject numbers
subjlist = 100+get_sesstype(1,true);
nsubj = length(subjlist);

% get list of session numbers
sessname = {'placebo','ketamine'};
sesslist = zeros(nsubj,2);
[~,sesslist(:,1)] = get_sesstype(1,true); % placebo session number
[~,sesslist(:,2)] = get_sesstype(2,true); % ketamine session number

% configure data reading
cfg      = [];
cfg.toff = [-3.0,+1.0]; % time window around response (in secs)
cfg.istp = 1; % time decimation step size (from 512 Hz)
cfg.trig = 'resp'; % trigger event of interest

% define frequencies of interest
% => log-spaced from 8 to 32 Hz (21 frequencies)
ftarg = 8*2.^(0:0.05:2);
flist = 1:1/(cfg.toff(2)-cfg.toff(1)):128;
itarg = [];
for ifreq = 1:length(ftarg)
    [e,i] = min(abs(flist-ftarg(ifreq)));
    itarg(ifreq) = i;
end
freqlist = flist(itarg)
nfreqs = length(freqlist);

% define pruning of discriminant features
% => log-spaced from 1 to 64 features (7 prunings)
prunlist = 8; %2.^(0:6);
npruns = numel(prunlist);

nfolds  = 10; % number of cross-validation folds (default=10)
roctype = 'onevsall'; % ROC type (unimportant)
flag_ci = 'sum'; % flag for analyzing consistent/inconsistent evidence

hbar = waitbar(0,'');
set(get(findobj(hbar,'Type','Axes'),'Title'),'FontSize',16);

% pre-initialize analysis output
xpow_all = cell(nsubj,2); % mean spectral power
zprj_all = cell(nsubj,2); % parametric (z-score) decoding metric (chance = 0)
proc_all = cell(nsubj,2); % non-parametric (ROC) decoding metric (chance = 0.5)
breg_all = cell(nsubj,2); % parametric DV encoding metric (chance = 0)
bcon_all = cell(nsubj,2); % parametric DV encoding metric for consistent evidence
binc_all = cell(nsubj,2); % parametric DV encoding metric for inconsistent evidence

bcon_all_1st = cell(nsubj,2);
bcon_all_lst = cell(nsubj,2);
binc_all_1st = cell(nsubj,2);
binc_all_lst = cell(nsubj,2);

col = @(x)x(:);

for isubj = 1:nsubj
    
    for isess = 1:2
        waitbar(0,hbar,sprintf('Processing S%03d (%s)',subjlist(isubj),sessname{isess}));
        
        % read data
        dat = read_dat_resp(subjlist(isubj),sesslist(isubj,isess),cfg);
        
        % get sequence-wise LLRs
        seqllr = cellfun(@minus,dat.seqang,num2cell(dat.catdir),'UniformOutput',false);
        seqllr = cellfun(@(x)sin(2*x),seqllr,'UniformOutput',false);
        smpllr = seqllr;
        seqllr = cellfun(@sum,seqllr);
        
        % get subject responses
        resp = dat.resp;
        
        % split sequence-wise LLRs into consistent/inconsistent wrt upcoming response
        nseq = numel(smpllr);
        switch flag_ci
            case 'sum'
                seqllr_con = zeros(nseq,1); seqllr_con_1st = zeros(nseq,1); seqllr_con_lst = zeros(nseq,1);
                seqllr_inc = zeros(nseq,1); seqllr_inc_1st = zeros(nseq,1); seqllr_inc_lst = zeros(nseq,1);
                for iseq = 1:nseq
                    x = smpllr{iseq}.*(3-2*resp(iseq));
                    seqllr_con(iseq) = sum(x(x > 0));
                    seqllr_inc(iseq) = sum(x(x < 0));
                    x1st = x(1:4);
                    seqllr_con_1st(iseq) = sum(x1st(x1st > 0));
                    seqllr_inc_1st(iseq) = sum(x1st(x1st < 0));
                    xlst = x(end:-1:end-3);
                    seqllr_con_lst(iseq) = sum(xlst(xlst > 0));
                    seqllr_inc_lst(iseq) = sum(xlst(xlst < 0));
                end
            case 'sep'
                seqllr_con = []; seqpos_con = []; seqneg_con = []; seqind_con = [];
                seqllr_inc = []; seqpos_inc = []; seqneg_inc = []; seqind_inc = [];
                for iseq = 1:nseq
                    x = smpllr{iseq}.*(3-2*resp(iseq));
                    ncue = numel(x);
                    icon = find(x > 0);
                    seqllr_con = cat(1,seqllr_con,col(x(icon)));
                    seqpos_con = cat(1,seqpos_con,icon(:));
                    seqneg_con = cat(1,seqneg_con,ncue-icon(:)+1);
                    seqind_con = cat(1,seqind_con,iseq(ones(numel(icon),1)));
                    iinc = find(x < 0);
                    seqllr_inc = cat(1,seqllr_inc,col(x(iinc)));
                    seqpos_inc = cat(1,seqpos_inc,iinc(:));
                    seqneg_inc = cat(1,seqneg_inc,ncue-iinc(:)+1);
                    seqind_inc = cat(1,seqind_inc,iseq(ones(numel(iinc),1)));
                end
            otherwise
                error('Undefined flag!');
        end
        
        % get lists of channels and time samples
        if isubj == 1 && isess == 1
            chanlist = dat.chans; % electrode names
            timelist = dat.toff([find(dat.toff >= -2,1),find(dat.toff >= 0,1)]);
        end
        
        % configure multi-tapering time-frequency decomposition
        % => 8 cycles per window, 25% frequency smoothing (3 tapers)
        varglist = {};
        varglist = cat(2,varglist,{'taper','dpss'});
        varglist = cat(2,varglist,{'timeoi',timelist});
        varglist = cat(2,varglist,{'freqoi',freqlist});
        varglist = cat(2,varglist,{'timwin',8./freqlist});
        varglist = cat(2,varglist,{'tapsmofrq',freqlist./4});
        varglist = cat(2,varglist,{'dimord','tap_chan_freq_time'});
        varglist = cat(2,varglist,{'verbose',0});
        
        % perform multi-tapering time-frequency decomposition
        ntrl = size(dat.dat,1);
        xpow = nan(ntrl,numel(chanlist),numel(freqlist),numel(timelist),'single');
        for itrl = 1:ntrl
            if mod(itrl,10) == 0
                waitbar(itrl/ntrl,hbar);
            end
            powspctrm = ft_specest_mtmconvol(squeeze(dat.dat(itrl,:,:)),dat.toff,varglist{:});
            powspctrm = nanmean(abs(powspctrm).^2,1); % average raw power across tapers
            xpow(itrl,:,:,:) = 10*log10(powspctrm); % store power in dB units
        end
        
        % store mean power
        xpow_all{isubj,isess} = squeeze(mean(xpow,1));
        
        % decode responses from spectral power using pruned discriminant features
        zprj = nan(npruns,nfreqs); % z-score decoding metric
        proc = nan(npruns,nfreqs); % ROC decoding metric
        breg = nan(npruns,nfreqs); % DV encoding metric
        bcon = nan(npruns,nfreqs); % DV encoding metric for consistent evidence
        binc = nan(npruns,nfreqs); % DV encoding metric for inconsistent evidence
        % ! only defined if flag_ci = 'sep'
        bcon_1st = nan(npruns,nfreqs);
        bcon_lst = nan(npruns,nfreqs);
        binc_1st = nan(npruns,nfreqs);
        binc_lst = nan(npruns,nfreqs);
        xlab = dat.resp(:);
        for ifreq = 1:nfreqs
            xdat = xpow(:,:,ifreq,2);
            for iprun = 1:npruns
                out = classify_proj(xlab,xdat,nfolds,roctype,prunlist(iprun));
                % compute z-score decoding metric
                xprj = out.xprj.*(3-2*out.lab);
                zprj(iprun,ifreq) = mean(xprj)/std(xprj);
                % get ROC decoding metric
                proc(iprun,ifreq) = out.proc;
                % compute DV encoding metric
                xreg = [seqllr.*(3-2*out.lab),ones(numel(out.lab),1)];
                yreg = xprj/std(xprj);
                b = xreg\yreg;
                breg(iprun,ifreq) = b(1);
                % split consistent/inconsistent wrt upcoming response
                switch flag_ci
                    case 'sum'
                        xreg = [seqllr_con,seqllr_inc,ones(size(seqllr))];
                        b = xreg\yreg;
                        bcon(iprun,ifreq) = b(1);
                        binc(iprun,ifreq) = b(2);
                        xreg = [seqllr_con_1st,seqllr_con_lst, ...
                            seqllr_inc_1st,seqllr_inc_lst,ones(size(seqllr))];
                        b = xreg\yreg;
                        bcon_1st(iprun,ifreq) = b(1);
                        bcon_lst(iprun,ifreq) = b(2);
                        binc_1st(iprun,ifreq) = b(3);
                        binc_lst(iprun,ifreq) = b(4);
                    case 'sep'
                        % analyze consistent evidence
                        xtmp = [seqllr_con,ones(size(seqllr_con))];
                        ytmp = yreg(seqind_con);
                        b = xtmp\ytmp;
                        bcon(iprun,ifreq) = b(1);
                        i = find(seqpos_con <= 4); b = xtmp(i,:)\ytmp(i);
                        bcon_1st(iprun,ifreq) = b(1);
                        i = find(seqneg_con <= 4); b = xtmp(i,:)\ytmp(i);
                        bcon_lst(iprun,ifreq) = b(1);
                        % analyze inconsistent evidence
                        xtmp = [seqllr_inc,ones(size(seqllr_inc))];
                        ytmp = yreg(seqind_inc);
                        b = xtmp\ytmp;
                        binc(iprun,ifreq) = b(1);
                        i = find(seqpos_inc <= 4); b = xtmp(i,:)\ytmp(i);
                        binc_1st(iprun,ifreq) = b(1);
                        i = find(seqneg_inc <= 4); b = xtmp(i,:)\ytmp(i);
                        binc_lst(iprun,ifreq) = b(1);
                end
            end
        end
        zprj_all{isubj,isess} = zprj;
        proc_all{isubj,isess} = proc;
        breg_all{isubj,isess} = breg;
        bcon_all{isubj,isess} = bcon;
        binc_all{isubj,isess} = binc;
        bcon_all_1st{isubj,isess} = bcon_1st;
        bcon_all_lst{isubj,isess} = bcon_lst;
        binc_all_1st{isubj,isess} = binc_1st;
        binc_all_lst{isubj,isess} = binc_lst;
        
        % save results to disk
        filename = sprintf('./data_analysis_mu_step1_%s.mat',flag_ci);
        save(filename, ...
            'subjlist','sessname','chanlist','freqlist','timelist','prunlist', ...
            'xpow_all','zprj_all','proc_all','breg_all','bcon_all','binc_all', ...
            'bcon_all_1st','bcon_all_lst','binc_all_1st','binc_all_lst');

    end
    
end
close(hbar);
