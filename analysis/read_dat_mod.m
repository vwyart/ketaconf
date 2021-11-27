function [out] = read_dat_mod(subj,sess,cfg)
%  READ_DAT_MOD  Read stimulus-locked EEG data of KETACONF study
%
%  Usage: [out] = READ_DAT_MOD(subj,sess,cfg)
%
%  where subj is the subject number
%        sess is the session number
%        cfg is a configuration structure
%
%  The configuration structure should include the following fields:
%    * toff = the time window (offset from the response onset)
%    * istp = decimation factor (1=no-decimation)
%    * flowpass = low-pass frequency cutoff (Hz)
%    * fhighpass = high-pass frequency cutoff (Hz, optional)
%    * use_analytic_signal = use analytic signal instead of real signal?
%
%  When use_analysis_signal == true, the function performs the Hilbert transform
%  of the band-pass EEG signal to extract the real and imaginary signals that
%  will be later used for multivariate pattern analyses.
%
%  Valentin Wyart <valentin.wyart@inserm.fr>

% check input arguments
if nargin < 3
    error('missing input arguments!');
end
if ~all(isfield(cfg,{'toff','istp','flowpass','use_analytic_signal'}))
    error('incomplete configuration structure!');
end
if ~isfield(cfg,'fhighpass')
    cfg.fhighpass = [];
end

rootfolder = '~/Dropbox/MATLAB/KETACONF'; % root folder

% define event-related window
ioff = fix(cfg.toff(1)*512):cfg.istp:fix(cfg.toff(2)*512);
toff = ioff/512;
ntime = length(toff);

% locate files
foldname = sprintf('%s/Data/EEG/S%03d',rootfolder,subj);
filename = {};
d = dir(foldname);
for ifile = 1:length(d)
    if ~isempty(strfind(d(ifile).name,'.mat'))
        s = sscanf(d(ifile).name,'KETACONF_S%03d_sess%d_b%d.mat');
        if length(s) == 3 && s(2) == sess
            filename{end+1} = d(ifile).name;
        end
    end
end
nblck = numel(filename);

% read data
seqcat = nan(500,1);
seqang = cell(500,1);
catdir = nan(500,1);
resp   = nan(500,1);
conf   = nan(500,1);
tresp  = nan(500,1);
tconf  = nan(500,1);
kseq = 1;
kcue = 1;
for iblck = 1:nblck
    % load block
    load(fullfile(foldname,filename{iblck}));
    nseq = length(trg_c.blck.seqang);
    for iseq = 1:nseq
        ncue = length(trg_c.blck.seqang{iseq});
        if length(trg_c.ipos.card{iseq}) ~= ncue
            error('mismatching event information!');
        end
        dat_tmp = dat_c.trial{iseq};
        if kcue == 1
            nchan = size(dat_tmp,1);
        end
        if ~isempty(cfg.flowpass) && isempty(cfg.fhighpass)
            % low-pass filter using 6th-order Butterworth IIR filter
            dat_tmp = ft_preproc_lowpassfilter(dat_tmp,dat_c.fsample,cfg.flowpass,6,'but');
        elseif ~isempty(cfg.flowpass) && ~isempty(cfg.fhighpass)
            dat_tmp = ft_preproc_bandpassfilter(dat_tmp,dat_c.fsample, ...
                [cfg.flowpass,cfg.fhighpass],4,'but');
        elseif isempty(cfg.lowpass) && ~isempty(cfg.highpass)
            dat_tmp = ft_preproc_highpassfilter(dat_tmp,dat_c.fsample,cfg.fhighpass,6,'but');
        end
        
        % compute analytic signal using Hilbert transform
        if cfg.use_analytic_signal
            dat_tmp = cat(1, ...
                ft_preproc_hilbert(dat_tmp,'real'), ...
                ft_preproc_hilbert(dat_tmp,'imag'));
        end
        % initialize data array
        if kcue == 1
            nfeat = size(dat_tmp,1);
            dat = nan(5000,nfeat,ntime);
            dat_seq = nan(5000,1);
            dat_pos = nan(5000,1);
            dat_len = nan(5000,1);
            dir_seq = nan(5000,1);
            ang_prr = nan(5000,1);
            ang_prv = nan(5000,1);
            ang_cur = nan(5000,1);
            ang_nxt = nan(5000,1);
        end
        seqcat(kseq) = trg_c.blck.seqcat(iseq);
        seqang{kseq} = trg_c.blck.seqang{iseq};
        catdir(kseq) = trg_c.blck.catdir(iseq);
        resp(kseq) = trg_c.blck.resp(iseq);
        conf(kseq) = trg_c.blck.conf(iseq);
        tresp(kseq) = trg_c.time.resp(iseq)-trg_c.time.probe(iseq);
        tconf(kseq) = trg_c.time.conf(iseq)-trg_c.time.offer(iseq);
        for icue = 1:ncue
            % crop data around event of interest
            ipos = trg_c.ipos.card{iseq}(icue)+ioff;
            dat(kcue,:,:) = dat_tmp(:,ipos);
            dat_seq(kcue) = kseq;
            dat_pos(kcue) = icue;
            dat_len(kcue) = ncue;
            dir_seq(kcue) = trg_c.blck.catdir(iseq);
            if icue > 2
                ang_prr(kcue) = trg_c.blck.seqang{iseq}(icue-2);
            end
            if icue > 1
                ang_prv(kcue) = trg_c.blck.seqang{iseq}(icue-1);
            end
            ang_cur(kcue) = trg_c.blck.seqang{iseq}(icue);
            if icue < ncue
                ang_nxt(kcue) = trg_c.blck.seqang{iseq}(icue+1);
            end
            kcue = kcue+1;
        end
        kseq = kseq+1;
    end
end
dat(kcue:end,:,:) = [];
dat_seq(kcue:end) = [];
dat_pos(kcue:end) = [];
dat_len(kcue:end) = [];
dir_seq(kcue:end) = [];
ang_prr(kcue:end) = [];
ang_prv(kcue:end) = [];
ang_cur(kcue:end) = [];
ang_nxt(kcue:end) = [];
seqcat(kseq:end)  = [];
seqang(kseq:end)  = [];
catdir(kseq:end)  = [];
resp(kseq:end)    = [];
conf(kseq:end)    = [];
tresp(kseq:end)   = [];
tconf(kseq:end)   = [];

% create output structure
out         = [];
out.toff    = toff;
out.chans   = dat_c.label;
out.nchan   = nchan;
out.nfeat   = nfeat;
out.dat     = dat;
out.dat_seq = dat_seq;
out.dat_pos = dat_pos;
out.dat_len = dat_len;
out.dir_seq = dir_seq;
out.ang_prr = ang_prr;
out.ang_prv = ang_prv;
out.ang_cur = ang_cur;
out.ang_nxt = ang_nxt;
out.seqcat  = seqcat;
out.seqang  = seqang;
out.catdir  = catdir;
out.mu      = mod(bsxfun(@plus,catdir,pi/4*[+1,-1]),pi);
out.resp    = resp;
out.conf    = conf;
out.tresp   = tresp;
out.tconf   = tconf;

end