function [out] = read_dat_resp(subj,sess,cfg);
%  READ_DAT_RESP  Read response-locked EEG data of KETACONF study
%
%  Usage: [out] = READ_DAT_RESP(subj,sess,cfg)
%
%  where subj is the subject number
%        sess is the session number
%        cfg is a configuration structure
%
%  The configuration structure should include the following fields:
%    * toff = the time window (offset from the response onset)
%    * istp = decimation factor (1=no-decimation)
%
%  Valentin Wyart <valentin.wyart@inserm.fr>

% check input arguments
if nargin < 3
    error('missing input arguments!');
end
if ~all(isfield(cfg,{'toff','istp'}))
    error('incomplete configuration structure!');
end
if ~isfield(cfg,'trig')
    cfg.trig = 'resp';
end

rootfolder = '~/Dropbox/MATLAB/KETACONF'; % root folder

% define event-related window
ioff = fix(cfg.toff(1)*512)+(0:cfg.istp:fix(diff(cfg.toff)*512)-1);
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
k = 1;
for iblck = 1:nblck
    % load block
    load(fullfile(foldname,filename{iblck}));
    nseq = length(trg_c.blck.seqang);
    for iseq = 1:nseq
        dat_tmp = dat_c.trial{iseq};
        % initialize data array
        if k == 1
            nchan = size(dat_tmp,1);
            dat = nan(1000,nchan,ntime);
            seqang = cell(1000,1);
            catdir = nan(1000,1);
            resp   = nan(1000,1);
        end
        % crop data around event of interest
        ipos = trg_c.ipos.(cfg.trig)(iseq)+ioff;
        dat(k,:,:) = dat_tmp(:,ipos);
        seqang{k} = trg_c.blck.seqang{iseq};
        catdir(k) = trg_c.blck.catdir(iseq);
        resp(k)   = trg_c.blck.resp(iseq);
        k = k+1;
    end
end
dat(k:end,:,:) = [];
seqang(k:end) = [];
catdir(k:end) = [];
resp(k:end) = [];

% create output structure
out        = [];
out.toff   = toff;
out.chans  = dat_c.label;
out.nchan  = nchan;
out.dat    = dat;
out.seqang = seqang;
out.catdir = catdir;
out.mu     = mod(bsxfun(@plus,catdir,pi/4*[+1,-1]),pi);
out.kappa  = 0.5*ones(size(catdir));
out.resp   = resp;

end