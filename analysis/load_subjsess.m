%  load_subjsess  Load subject/session behavioral data for KETACONF study
%
%  This script loads the subject- and session-wise behavioral data using the
%  subj and sess variables, respectively. Both need to be found in the calling
%  workspace for this script to run.
%
%  Valentin Wyart <valentin.wyart@ens.fr>

% check for variables
if ~exist('subj','var') || ~exist('sess','var')
    error('Missing required variable(s) in workspace!');
end

% locate and load behavior
foldname = sprintf('../data/behavior/S%03d',subj+100);
d = dir(foldname);
i = find( ...
    cellfun(@(s)~isempty(strfind(s,sprintf('sess%d',sess))),{d.name}) & ...
    cellfun(@(s) isempty(strfind(s,'blck')),{d.name}) & ...
    cellfun(@(s)~isempty(strfind(s,'.mat')),{d.name}));
if numel(i) ~= 1
    error('Unable to find datafile for S%03d!',subj+100);
end
filename = d(i).name;
load(fullfile(foldname,filename));
% clear temporary variables
clear foldname filename

% locate blocks of interest
iblck = find(cellfun(@(x)nnz(x == 0),{rslt.resp}) == 0);
iblck = setdiff(iblck,1); % exclude titration block
nblck = length(iblck); % number of blocks

% get sequence-wise information
condtn = cat(2,blck(iblck).condtn); % delay condition => 1:short 2:long
seqcat = cat(2,blck(iblck).seqcat); % sequence category => 1:orange 2:blue
seqsig = cat(2,blck(iblck).seqsig); % sequence signal strength => 1:low 2:middle 3:high
seqlen = cat(2,blck(iblck).seqlen); % sequence length (in samples)
lotval = cat(2,blck(iblck).lotval); % lottery probability value
seqang = cat(2,blck(iblck).seqang); % sequence angles (in radians)
catdir = cat(2,blck(iblck).catdir); % categorization direction (in radians)
mu_cue = cat(2,blck(iblck).mu_cue); % titration level
n_titr = kron(1:2*nblck,ones(1,36)); % titration step

% get responses and obtained rewards
resp   = cat(2,rslt(iblck).resp);
conf   = cat(2,rslt(iblck).conf);
rwrd   = cat(2,rslt(iblck).rwrd);

% recompute response/confirmation latencies
tresp  = [];
tconf  = [];
ntoofast = 0;
for jblck = iblck
    % get response latencies
    t0 = zeros(1,72);
    te = zeros(1,72);
    i0 = evnt(jblck).iseq(strcmp(evnt(jblck).type,'PROBE'));
    ie = evnt(jblck).iseq(strcmp(evnt(jblck).type,'RESP'));
    t0(i0) = evnt(jblck).time(strcmp(evnt(jblck).type,'PROBE'));
    te(ie) = evnt(jblck).time(strcmp(evnt(jblck).type,'RESP'));
    tresp = cat(2,tresp,(te-t0).*(t0>0 & te>0));
    % get confirmation latencies
    t0 = zeros(1,72);
    te = zeros(1,72);
    i0 = evnt(jblck).iseq(strcmp(evnt(jblck).type,'OFFER'));
    ie = evnt(jblck).iseq(strcmp(evnt(jblck).type,'CONF'));
    t0(i0) = evnt(jblck).time(strcmp(evnt(jblck).type,'OFFER'));
    te(ie) = evnt(jblck).time(strcmp(evnt(jblck).type,'CONF'));
    tconf = cat(2,tconf,(te-t0).*(t0>0 & te>0));
    ntoofast = ntoofast+nnz(cellfun(@(s)strcmp(s,'RESP_EARLY'),evnt(jblck).type));
end
% clear temporary variables
clear jblck t0 te i0 ie

% exclude bad trials
isbad = conf <= 0;
condtn(isbad) = [];
seqcat(isbad) = [];
seqsig(isbad) = [];
seqlen(isbad) = [];
lotval(isbad) = [];
seqang(isbad) = [];
catdir(isbad) = [];
mu_cue(isbad) = [];
n_titr(isbad) = [];
resp(isbad)   = [];
tresp(isbad)  = [];
conf(isbad)   = [];
tconf(isbad)  = [];
rwrd(isbad)   = [];
% clear temporary variables
clear isbad

% is response correct?
iscor = resp == seqcat;

% get number of sequences
nseq = length(condtn);
