function [expe] = KETACONF_gen_expe(nblck)

% check input arguments
if nargin < 1
    nblck = 6;
end
if mod(nblck,2) > 0
    error('number of blocks must be even!');
end

% add toolbox to path if necessary
if exist('HasConsecutiveValues') ~= 2
    addpath('./Toolboxes/Rand');
end

% number of repetitions of unique trials
nrep = nblck/2;

% define experimental factors
condtn_lst = [1,2]; % list of conditions (n=2)
seqcat_lst = [1,2]; % list of sequence categories (n=2)
seqsig_lst = [1,2,3]; % list of sequence signal strengths (n=3)
seqlen_lst = [4,8,12]; % list of sequence lengths (n=3)
lotval_lst = [0.6,0.7,0.8,0.9]; % list of lottery probability values (n=4)

% average of lottery probability values
lotval_avg = mean(lotval_lst);

% create block structure
blck = [];
for irep = 1:nrep
    % define set of unique trials
    condtn = kron(condtn_lst,ones(1,72));
    seqcat = repmat(kron(seqcat_lst,ones(1,36)),[1,2]);
    seqsig = repmat(kron(seqsig_lst,ones(1,12)),[1,4]);
    seqlen = repmat(kron(seqlen_lst,ones(1,4)),[1,12]);
    lotval = repmat(lotval_lst,[1,36]);
    % shuffle trials across titration steps
    while true
        iseq_inp = 1:numel(condtn);
        iseq_out = [];
        for i = 1:4
            found = false;
            for itry = 1:1+99*(i < 6)
                jseq_inp = iseq_inp;
                iseq_sup = [];
                for j = 1:2
                    while true
                        iseq_sub = [];
                        for icat = 1:2
                            for isig = 1:3
                                for ilen = 1:3
                                    iseq = find( ...
                                        seqcat == seqcat_lst(icat) & ...
                                        seqsig == seqsig_lst(isig) & ...
                                        seqlen == seqlen_lst(ilen));
                                    iseq = intersect(iseq,jseq_inp);
                                    iseq_sub = cat(2,iseq_sub,iseq(ceil(numel(iseq)*rand)));
                                end
                            end
                        end
                        if all(hist(condtn(iseq_sub),condtn_lst) == 9)
                            break
                        end
                    end
                    iseq_sup = cat(2,iseq_sup,iseq_sub);
                    jseq_inp = setdiff(jseq_inp,iseq_sub);
                end
                if ...
                        range(hist(lotval(iseq_sup),lotval_lst)) <= 6 && ...
                        abs(mean(lotval(iseq_sup))-mean(lotval_lst)) <= 0.05 && ...
                        range(grpstats(lotval(iseq_sup),condtn(iseq_sup),@mean)) <= 0.1
                    found = true;
                    break
                end
            end
            if ~found
                break
            end
            iseq_inp = jseq_inp;
            iseq_out = cat(2,iseq_out,iseq_sup);
        end
        if found
            break
        end
    end
    condtn = condtn(iseq_out);
    seqcat = seqcat(iseq_out);
    seqsig = seqsig(iseq_out);
    seqlen = seqlen(iseq_out);
    lotval = lotval(iseq_out);
    % shuffle trials within each titration step
    for i = 1:4
        iseq = (i-1)*36+(1:36);
        while true
            jseq = iseq(randperm(36));
            if ...
                    ~HasConsecutiveValues(condtn(jseq),4) && ...
                    ~HasConsecutiveValues(seqcat(jseq),4) && ...
                    ~HasConsecutiveValues(seqsig(jseq),4) && ...
                    ~HasConsecutiveValues(seqlen(jseq),4) && ...
                    ~HasConsecutiveValues(lotval(jseq),3) && ...
                    ~HasConsecutiveValues(lotval(jseq) < lotval_avg,4)
                break
            end
        end
        condtn(iseq) = condtn(jseq);
        seqcat(iseq) = seqcat(jseq);
        seqsig(iseq) = seqsig(jseq);
        seqlen(iseq) = seqlen(jseq);
        lotval(iseq) = lotval(jseq);
    end
    % assign shuffled trials to blocks
    n = length(blck);
    for i = 1:2
        iseq = (i-1)*72+(1:72);
        blck(n+i).condtn = condtn(iseq);
        blck(n+i).seqcat = seqcat(iseq);
        blck(n+i).seqsig = seqsig(iseq);
        blck(n+i).seqlen = seqlen(iseq);
        blck(n+i).lotval = lotval(iseq);
    end
end

% use final block as practice block
n = length(blck);
blck = blck([n,1:n]);
% no lottery comparison in 1st half of practice block
blck(1).condtn(1:36) = 0;

% create experiment structure
expe = [];
% store configuration information
expe.cfg = [];
expe.cfg.condtn_lst = condtn_lst;
expe.cfg.seqcat_lst = seqcat_lst;
expe.cfg.seqsig_lst = seqsig_lst;
expe.cfg.seqlen_lst = seqlen_lst;
expe.cfg.lotval_lst = lotval_lst;
% store block-wise information
expe.blck = blck;

end