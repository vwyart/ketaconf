
%% Run simulations of model of inference and confidence with pre-commitments
%
%  Multiple simulations can be run using different parameter sets. Different
%  parameters need to be provided per simulation:
%    * cs_prb = pre-commitment probability
%    * cs_1st = 1st stimulus for possible pre-commitment (default = 5)
%    * sd_inf = inference noise s.d.
%    * el_inf = inference leak (in log-units, 0=leak-free)
%    * th_cnf = confidence threshold
%    * sd_obs = EEG observation noise (s.d.)
%
%  In the simulations below, the confidence threshold th_cnf (which controls
%  the overall decision uncertainty of the model) is increased for parameter
%  sets which include pre-commitments (cs_prb > 0).
%
%  sd_obs only controls the EEG observation noise for 'decoding' variables from
%  the model. It is set to 1.50 to reach coding precisions that match the ones
%  obtained when decoding the same variables from EEG data.
%
%  Valentin Wyart <valentin.wyart@inserm.fr>

% clear workspace
clear all java
close all hidden
clc

% define parameter sets to be tested
cs_prb_set = [0.00,0.05,0.10,0.15]; % pre-commitment probability
cs_1st_set = [5   ,5   ,5   ,5   ]; % 1st stimulus for possible pre-commitment
sd_inf_set = [0.40,0.40,0.40,0.40]; % inference noise s.d.
el_inf_set = [0.15,0.15,0.15,0.15]; % inference leak
th_cnf_set = [0.80,1.00,1.00,1.00]; % confidence threshold
sd_obs_set = [1.50,1.50,1.50,1.50]; % EEG observation noise s.d.

% define other required parameters
nsim = 1e2; % number of simulations per parameter set
nfit = 1; % number of simulations to be used for model fits (nfit <= nsim)
fit_sd_inf = false; % fit inference noise?

nset = numel(cs_prb_set);

% get list of subject numbers
subjlist = get_sesstype(1,true);
nsubj = length(subjlist);

% get list of session numbers
sesslist = zeros(nsubj,2);
[~,sesslist(:,1)] = get_sesstype(1,true); % placebo session number
[~,sesslist(:,2)] = get_sesstype(2,true); % ketamine session number

r2z = @(r)0.5*log((1+r)./(1-r)); % Fisher transform for correlation coefficient

% create output:
% * coding of consistent/conflicting evidence
zdec_cons      = nan(nsubj,6,2,nset); zerr_cons      = nan(nsubj,6,2,nset);
zdec_cons_less = nan(nsubj,6,2,nset); zerr_cons_less = nan(nsubj,6,2,nset);
zdec_cons_more = nan(nsubj,6,2,nset); zerr_cons_more = nan(nsubj,6,2,nset);
% * response preparation activity
xresp      = nan(nsubj,12,nset); eresp = nan(nsubj,12,nset);
xresp_less = nan(nsubj,12,nset);
xresp_more = nan(nsubj,12,nset);
% * best-fitting inference parameters
sd_inf_hat = nan(nsubj,nset);
el_inf_hat = nan(nsubj,nset);

xtot_com = []; % accumulated evidence for committed trials
xtot_reg = []; % accumulated evidence for regular trials

for isubj = 1:nsubj
    fprintf('Processing S%03d...\n',100+subjlist(isubj));
    
    % read data
    dat        = [];
    dat.stype  = []; % session type (1:placebo or 2:ketamine)
    dat.seqcat = []; % sequence category (1 or 2)
    dat.seqang = {}; % sequence angles
    dat.catdir = []; % categorization axis
    dat.resp   = []; % response (1 or 2)
    for stype = 1:2 % 1:placebo 2:ketamine
        % load data
        subj = subjlist(isubj); % subject number
        sess = sesslist(isubj,stype); % session number
        load_subjsess;
        % append session-wise data
        dat.stype  = cat(2,dat.stype,stype(ones(1,nseq)));
        dat.seqcat = cat(2,dat.seqcat,seqcat);
        dat.seqang = cat(2,dat.seqang,seqang);
        dat.catdir = cat(2,dat.catdir,catdir);
        dat.resp   = cat(2,dat.resp,resp);
    end
    
    if fit_sd_inf
        % fit inference noise
        cfg_fit = dat;
        cfg_fit.parmod = {'sd_inf'};
        cfg_fit.criter = 0;
        cfg_fit.ec_inf = 0;
        out_fit = fit_model_resp(cfg_fit);
        sd_inf_ref = out_fit.sd_inf(1);
    end
    
    % get useful fields
    stype  = dat.stype;
    seqang = dat.seqang;
    catdir = dat.catdir;
    seqcat = dat.seqcat;

    % duplicate sequences for multiple simulations
    stype  = repmat(stype,[1,nsim]);
    seqang = repmat(seqang,[1,nsim]);
    catdir = repmat(catdir,[1,nsim]);
    seqcat = repmat(seqcat,[1,nsim]);
    
    % total number of sequences
    nseq = numel(seqang);
    
    for iset = 1:nset
        fprintf('  * testing parameter set %d/%d...\n',iset,nset);
        
        % get parameter values
        cs_prb = cs_prb_set(iset);
        cs_1st = cs_1st_set(iset);
        sd_inf = sd_inf_set(iset);
        el_inf = el_inf_set(iset);
        th_cnf = th_cnf_set(iset);
        sd_obs = sd_obs_set(iset);
        
        if fit_sd_inf
            % use fitted inference noise
            sd_inf = sd_inf_ref;
        end

        comd = nan(nseq,12);
        comt = false(nseq,12);
        xsmp = nan(nseq,12);
        xobs = nan(nseq,12);
        xacc = nan(nseq,12);
        xtot = nan(nseq,1);
        resp = nan(nseq,1);
        comf = nan(nseq,1);
        
        rcnf = nan(nseq,12);
        
        xcons_tru = nan(nseq,12);
        xconf_tru = nan(nseq,12);
        xcons_obs = nan(nseq,12);
        xconf_obs = nan(nseq,12);
        xprep_obs = nan(nseq,12);
        
        for iseq = 1:nseq
            
            % objective stimulus-wise evidence
            seqllr = sin(2*(seqang{iseq}-catdir(iseq)));
            
            cd = nan;
            ct = false;
            
            nsmp = numel(seqllr);
            for ismp = 1:nsmp
                
                % pre-commit
                if ~ct && ismp >= cs_1st && rand < cs_prb
                    cd = sign(xacc(iseq,ismp-1));
                    ct = true;
                end
                comd(iseq,ismp) = cd;
                comt(iseq,ismp) = ct;
                
                % estimate evidence
                if ~ct || (ct && sign(seqllr(ismp)) == cd)
                    % no pre-commitment or consistent stimulus
                    xsmp(iseq,ismp) = seqllr(ismp)+sd_inf*randn;
                else
                    % pre-commitment and conflicting stimulus
                    xsmp(iseq,ismp) = sd_inf*randn;
                end
                
                % observe evidence
                xobs(iseq,ismp) = abs(xsmp(iseq,ismp))+sd_obs*randn;
                
                % accumulate evidence
                if ismp == 1
                    xacc(iseq,1) = xsmp(iseq,1);
                else
                    xacc(iseq,ismp) = exp(-el_inf)*xacc(iseq,ismp-1)+xsmp(iseq,ismp);
                end
                
            end
            
            % decide
            if ~ct
                % no pre-commitment => choose based on accumulated evidence
                resp(iseq) = sign(xacc(iseq,nsmp));
            else
                % pre-commitment => choose based on pre-commitment direction
                resp(iseq) = cd;
            end
            
            % total accumulated evidence in favor of provided response
            xtot(iseq,1) = xacc(iseq,nsmp)*resp(iseq);
            
            % has pre-committed before end of sequence?
            comf(iseq,1) = comt(iseq,nsmp);
            
            % sort stimuli wrt consistency with upcoming decision
            icons = sign(seqllr) == sign(resp(iseq));
            iconf = ~icons;
            xcons_tru(iseq,icons) = abs(seqllr(icons));
            xconf_tru(iseq,iconf) = abs(seqllr(iconf));
            xcons_obs(iseq,icons) = xobs(iseq,icons);
            xconf_obs(iseq,iconf) = xobs(iseq,iconf);
            
            % response preparation signal (only if pre-committed)
            xprep_obs(iseq,:) = comt(iseq,:).*xacc(iseq,:)*resp(iseq)+sd_obs*2*randn(1,12);
            
            % make confidence report (not compared to lottery here)
            rcnf(iseq,1:nsmp) = abs(xacc(iseq,1:nsmp)) > th_cnf;
            
        end
        
        seqlen = cellfun(@length,seqang);
        
        % get overall fraction confident
        pcnf(isubj,:,iset) = mean(rcnf(seqlen == 12,:),1);
        
        % get conditional fractions confident (with/without pre-commitment)
        pcnf_comf(isubj,:,iset) = mean(rcnf(seqlen(:) == 12 & comf(:) == true,:),1);
        pcnf_comt(isubj,:,iset) = mean(rcnf(seqlen(:) == 12 & comf(:) == false,:),1);

        % compute confidence sensitivity
        for ipos = 1:12
            ifilt = comf(:) == false & seqlen(:) == 12;
            aroc_comf(isubj,ipos,iset) = ...
                get_aroc(1+(xacc(ifilt,ipos) < 0) == seqcat(ifilt)',abs(xacc(ifilt,ipos)));
            ifilt = comf(:) == true & seqlen(:) == 12;
            if nnz(ifilt) > 0
                aroc_comt(isubj,ipos,iset) = ...
                    get_aroc(1+(xacc(ifilt,ipos) < 0) == seqcat(ifilt)',abs(xacc(ifilt,ipos)));
            end
        end
        
        seqlen = cellfun(@length,seqang);
        iseq = comf == false & seqlen(:) == 12;
        xaccu_comf(isubj,:,1,iset) = mean(abs(xacc(iseq,:)),1);
        iseq = comf == true & seqlen(:) == 12;
        xaccu_comf(isubj,:,2,iset) = mean(abs(xacc(iseq,:)),1);
        
        if iset == 3
            xtot_com = cat(1,xtot_com,xtot( any(comt,2)));
            xtot_reg = cat(1,xtot_reg,xtot(~any(comt,2)));
        end
        
        % all trials
        for ipos = 1:6
            % group stimuli by pairs
            jpos = 2*(ipos-1)+(1:2);
            % correlate observed and true evidence values for consistent stimuli
            xt = reshape(xcons_tru(:,jpos),[],1);
            xo = reshape(xcons_obs(:,jpos),[],1);
            zdec_cons(isubj,ipos,1,iset) = r2z(corr(xt,xo,'rows','complete'));
            zerr_cons(isubj,ipos,1,iset) = 1/sqrt(numel(xt)/nsim-3);
            % correlate observed and true evidence values for conflicting stimuli
            xt = reshape(xconf_tru(:,jpos),[],1);
            xo = reshape(xconf_obs(:,jpos),[],1);
            zdec_cons(isubj,ipos,2,iset) = r2z(corr(xt,xo,'rows','complete'));
            zerr_cons(isubj,ipos,2,iset) = 1/sqrt(numel(xt)/nsim-3);
        end
        for j = 1:2
            x = linspace(0,1,6)';
            y = zdec_cons(isubj,:,j,iset)';
            b = [x,ones(6,1)]\y;
            zchg_cons(isubj,j,iset) = b(1);
        end
        xresp(isubj,:,iset) = nanmean(xprep_obs,1);
        eresp(isubj,:,iset) = nanstd(xprep_obs,[],1)./sqrt(sum(~isnan(xprep_obs),1)/nsim);

        % trials with less accumulated evidence
        iseq = xtot < median(xtot);
        for ipos = 1:6
            % group stimuli by pairs
            jpos = 2*(ipos-1)+(1:2);
            indx = sub2ind([nseq,12],repmat(find(iseq),[2,1]),kron(jpos(:),ones(nnz(iseq),1)));
            % correlate observed and true evidence values for consistent stimuli
            xt = reshape(xcons_tru(indx),[],1);
            xo = reshape(xcons_obs(indx),[],1);
            zdec_cons_less(isubj,ipos,1,iset) = r2z(corr(xt,xo,'rows','complete'));
            zerr_cons_less(isubj,ipos,1,iset) = 1/sqrt(nnz(~isnan(xt(:)))/nsim-3);
            % correlate observed and true evidence values for conflicting stimuli
            xt = reshape(xconf_tru(indx),[],1);
            xo = reshape(xconf_obs(indx),[],1);
            zdec_cons_less(isubj,ipos,2,iset) = r2z(corr(xt,xo,'rows','complete'));
            zerr_cons_less(isubj,ipos,2,iset) = 1/sqrt(nnz(~isnan(xt(:)))/nsim-3);
        end
        xresp_less(isubj,:,iset) = nanmean(xprep_obs(iseq,:),1);
        eresp_less(isubj,:,iset) = nanstd(xprep_obs(iseq,:),[],1)./sqrt(sum(~isnan(xprep_obs(iseq,:)),1)/nsim);
        comf_less(isubj,iset) = mean(comf(iseq));
        xt = reshape(xcons_tru(iseq,:),[],1);
        xo = reshape(xcons_obs(iseq,:),[],1);
        zdec_cons_less_all(isubj,1,iset) = r2z(corr(xt,xo,'rows','complete'));
        zerr_cons_less_all(isubj,1,iset) = 1/sqrt(nnz(~isnan(xt(:)))/nsim-3);
        xt = reshape(xconf_tru(iseq,:),[],1);
        xo = reshape(xconf_obs(iseq,:),[],1);
        zdec_cons_less_all(isubj,2,iset) = r2z(corr(xt,xo,'rows','complete'));
        zerr_cons_less_all(isubj,2,iset) = 1/sqrt(nnz(~isnan(xt(:)))/nsim-3);
        
        % trials with more accumulated evidence
        iseq = xtot > median(xtot);
        for ipos = 1:6
            % group stimuli by pairs
            jpos = 2*(ipos-1)+(1:2);
            indx = sub2ind([nseq,12],repmat(find(iseq),[2,1]),kron(jpos(:),ones(nnz(iseq),1)));
            % correlate observed and true evidence values for consistent stimuli
            xt = reshape(xcons_tru(indx),[],1);
            xo = reshape(xcons_obs(indx),[],1);
            zdec_cons_more(isubj,ipos,1,iset) = r2z(corr(xt,xo,'rows','complete'));
            zerr_cons_more(isubj,ipos,1,iset) = 1/sqrt(nnz(~isnan(xt(:)))/nsim-3);
            % correlate observed and true evidence values for conflicting stimuli
            xt = reshape(xconf_tru(indx),[],1);
            xo = reshape(xconf_obs(indx),[],1);
            zdec_cons_more(isubj,ipos,2,iset) = r2z(corr(xt,xo,'rows','complete'));
            zerr_cons_more(isubj,ipos,2,iset) = 1/sqrt(nnz(~isnan(xt(:)))/nsim-3);
        end
        xresp_more(isubj,:,iset) = nanmean(xprep_obs(iseq,:),1);
        eresp_more(isubj,:,iset) = nanstd(xprep_obs(iseq,:),[],1)./sqrt(sum(~isnan(xprep_obs(iseq,:)),1)/nsim);
        comf_more(isubj,iset) = mean(comf(iseq));
        xt = reshape(xcons_tru(iseq,:),[],1);
        xo = reshape(xcons_obs(iseq,:),[],1);
        zdec_cons_more_all(isubj,1,iset) = r2z(corr(xt,xo,'rows','complete'));
        zerr_cons_more_all(isubj,1,iset) = 1/sqrt(nnz(~isnan(xt(:)))/nsim-3);
        xt = reshape(xconf_tru(iseq,:),[],1);
        xo = reshape(xconf_obs(iseq,:),[],1);
        zdec_cons_more_all(isubj,2,iset) = r2z(corr(xt,xo,'rows','complete'));
        zerr_cons_more_all(isubj,2,iset) = 1/sqrt(nnz(~isnan(xt(:)))/nsim-3);

        % define sequences to be used for model fit
        iseq = 1:min(nseq,nseq/nsim*nfit);
        % fit model
        cfg        = [];
        cfg.stype  = stype(iseq);
        cfg.seqang = seqang(iseq);
        cfg.catdir = catdir(iseq);
        cfg.resp   = 1+(resp(iseq) < 0);
        cfg.parmod = {};
        cfg.ec_inf = 0;
        cfg.criter = 0;
        out_fit = fit_model_resp(cfg);
        % store best-fitting parameter values
        sd_inf_hat(isubj,iset) = out_fit.sd_inf; 
        el_inf_hat(isubj,iset) = out_fit.el_inf;
        
    end

end

%% Plot inference noise across parameter sets
%
%  This measures the effect of pre-commitment probability on the effective
%  inference noise measured using the computational model of behavior.

% compute group-level statistics of best-fitting inference noise
zavg = mean(sd_inf_hat,1);
zstd = std(sd_inf_hat,[],1);

% define colors
rgb = [1,0,0.2];
rgb = cat(1,rgb*0.25+0.75,rgb*0.50+0.50,rgb*0.75+0.25,rgb);

% plot inference noise across parameter sets
pbar = 1/2*(4.2/2.2);
figure;
hold on
xlim([0.4,4.6]);
ylim([0,1]);
for icond = 1:4
    pos = icond;
    wid = 0.4;
    xk = linspace(zavg(icond)-2*zstd(icond),zavg(icond)+2*zstd(icond),100);
    pk = normpdf(xk,zavg(icond),zstd(icond));
    pk = pk/max(pk);
    str = interp1(xk,pk,x);
    jit = linspace(-1,+1,numel(x));
    patch([pos+pk*wid,fliplr(pos-pk*wid)],[xk,fliplr(xk)],0.5*(rgb(icond,:)+1),'EdgeColor','none');
    plot(pos+pk*wid,xk,'-','Color',rgb(icond,:),'LineWidth',0.5);
    plot(pos-pk*wid,xk,'-','Color',rgb(icond,:),'LineWidth',0.5);
end
plot(xlim,[0,0],'k-');
plot(xlim,zavg(1)*[1,1],'b-');
for icond = 1:4
    pos = icond;
    plot(pos,zavg(icond),'ko','MarkerFaceColor','k','MarkerSize',5,'LineWidth',0.5);
end
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',[1,2,3,4],'XTickLabel',{'0','0.05','0.10','0.15'});
set(gca,'YTick',0:0.2:1);
xlabel('p(commitment)','FontSize',8);
ylabel('inference noise','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

%% Plot coding precision and change as a function of stimulus position

% set parameters
iset = 3; % parameter set of interest
k    = 1; % 1=trials without pre-commitment, 2=trials with pre-commitment

% define colors
rgb = [1,0,0.2];
rgb = cat(1,rgb*0.25+0.75,rgb*0.50+0.50,rgb*0.75+0.25,rgb);

% plot coding precision as a function of stimulus position
pbar = 1.25;
figure('Color','white');
xlim([0.5,6.5]);
ylim([-1/5,1]*0.2);
hold on
for iset = 1:nset
    zavg = mean(zdec_cons(:,:,k,iset),1);
    zerr = mean(zerr_cons(:,:,k,iset),1);
    patch([1:6,fliplr(1:6)],[zavg+zerr,fliplr(zavg-zerr)],rgb(iset,:)*0.50+0.50,'EdgeColor','none');
end
plot(xlim,[0,0],'k-','Color',rand(1,3));
for iset = 1:nset
    zavg = mean(zdec_cons(:,:,k,iset),1);
    plot(1:6,zavg,'-','Color',rgb(iset,:),'LineWidth',2);
    plot(1:6,zavg,'ko','MarkerFaceColor',rgb(iset,:),'MarkerSize',5);
end
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',1:6,'XTickLabel',{'2','4','6','8','10','12'});
set(gca,'YTick',0:0.1:1);
set(gca,'YTickLabel',arrayfun(@(x)num2str(x,'%.1f'),get(gca,'YTick'),'UniformOutput',false));
xlabel('stimulus position','FontSize',8);
ylabel('coding precision (z)','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

% plot coding change as a function of stimulus position
zavg = squeeze(mean(zchg_cons,1));
zstd = squeeze(std(zchg_cons,[],1))*sqrt(nsim);
pbar = 1/2*(4.2/2.2);
figure;
hold on
xlim([0.4,4.6]);
ylim([-0.4,+0.4]);
for icond = 1:4
    pos = icond;
    wid = 0.4;
    xk = linspace(zavg(k,icond)-2*zstd(k,icond),zavg(k,icond)+2*zstd(k,icond),100);
    pk = normpdf(xk,zavg(k,icond),zstd(k,icond));
    pk = pk/max(pk);
    str = interp1(xk,pk,x);
    jit = linspace(-1,+1,numel(x));
    patch([pos+pk*wid,fliplr(pos-pk*wid)],[xk,fliplr(xk)],0.5*(rgb(icond,:)+1),'EdgeColor','none');
    plot(pos+pk*wid,xk,'-','Color',rgb(icond,:),'LineWidth',0.5);
    plot(pos-pk*wid,xk,'-','Color',rgb(icond,:),'LineWidth',0.5);
end
plot(xlim,[0,0],'k-');
plot(1:4,zavg(k,:),'-','LineWidth',1,'Color','k');
for icond = 1:4
    pos = icond;
    plot(pos,zavg(k,icond),'ko','MarkerFaceColor','k','MarkerSize',5,'LineWidth',0.5);
end
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',[1,2,3,4],'XTickLabel',{'0','0.05','0.10','0.15'});
set(gca,'YTick',-0.4:0.2:+0.4,'YTickLabel',{'0.4','0.2','0','0.2','0.4'});
xlabel('p(commitment)','FontSize',8);
ylabel('coding change (z)','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

%% Plot coding precision for consistent and conflicting stimuli
%  left=trials which ended with more decision evidence (validation more likely)
%  right=trials which ended with less decision evidence (opt-out more likely)

iset = 1; % parameter set of interest (1=without pre-commitment)

% get group-level variables of interest
zl = squeeze(mean(zdec_cons_less_all,1));
el = squeeze(mean(zerr_cons_less_all,1));
zm = squeeze(mean(zdec_cons_more_all,1));
em = squeeze(mean(zerr_cons_more_all,1));

% define colors
rgb = [0.5,0.5,0.5;1.0,0.0,0.2];
if iset == 1
    rgb = [rgb(1,:);rgb(1,:)*0.5+0.5];
else
    rgb = [rgb(2,:);rgb(2,:)*0.5+0.5];
end

% plot coding precision for consistent and conflicting stimuli
pbar = 1/2*4.2/2.2;
bavg = [zm(:,iset);zl(:,iset)];
berr = [em(:,iset);el(:,iset)];
figure('Color','white');
hold on
xlim([0.4,4.6]);
ylim([0,0.2]);
for i = 1:4
    k = mod(i-1,2)+1;
    bar(i,bavg(i),0.8,'FaceColor',0.5*(rgb(k,:)+1),'LineWidth',1,'BaseValue',0,'EdgeColor',rgb(k,:));
end
for i = 1:4
    plot(i*[1,1],bavg(i)+berr(i)*[-1,+1],'k-');
end
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',1:4,'XTickLabel',{'cons.','conf.','cons.','conf.'});
set(gca,'YTick',0:0.1:1);
xlabel('session','FontSize',8);
ylabel('precision (z)','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

% get group-level variables of interest
zdif = zdec_cons_more_all(:,1,[1,3])-zdec_cons_more_all(:,2,[1,3]);
zerr = sqrt(zerr_cons_more_all(:,1,[1,3]).^2+zerr_cons_more_all(:,2,[1,3]).^2);
zm = squeeze(mean(zdif,1));
ze = squeeze(mean(zerr,1));
zdif = zdec_cons_less_all(:,1,[1,3])-zdec_cons_less_all(:,2,[1,3]);
zerr = sqrt(zerr_cons_less_all(:,1,[1,3]).^2+zerr_cons_less_all(:,2,[1,3]).^2);
zm = cat(2,zm,squeeze(mean(zdif,1)));
ze = cat(2,ze,squeeze(mean(zerr,1)));
zm = zm';
ze = ze';

% define colors
rgb = [0.5,0.5,0.5;1.0,0.0,0.2];
rgb = [rgb(1,:);rgb(1,:)*0.5+0.5;rgb(2,:);rgb(2,:)*0.5+0.5];

% plot consistency effect on coding precision
pbar = 1/2*(4.2/2.2);
figure;
hold on
xlim([0.4,4.6]);
ylim([-0.15,+0.15]);
for icond = 1:4
    pos = icond;
    wid = 0.4;
    xk = linspace(zm(icond)-2*ze(icond),zm(icond)+2*ze(icond),100);
    pk = normpdf(xk,zm(icond),ze(icond));
    pk = pk/max(pk);
    str = interp1(xk,pk,x);
    jit = linspace(-1,+1,numel(x));
    patch([pos+pk*wid,fliplr(pos-pk*wid)],[xk,fliplr(xk)],0.5*(rgb(icond,:)+1),'EdgeColor','none');
    plot(pos+pk*wid,xk,'-','Color',rgb(icond,:),'LineWidth',0.5);
    plot(pos-pk*wid,xk,'-','Color',rgb(icond,:),'LineWidth',0.5);
end
plot(xlim,[0,0],'k-');
plot(xlim,zavg(1)*[1,1],'b-');
for icond = 1:4
    pos = icond;
    plot(pos,zm(icond),'ko','MarkerFaceColor','k','MarkerSize',5,'LineWidth',0.5);
end
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',[1,2,3,4],'XTickLabel',{'a','b','c','d'});
set(gca,'YTick',-0.4:0.1:+0.4);
xlabel('condition','FontSize',8);
ylabel('coding change','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

%% Plot fraction of pre-committed trials as a function of decision evidence

% define colors
rgb = [0.5,0.5,0.5;1.0,0.0,0.2];
rgb = [rgb(2,:);rgb(2,:)*0.5+0.5];

% plot fraction committed for trials with more or less decision evidence
pbar = 1/2;
bavg = [mean(comf_more(:,3)),mean(comf_less(:,3))];
berr = [std(comf_more(:,3)),std(comf_less(:,3))];
figure('Color','white');
hold on
xlim([0.4,2.6]);
ylim([0,0.5]);
for i = 1:2
    bar(i,bavg(i),0.8,'FaceColor',0.5*(rgb(i,:)+1),'LineWidth',1,'BaseValue',0,'EdgeColor',rgb(i,:));
end
for i = 1:2
    plot(i*[1,1],bavg(i)+berr(i)*[-1,+1],'k-');
end
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',1:2,'XTickLabel',{'strong','weak'});
set(gca,'YTick',0:0.1:1);
xlabel('belief','FontSize',8);
ylabel('fraction committed','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

%% Plot time course of belief strength for trials with/without pre-commitments

% plot absolute belief strength
rgb = [0.5,0.5,0.5;1.0,0.0,0.2];
pbar = 1.15;
figure;
hold on
xlim([0.5,12.5]);
ylim([-0.4,+0.4]);
plot(xlim,[0,0]);
plot(squeeze(mean(xaccu_comf(:,:,1,3)-0.2,1))-squeeze(mean(xaccu_comf(:,:,1,1),1)),'k-');
plot(squeeze(mean(xaccu_comf(:,:,2,3)-0.2,1))-squeeze(mean(xaccu_comf(:,:,1,1),1)),'k-');
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',[1],'XTickLabel',{'start'});
set(gca,'YTick',-0.4:0.2:+0.4,'YTickLabel',{'0.4','0.2','0','0.2','0.4'});
xlabel('number of presented cues','FontSize',8);
ylabel('relative belief strength','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

% plot relative belief strength
rgb = [0.5,0.5,0.5;1.0,0.0,0.2];
pbar = 1.15;
figure;
hold on
xlim([0.5,12.5]);
ylim([0,1.5]);
plot(squeeze(mean(xaccu_comf(:,:,1,1),1)),'Color',rgb(1,:),'LineWidth',2);
plot(squeeze(mean(xaccu_comf(:,:,1,3)-0.2,1)),'Color',rgb(2,:),'LineWidth',1);
plot(squeeze(mean(xaccu_comf(:,:,2,3)-0.2,1)),'Color',rgb(2,:),'LineWidth',2);
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',[1],'XTickLabel',{'start'});
set(gca,'YTick',0:0.5:1.5,'YTickLabel',{'0','0.5','1.0','1.5'});
xlabel('number of presented cues','FontSize',8);
ylabel('belief strength','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

%% Plot time course of metacognitive sensitivity

% define parameter set of interest
iset = 3;

% plot metacognitive sensitivity (AROC) as a function of number of cues
rgb = [0.5,0.5,0.5;1.0,0.0,0.2];
pbar = 1.15;
figure;
hold on
xlim([0.5,12.5]);
ylim([0.5,0.75]);
plot(xlim,[0,0]);
plot(mean(aroc_comf(:,:,iset),1),'-','Color',rgb(1,:),'LineWidth',2); % trials without pre-commitment
plot(mean(aroc_comt(:,:,iset),1),'-','Color',rgb(2,:),'LineWidth',2); % trials with pre-commitment
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',[1],'XTickLabel',{'start'});
set(gca,'YTick',0.5:0.1:1);
xlabel('number of presented cues','FontSize',8);
ylabel('metacognitive sensitivity (AROC)','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end
