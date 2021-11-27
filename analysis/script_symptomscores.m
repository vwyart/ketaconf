
%% Analyze symptom scores and ketamine plasma concentration for KETACONF study
%
%  Data are available upon request.
%
%  Valentin Wyart <valentin.wyart@inserm.fr>

clear all
close all
clc

% get list of subjects
subjlist = get_sesstype(1,true);
nsubj = numel(subjlist);

% data available upon request
load ./dat_symptomscores.mat

% read data
t = [0,30,90];
for isubj = 1:nsubj
    subj = 100+subjlist(isubj);
    for stype = 1:2
        for it = 1:3
            i = dat.subj == subj & dat.stype == stype & dat.time == t(it);
            j = i & strcmp(dat.snam,'BPRS');
            vBPRS(isubj,it,stype) = dat.sval(j);
            j = i & strcmp(dat.snam,'CADSS');
            vCADSS(isubj,it,stype) = dat.sval(j);
            vketa(isubj,it,stype) = dat.keta(j);
        end
    end
end

% select ketamine session
stype = 2;

% define colors
rgb = [0.5,0.5,0.5;1.0,0.0,0.2];
rgb = [0.33*rgb(stype,:)+0.67;rgb(stype,:);0.67*rgb(stype,:)+0.33];

% plot BPRS scores
pbar = 1/2*3.2/2.2;
figure;
hold on
xlim([0.4,3.6]);
ylim([22.5,42.5]);
for icond = 1:3
    x = vBPRS(:,icond,stype);
    pos = icond;
    wid = 0.4;
    xk = linspace(min(x),max(x),100);
    pk = ksdensity(x,xk);
    pk = pk/max(pk);
    str = interp1(xk,pk,x);
    jit = linspace(-1,+1,numel(x));
    patch([pos+pk*wid,fliplr(pos-pk*wid)],[xk,fliplr(xk)],0.5*(rgb(icond,:)+1),'EdgeColor','none');
    plot(pos+pk*wid,xk,'-','Color',rgb(icond,:),'LineWidth',0.5);
    plot(pos-pk*wid,xk,'-','Color',rgb(icond,:),'LineWidth',0.5);
end
plot(xlim,24*[1,1],'b-');
for icond = 1:3
    x = vBPRS(:,icond,stype);
    pos = icond;
    wid = 0.4;
    xk = linspace(min(x),max(x),100);
    pk = ksdensity(x,xk);
    pk = pk/max(pk);
    str = interp1(xk,pk,x);
    jit = linspace(-1,+1,numel(x));
    for i = 1:numel(x)
        plot(pos+jit(i)*str(i)*wid,x(i),'wo','MarkerFaceColor',rgb(icond,:),'MarkerSize',4,'LineWidth',0.5);
        xsub(i,icond) = pos+jit(i)*str(i)*wid;
    end
    plot(pos*[1,1],mean(x)+std(x)/sqrt(numel(x))*[-1,+1],'k-');
    plot(pos,mean(x),'ko','MarkerFaceColor','k','MarkerSize',5,'LineWidth',0.5);
end
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',[1,2,3],'XTickLabel',{'0','30','90'});
set(gca,'YTick',0:5:40);
xlabel('time','FontSize',8);
ylabel('BPRS score','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

% plot CADSS scores
pbar = 1/2*3.2/2.2;
figure;
hold on
xlim([0.4,3.6]);
ylim([-5,60]);
for icond = 1:3
    x = vCADSS(:,icond,stype);
    pos = icond;
    wid = 0.4;
    xk = linspace(min(x),max(x),100);
    pk = ksdensity(x,xk);
    pk = pk/max(pk);
    str = interp1(xk,pk,x);
    jit = linspace(-1,+1,numel(x));
    patch([pos+pk*wid,fliplr(pos-pk*wid)],[xk,fliplr(xk)],0.5*(rgb(icond,:)+1),'EdgeColor','none');
    plot(pos+pk*wid,xk,'-','Color',rgb(icond,:),'LineWidth',0.5);
    plot(pos-pk*wid,xk,'-','Color',rgb(icond,:),'LineWidth',0.5);
end
plot(xlim,0*[1,1],'b-');
for icond = 1:3
    x = vCADSS(:,icond,stype);
    pos = icond;
    wid = 0.4;
    xk = linspace(min(x),max(x),100);
    pk = ksdensity(x,xk);
    pk = pk/max(pk);
    str = interp1(xk,pk,x);
    jit = linspace(-1,+1,numel(x));
    for i = 1:numel(x)
        plot(pos+jit(i)*str(i)*wid,x(i),'wo','MarkerFaceColor',rgb(icond,:),'MarkerSize',4,'LineWidth',0.5);
        xsub(i,icond) = pos+jit(i)*str(i)*wid;
    end
    plot(pos*[1,1],mean(x)+std(x)/sqrt(numel(x))*[-1,+1],'k-');
    plot(pos,mean(x),'ko','MarkerFaceColor','k','MarkerSize',5,'LineWidth',0.5);
end
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',[1,2,3],'XTickLabel',{'0','30','90'});
set(gca,'YTick',0:20:60);
xlabel('time','FontSize',8);
ylabel('CADSS score','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

% plot ketamine plasma concentrations
xall = [];
pbar = 1/2;
figure;
hold on
xlim([0.4,2.6]);
ylim([0,300]);
for icond = 1:2
    x = vketa(:,icond+1,stype);
    xall(:,icond) = x;
    pos = icond;
    wid = 0.4;
    xk = linspace(min(x),max(x),100);
    pk = ksdensity(x,xk);
    pk = pk/max(pk);
    str = interp1(xk,pk,x);
    jit = linspace(-1,+1,numel(x));
    patch([pos+pk*wid,fliplr(pos-pk*wid)],[xk,fliplr(xk)],0.5*(rgb(icond+1,:)+1),'EdgeColor','none');
    plot(pos+pk*wid,xk,'-','Color',rgb(icond+1,:),'LineWidth',0.5);
    plot(pos-pk*wid,xk,'-','Color',rgb(icond+1,:),'LineWidth',0.5);
    for i = 1:numel(x)
        plot(pos+jit(i)*str(i)*wid,x(i),'wo','MarkerFaceColor',rgb(icond+1,:),'MarkerSize',4,'LineWidth',0.5);
    end
    plot(pos*[1,1],mean(x)+std(x)/sqrt(numel(x))*[-1,+1],'k-');
    plot(pos,mean(x),'ko','MarkerFaceColor','k','MarkerSize',5,'LineWidth',0.5);
end
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',[1,2],'XTickLabel',{'30','90'});
set(gca,'YTick',0:100:300);
xlabel('time','FontSize',8);
ylabel('plasma concentration','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

%% Plot pharmacokinetics time course

t = -30:0.01:100;
k = 1; % do not change!

% define ketamine infusion (mg/kg/min)
v = zeros(size(t));
v(t >  0 & t <= 10) = 0.23*k/10;
v(t > 10 & t <= 30) = 0.00967*k;
v(t > 30 & t <= 90) = 0.00483*k;

% define colors
rgb = [0.5,0.5,0.5;1.0,0.0,0.2];

% plot pharmacokinetics time course
pbar = 3;
figure;
hold on
xlim([-30,+100]);
ylim([0,0.025]);
i = t >= 0 & t <= 10;
patch([t(i),fliplr(t(i))],[v(i),zeros(size(v(i)))],0.50*(rgb(2,:)+1),'EdgeColor','none');
i = t >= 10 & t <= 30;
patch([t(i),fliplr(t(i))],[v(i),zeros(size(v(i)))],0.33*rgb(2,:)+0.67,'EdgeColor','none');
i = t >= 30 & t <= 90;
patch([t(i),fliplr(t(i))],[v(i),zeros(size(v(i)))],0.17*rgb(2,:)+0.83,'EdgeColor','none');
plot(t,v,'-','Color',rgb(2,:),'LineWidth',2);
hold off
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
set(gca,'XTick',0:30:90);
set(gca,'YTick',0:0.01:0.02);
xlabel('time (min)','FontSize',8);
ylabel('ketamine infusion (mg/kg/min)','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

%% Get subject information

% get subject information
sage = []; % subject age
ssiz = []; % subject height (cm)
ssex = []; % subject gender (1:male 2:female)
swgh = []; % subject weight (kg)
for isubj = 1:nsubj
    subj = 100+subjlist(isubj);
    i = find(dat.subj == subj,1);
    sage(isubj,1) = dat.age(i);
    ssiz(isubj,1) = dat.size(i);
    ssex(isubj,1) = dat.sex(i);
    swgh(isubj,1) = dat.weight(i);
end
