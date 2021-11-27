function [img] = make_pie(cfg)
%  MAKE_PIE  Make pie chart for lottery
%
%  Usage: [img] = MAKE_PIE(cfg)
%
%  This worker function should not be called directly unless you know what you
%  are doing. See KETACONF_run_expe.m for examples.
%
%  Valentin Wyart <valentin.wyart@ens.fr>

% check configuration structure
if ~all(isfield(cfg,{'ppd','pfill'}))
    error('invalid configuration structure!');
end

% get required configuration parameters
ppd   = cfg.ppd;
pfill = cfg.pfill;

% set internal configuration parameters
carddmtr = deg2pix(8,2);
cardwdth = deg2pix(8/60,1);
axiswdth = deg2pix(4/60,1);
lumibg   = 128/255;
lumifill = [112,160]/255;

% define function handles
fsmooth = @(x,r)r(1)+diff(r)./(1+99.^(-x)); % smoothing function
modzero = @(x,m)mod(x+m/2,m)-m/2; % zero-centered modulus

% set outer diameter
outrdmtr = carddmtr+4*cardwdth;

% set polar coordinates
[x,y] = meshgrid([1:outrdmtr]-(outrdmtr+1)/2);
r = sqrt(x.^2+y.^2);
t = -atan2(y,x);
t = mod(t+pi/2+pi*pfill,2*pi);

% start with background patch
patchfg = lumibg(ones(outrdmtr,outrdmtr));

% draw darker/filled portion of pie
patchtp = ones(outrdmtr,outrdmtr);
patchtp = min(patchtp,fsmooth(r-carddmtr/2+2,[1,0]));
patchtp = min(patchtp,fsmooth(carddmtr*(t-2*pi*pfill),[1,0]));
patchfg = lumifill(1)*patchtp+patchfg.*(1-patchtp);

% draw lighter/empty portion of pie
patchtp = ones(outrdmtr,outrdmtr);
patchtp = min(patchtp,fsmooth(r-carddmtr/2+2,[1,0]));
patchtp = min(patchtp,fsmooth(carddmtr*(t-2*pi*pfill),[0,1]));
patchfg = lumifill(2)*patchtp+patchfg.*(1-patchtp);

% draw filling borders
% 1st border
axis = pi/2-pi*pfill;
u = sin(axis)*x+cos(axis)*y;
v = cos(axis)*x-sin(axis)*y;
patchtp = ones(outrdmtr,outrdmtr);
patchtp = min(patchtp,fsmooth(u+axiswdth/2,[0,1]));
patchtp = min(patchtp,fsmooth(u-axiswdth/2,[1,0]));
patchtp = min(patchtp,fsmooth(v,[1,0]));
patchtp = min(patchtp,fsmooth(r-carddmtr/2+2,[1,0]));
patchfg = bsxfun(@times,patchfg,1-patchtp);
% 2nd border
axis = pi/2+pi*pfill;
u = sin(axis)*x+cos(axis)*y;
v = cos(axis)*x-sin(axis)*y;
patchtp = ones(outrdmtr,outrdmtr);
patchtp = min(patchtp,fsmooth(u+axiswdth/2,[0,1]));
patchtp = min(patchtp,fsmooth(u-axiswdth/2,[1,0]));
patchtp = min(patchtp,fsmooth(v,[1,0]));
patchtp = min(patchtp,fsmooth(r-carddmtr/2+2,[1,0]));
patchfg = bsxfun(@times,patchfg,1-patchtp);

% draw card border
patchtp = ones(outrdmtr,outrdmtr);
patchtp = min(patchtp,fsmooth(r-carddmtr/2+cardwdth+2,[0,1]));
patchtp = min(patchtp,fsmooth(r-carddmtr/2+2,[1,0]));
patchfg = bsxfun(@times,patchfg,1-patchtp);

% return image
img = patchfg;

    function [p] = deg2pix(d,b)
        p = d*ppd;
        if nargin > 1 && b > 0
            p = max(round(p/b),1)*b;
        end
    end

end