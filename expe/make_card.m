function [img] = make_card(cfg)
%  MAKE_CARD  Make card
%
%  Usage: [img] = MAKE_CARD(cfg)
%
%  This worker function should not be called directly unless you know what you
%  are doing. See KETACONF_run_expe.m for examples.
%
%  Valentin Wyart <valentin.wyart@ens.fr>

% check configuration structure
if ~all(isfield(cfg,{'ppd','cardtype'}))
    error('invalid configuration structure!');
end
if ~isfield(cfg,'angl')
    cfg.angl = [];
end
if ~isfield(cfg,'axis')
    cfg.axis = [];
end

% get required configuration parameters
ppd      = cfg.ppd;
cardtype = cfg.cardtype;
angl     = cfg.angl;
axis     = cfg.axis;

% set internal configuration parameters
carddmtr = deg2pix(8,2);
cardwdth = deg2pix(8/60,1);
anglwdth = deg2pix(8/60,1);
axiswdth = deg2pix(4/60,1);
lumibg   = 128/255;
showaxis = false;

% set type-dependent configuration parameters
switch cardtype
    case 0 % empty card
        axis     = 0;
        angl     = 0;
        lumifg   = 160/255;
        opctfg   = 0;
        opctangl = 0;
    case 1 % mapping/response card
        if isempty(axis)
            error('invalid configuration structure!');
        end
        angl     = 0;
        lumifg   = 160/255;
        opctfg   = 2/3;
        opctangl = 0;
    case 2 % flashing/blinking card
        if isempty(axis)
            error('invalid configuration structure!');
        end
        angl     = 0;
        lumifg   = 176/255;
        opctfg   = 2/3;
        opctangl = 0;
    case 3 % angle card
        if isempty(axis) || isempty(angl)
            error('invalid configuration structure!');
        end
        lumifg   = 160/255;
        opctfg   = 2/3;
        opctangl = 1;
    otherwise
        error('undefined card type!');
end

% define function handles
fsmooth = @(x,r)r(1)+diff(r)./(1+99.^(-x)); % smoothing function
modzero = @(x,m)mod(x+m/2,m)-m/2; % zero-centered modulus

% set outer diameter
outrdmtr = carddmtr+4*cardwdth;

% set polar coordinates
[x,y] = meshgrid([1:outrdmtr]-(outrdmtr+1)/2);
r = sqrt(x.^2+y.^2);
t = -atan2(y,x);

% draw card foreground
colrfg = [1,lumifg,lumifg*2-1;lumifg*2-1,lumifg,1];
p = 0.5*(1+sin(2*(t-axis)));
patchfg = ...
    bsxfun(@times,  p,reshape(colrfg(1,:),1,1,[]))+ ...
    bsxfun(@times,1-p,reshape(colrfg(2,:),1,1,[]));
patchfg = opctfg*patchfg+(1-opctfg)*lumifg;
patchtp = ones(outrdmtr,outrdmtr);
patchtp = min(patchtp,fsmooth(r-carddmtr/2+2,[1,0]));
patchfg = bsxfun(@plus,bsxfun(@times,patchfg,patchtp),lumibg*(1-patchtp));

if showaxis % draw category boundaries?
    % primary boundary
    u = sin(axis)*x+cos(axis)*y;
    patchtp = ones(outrdmtr,outrdmtr);
    patchtp = min(patchtp,fsmooth(u+axiswdth/2,[0,1]));
    patchtp = min(patchtp,fsmooth(u-axiswdth/2,[1,0]));
    patchtp = min(patchtp,fsmooth(r-carddmtr/2+2,[1,0]));
    patchtp = patchtp*opctfg;
    patchfg = bsxfun(@plus,bsxfun(@times,patchfg,1-patchtp),lumibg*patchtp);
    % secondary boundary
    u = cos(axis)*x-sin(axis)*y;
    patchtp = ones(outrdmtr,outrdmtr);
    patchtp = min(patchtp,fsmooth(u+axiswdth/2,[0,1]));
    patchtp = min(patchtp,fsmooth(u-axiswdth/2,[1,0]));
    patchtp = min(patchtp,fsmooth(r-carddmtr/2+2,[1,0]));
    patchtp = patchtp*opctfg;
    patchfg = bsxfun(@plus,bsxfun(@times,patchfg,1-patchtp),lumibg*patchtp);
end

% draw card angle
u = sin(angl)*x+cos(angl)*y;
patchtp = ones(outrdmtr,outrdmtr);
patchtp = min(patchtp,fsmooth(u+anglwdth/2,[0,1]));
patchtp = min(patchtp,fsmooth(u-anglwdth/2,[1,0]));
patchtp = min(patchtp,fsmooth(r-carddmtr/2+2,[1,0]));
patchtp = patchtp*opctangl;
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