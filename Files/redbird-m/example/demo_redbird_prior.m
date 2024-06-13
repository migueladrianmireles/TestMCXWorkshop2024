%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - A Diffusion Solver for Diffuse Optical Tomography, 
%      Copyright Qianqina Fang, 2018
%
% Continuous-Wave (CW) reconstruction of absorption (mua) target
% (streamlined version by calling rbrunrecon)
%
% This file is part of Redbird URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~exist('rbrun','file'))
    addpath(fullfile(pwd, '../matlab'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   prepare simulation input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg cfg0 recon

s0=[70, 50, 20];

[nobbx,fcbbx]=meshabox([40 0 0], [160, 120, 60], 25);
[nosp,fcsp]=meshasphere(s0, 5, 10);
[no,fc]=mergemesh(nobbx, fcbbx, nosp, fcsp);

[cfg0.node, cfg0.elem]=s2m(no,fc(:,1:3),1,100,'tetgen1.5',[41 1 1;s0]);
nn=size(cfg0.node,1);
cfg0.seg=cfg0.elem(:,5);
cfg0.srcdir=[0 0 1];

[xi,yi]=meshgrid(60:20:140,20:20:100);
cfg0.srcpos=[xi(:),yi(:),zeros(numel(yi),1)];
cfg0.detpos=[xi(:),yi(:),60*ones(numel(yi),1)];
cfg0.detdir=[0 0 -1];

cfg0.prop=[
    0 0 1 1
    0.008 1 0 1.37
    0.016 1 0 1.37
];

cfg0.omega=67.5e6*2*pi;

cfg=cfg0;

cfg0=rbmeshprep(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run forward for the heterogeneous domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detphi0=rbrun(cfg0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Reset the domain to a homogeneous medium for recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create forward mesh for reconstruction - if it has a different density
% compared to the forward mesh used for generating data, it may introduce
% numerical error - so, set density to 10 gives the best result, need
% to debug this further

[node,face,elem]=meshabox([40 0 0], [160, 120, 60], 25);
cfg=rbsetmesh(cfg,node,elem,cfg.prop,ones(size(node,1),1));

sd=rbsdmap(cfg);

% create coarse reconstruction mesh
[recon.node,face,recon.elem]=meshabox([40 0 0], [160, 120, 60], 50);
[recon.mapid, recon.mapweight]=tsearchn(recon.node,recon.elem,cfg.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Streamlined reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize reconstruction to homogeneous (label=1)
recon.prop=cfg.prop(ones(size(recon.node,1),1)+1,:);
cfg.prop=cfg.prop(ones(size(cfg.node,1),1)+1,:);
cfg=rmfield(cfg,'seg');

% run stream-lined image reconstruction
[newrecon,resid,newcfg,~,~,phi]=rbrun(cfg,recon,detphi0,sd,'mex',0,'maxiter',5,'lambda',5e-3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Create prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[vv, nn] = gaussiankernel(recon.node,s0,5);
% vv(vv > 1e-4) = 1;
% vv(vv <= 1e-4) = 0;
% prior = [abs(1 - vv) vv];
vv = vv./max(vv);
prior = [abs(1-vv) vv];

L = rbmakeL(cfg0,recon,prior);

%
[newreconL,residL,newcfgL,~,~,phiL]=rbrun(cfg,recon,detphi0,sd,'mex',0,'maxiter',5,'lambda',5e-3,'lmat',L);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plotting results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,
subplot(1,2,1)
plotmesh([newrecon.node,newrecon.prop(:,1)],newrecon.elem,'z=20','facecolor','interp','linestyle','none')
hold on;
plotmesh([newrecon.node,newrecon.prop(:,1)],newrecon.elem,'x=70','facecolor','interp','linestyle','none')
view(3);colorbar;

subplot(1,2,2)
plotmesh([newreconL.node,newreconL.prop(:,1)],newreconL.elem,'z=20','facecolor','interp','linestyle','none')
hold on;
plotmesh([newreconL.node,newreconL.prop(:,1)],newreconL.elem,'x=70','facecolor','interp','linestyle','none')
view(3);colorbar;