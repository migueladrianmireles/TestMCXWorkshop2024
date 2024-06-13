%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - A Diffusion Solver for Diffuse Optical Tomography, 
%      Copyright Qianqina Fang, 2018
%
% In this example, we show the most basic usage of Redbird.
%
% This file is part of Redbird URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~exist('rbrun','file'))
    addpath(fullfile(pwd, '../matlab'));
end

%%

dim = [60 60 60];
resgrid = 1; 

[cfg.node, cfg.elem] = meshgrid6(0:resgrid:dim(1),0:resgrid:dim(2),0:resgrid:dim(3)); % For homogeneous domain
cfg.elem(:,1:4)      = meshreorient(cfg.node(:,1:3),cfg.elem(:,1:4));
cfg.face             = volface(cfg.elem);

nn=size(cfg.node,1);
cfg.seg=ones(size(cfg.elem,1),1);

cfg.srcpos = dim.*[1 1 0]./2;
cfg.srcdir=[0 0 1];

cfg.detpos = [30 30 60];
cfg.detdir = [0 0 -1];

cfg.prop=[0 0 1 1;0.02 10 0.9 1.37];

cfg.omega=0;

cfg=rbmeshprep(cfg);

%%

[Amat,deldotdel] = rbfemlhs(cfg); % Build LHS
[rhs,loc,bary]   = rbfemrhs(cfg); % Build RHS
phi=rbfemsolve(Amat,rhs);
% phi(phi<0)=0;
detval=rbfemgetdet(phi, cfg, rhs); % or detval=rbfemgetdet(phi, cfg, loc, bary);

%%

figure,plotmesh([cfg.node log(abs(phi(:,1)))],cfg.elem,'y=29'),shading interp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   prepare simulation input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg

[cfg.node,cfg.face, cfg.elem]=meshabox([40 0 0], [160, 120, 60], 10);
cfg.seg=ones(size(cfg.elem,1),1);

[xi,yi]=meshgrid(60:20:140,20:20:100);
cfg.srcpos=[xi(:),yi(:),zeros(numel(yi),1)];
cfg.detpos=[xi(:),yi(:),60*ones(numel(yi),1)];
cfg.srcdir=[0 0 1];
cfg.detdir=[0 0 -1];

cfg.prop=[
    0 0 1 1
    0.008 1 0 1.37
    0.016 1 0 1.37
];


%cfg.omega=2*pi*70e6;
cfg.omega = 0 ;

tic
cfg=rbmeshprep(cfg);
fprintf('preparing mesh ... \t%f seconds\n',toc);

% save config.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Solve the forward problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
[detphi, phi]=rbrun(cfg);
fprintf('forward solution ... \t%f seconds\n',toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Visualization of fluence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rbplotforward(cfg,log10(abs(phi(:,13))),1,'y>60');
hold on
plotmesh([cfg.node log10(abs(phi(:,13)))],cfg.elem,'y=57','linestyle','none')
