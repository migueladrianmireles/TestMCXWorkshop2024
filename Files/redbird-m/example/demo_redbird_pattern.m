if(~exist('rbrun','file'))
    addpath(fullfile(pwd, '../matlab'));
end

%% Generate Source/Detector Patterns

% srcpattern = diag(ones(1,16))+diag(ones(1,15),-1);
% srcpattern(1,end)=1;
% srcpattern=permute(repmat(srcpattern,[1,1,16]),[2 3 1]);
% srcpattern=cat(3,srcpattern,permute(srcpattern,[2 1 3]));
% detpattern=srcpattern;


%% Gausian patterns

x = -10:0.2:10;
y = normpdf(x,0,1.5);
detpattern = repmat(y,41,1)';
srcpattern = detpattern;

%%
clear cfg cfg0 recon

% Bounding box w/ 2 inclusions
[cfg0.node,~,cfg0.elem]=meshabox([0 0 0], [120 60 40], 4);
cfg0.seg = ones(size(cfg0.node,1),1);

cfg0.srctype = 'pattern';   % Source setup
cfg0.srcpos = [10 10 0];
cfg0.srcparam1 = [100 0 0 0];
cfg0.srcparam2 = [0 40 0 0];
cfg0.srcdir = [0 0 1];
cfg0.srcpattern = srcpattern; % # src of patterns based upon  size(srcpattern,3)

cfg0.dettype = 'pattern';   % Detector setup
cfg0.detpos = [10 10 40];
cfg0.detparam1 = [100 0 0 0];
cfg0.detparam2 = [0 40 0 0];
cfg0.detdir = [0 0 -1];
cfg0.detpattern = srcpattern;   % Same patterns used as sources
cfg0.srcweight = 1;

cfg0.prop=[
    0 0 1 1
    0.008 1 0 1.37];

%cfg0.omega=2*pi*70e6;
cfg0.omega=0;   % CW (no frequency modulation)

cfg=cfg0;

cfg = rbmeshprep(cfg);

%% run forward simulation to get simulated data

[detphi,phi] = rbrun(cfg);

%%


if(exist('mcxlab','file'))
        xcfg.nphoton=1e9;
        xcfg.vol=uint8(ones(120,60,40));
        xcfg.srcdir=[0 0 1 0];
        xcfg.gpuid=1;
        xcfg.autopilot=1;
        xcfg.prop=cfg.prop;
        xcfg.tstart=0;
        xcfg.tend=5e-9;
        xcfg.tstep=5e-9;
        xcfg.seed=99999;
        xcfg.issrcfrom0=0;
        
        xcfg.isnormalized=1;

        % a uniform planar source outside the volume
        xcfg.srctype='pattern';
        xcfg.srcpos=[10 10 1];
        xcfg.srcparam1=[100 0 0 101];
        xcfg.srcparam2=[0 40 0 41];
        xcfg.srcpattern = squeeze(srcpattern)./max(max(srcpattern));

        flux=mcxlab(xcfg);
        fcw=flux.data*xcfg.tstep;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,

slice = 39;

subplot(211);
imagesc(rot90(log10(abs(squeeze(cfg.srcweight*fcw(:,:,slice))))))
axis equal; colorbar
title('MCX solution');

cl=get(subplot(211),'clim');
subplot(212);
plotmesh([cfg.node full(log10(phi(:,1)))],cfg.elem,sprintf('z=%d',slice+0.5))
view(2);
shading interp;
set(gca,'clim',cl);
colorbar;
title('Redbird solution');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clines = 0.5:-0.5:-10;
[xi,yi] = meshgrid(0.5:119.5,0.5:59.5);
[cutpos,cutvalue,facedata] = qmeshcut(cfg.elem,cfg.node,phi(:,1),sprintf('z=%d',slice+0.5));
vphi = griddata(cutpos(:,1),cutpos(:,2),cutvalue,xi+0.5,yi);

figure,[c,h] = contour(xi,yi,log10(vphi),clines,'r-','LineWidth',2);

cwf = squeeze(fcw(:,:,slice))';
hold on,contour(xi,yi,log10(cfg.srcweight*cwf),clines,'b-','LineWidth',2);
