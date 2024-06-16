function [recon, resid, cfg, updates, Jmua, detphi0iter, phi]=rbrunrecon(maxiter,cfg,recon,detphi0,sd,varargin)

% [newrecon, resid, newcfg]=rbrunrecon(maxiter,cfg,recon,detphi0,sd)
%   or
% [newrecon, resid, newcfg, updates, Jmua, detphi, phi]=rbrunrecon(maxiter,cfg,recon,detphi0,sd,'param1',value1,'param2',value2,...)
%
% Perform a single iteration of a Gauss-Newton reconstruction
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     maxiter: number of iterations
%     cfg: simulation settings stored as a redbird data structure
%     recon: reconstruction data structure, recon may have
%         node: reconstruction mesh node list
%         elem: reconstruction mesh elem list
%         bulk: a struct storing the initial guesses of the param
%              (wavelength-independent optical properties) and prop
%              (wavelength-dependent optical properties), accepted
%              subfields include
%
%              mua/musp/dcoeff/n/g: used to initialize recon/cfg.prop
%              hbo/hbr/scatamp/scatpow: used to initialize recon/cfg.param
%         param: wavelength-independent parameter on the recon mesh
%         prop: wavelength-dependent optical properties on the recon mesh
%         lambda: Tikhonov regularization parameter
%         mapid: the list of element indices of the reconstruction mesh where each forward mesh
%           node is enclosed
%         mapweight: the barycentric coordinates of the forward mesh nodes inside the
%           reconstruction mesh elements
%     detphi0: measurement data vector or matrix
%     sd (optional): source detector mapping table, if not provided, call
%         rbsdmap(cfg) to compute
%     param/value: acceptable optional parameters include
%         'lambda': Tikhonov regularization parameter (0.05), overwrite recon.lambda
%         'report': 1 (default) to print residual and runtimes; 0: silent
%         'tol': convergence tolerance, if relative residual is less than
%                this value, stop, default is 0, which runs maxiter
%                iterations
%         'reform': 'real': transform A*x=b so that A/x/b are all real
%                   'complex': do not transform A*x=b
%                   'logphase': transform Ax=b to [Alogamp,Aphase]*x=[log10(b),angle(b)]
%         'mex': 0 (default) use matlab native code rbjac to build Jacobian
%                1: use mex-file rbfemmatrix to rapidly compute Jacobian
%                  on forward (dense) mesh then interpolate to coarse mesh
%                2: call mex rbfemmatrix to build Jacobian directly on the
%                   recon mesh (coarse).
%                setting mex to 2 gives the fastest speed (2x faster than 0)
%         'prior': apply structure-prior-guided reconstruction,
%                supported methods include
%
%                'laplace': this is also known as the "soft-prior", where
%                    the L matrix used in (J'J+lambda*L'L)dx=dy is a
%                    Laplace smoothing matrix where l(i,j)=1 if i=j or
%                    -1/N_seg if i~=j, where N_seg is the total number of
%                    nodes/elems that are within each label or region;
%                    recon.seg must be a vector of integer labels
%                'comp': use compositional-priors, recon.seg must be a
%                    N-by-Nc matrix where N is the number of nodes, Nc is
%                    the number of tissue compositions, each element in the
%                    matrix must be a number between 0-1, denoting the
%                    volume fraction of each composition; the row-sum must
%                    be 1 for each node.
%
% output:
%     recon: the updated recon structure, containing recon mesh and
%          reconstructed values in recon.prop or recon.param
%     resid: the residual betweet the model and the measurement data for
%          each iteration
%     cfg: the updated cfg structure, containing forward mesh and
%          reconstructed values in cfg.prop or cfg.param
%     updates: a struct array, where the i-th element stores the update
%          vectors for each unknown block
%     Jmua: Jacobian in a struct form, each element is the Jacobian of an
%          unknown block
%     detphi: the final model prediction that best fits the data detphi0
%     phi: the final forward solutions resulting from the estimation
%     
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

if(maxiter==0 && nargin<3)
    % return detphi as recon and phi as resid
    [recon,resid]=rbrunforward(cfg);
    return;
end

resid=zeros(1,maxiter);
updates=repmat(struct,1,maxiter);

opt=varargin2struct(varargin{:});

lambda=0.05;
if(isfield(recon,'lambda'))
    lambda=recon.lambda;
end
lambda=jsonopt('lambda',lambda,opt);

doreport=jsonopt('report',1,opt);
convergetol=jsonopt('tol',0,opt);
reform=jsonopt('reform','real',opt);
ismexjac=jsonopt('mex',0,opt);
prior=jsonopt('prior','',opt);
rfcw = jsonopt('rfcw',1,opt);
solverflag=jsonopt('solverflag',{},opt);
blockscale = jsonopt('blockscale',1,opt);
musscale = jsonopt('musscale',0.5,opt);
mode = jsonopt('mode','image',opt);
tmesh = jsonopt('templatemesh',0,opt);
debugplot = jsonopt('debugplot',0,opt);
isreduced=0;

% create or accept regularization matrix

Aregu=struct;

if(isfield(opt,'lmat'))
    Aregu.lmat=opt.lmat;
elseif(isfield(opt,'ltl'))
    Aregu.ltl=opt.ltl;
elseif(isfield(opt,'lir'))
    Aregu.lir=opt.lir;
end

if(~isempty(prior) && isfield(recon,'seg') && ~isfield(Aregu,'lmat'))
    Aregu.lmat=rbprior(recon.seg,prior,opt);
end
if(nargin<5)
   sd = rbsdmap(cfg);
end

if isfield(opt,'rfcw')
    rfcw = rfcw;
else
    if isa(sd,'containers.Map')
        waves = sd.keys;
        sdwv = sd(waves{1});
        if (size(sdwv,2) > 3)
            rfcw = unique(sdwv(:,4));
            rfcw = rfcw(rfcw > 0)';
        else
            rfcw = 1;
        end
    end
end

if (length(rfcw) < 2 && (isstruct(detphi0) && length(detphi0) > 1))
    detphi0 = detphi0(rfcw).detphi;
end

% start iterative Gauss-Newton based reconstruction
for iter=1:maxiter
    tic
    
    if isfield(recon,'param')
        recon.param;
    elseif (size(recon.prop,1) < 6)
        recon.prop;
    end
    
    % update forward mesh prop/param using recon mesh prop/param if given
    % for rbsyncprop to work, one must provide initial values of cfg.prop
    % (or cfg.param) if recon.prop (or recon.param) is specified
    if((isfield(recon,'node') && isfield(recon,'elem')) || isfield(recon,'prop') || isfield(recon,'param'))
        [cfg,recon]=rbsyncprop(cfg,recon);
    end

    if(isfield(cfg,'param') && isstruct(cfg.param) && all(structfun(@isempty,cfg.param)==0))
        if(isfield(cfg,'prop') && isa(cfg.prop,'containers.Map') && ~isempty(keys(cfg.prop)))
            cfg.prop=rbupdateprop(cfg);
        end
    end
    
    if (isfield(recon,'param') && length(recon.param.hbo) == length(recon.node) && debugplot==1)
        figure(15),plotmesh([recon.node recon.param.hbo+recon.param.hbr],recon.elem, char(strcat('z=',num2str(round(max(recon.node(:,3))./2)))));colorbar;shading interp
    end
    
    % Run forward on forward mesh
    [detphi,phi] = rbrunforward(cfg,'solverflag',solverflag,'sd',sd,'rfcw',rfcw);
    
%     % Extract simulated fluence images (z = 1/musp) to compute SP and image residual
%     temp = OMCIEP3_CWGetFluenceImage(cfg,phi,round([size(cfg.phi0img,2) size(cfg.phi0img,1)]));
%     fwd_img(:,:,1:32,1) = temp(:,:, 1:32);
%     fwd_img(:,:,1:32,2) = temp(:,:,33:64);
%     clear temp
% 
%     tempdetpatternTarg(:,:,:,1) = rot90(cfg.widedetid('690').srcpattern{1});
%     tempdetpatternTarg(:,:,:,2) = rot90(cfg.widedetid('830').srcpattern{1});
%     %tempw = nansum(nansum(abs(tempdetpatternTarg),1),2)./max(nansum(nansum(abs(tempdetpatternTarg),1),2));
%     tempw = 1;%
%     tempSPfwd = CWSinglePixelMW(fwd_img,(tempdetpatternTarg./nansum(nansum(abs(tempdetpatternTarg),1),2)).*tempw);
%     
%     if isstruct(detphi) % RFCW bulk fitting
%         detphi(2).detphi('690') = reshape(tempSPfwd(:,1),[size(cfg.phi0img,3) size(cfg.phi0img,3)]);
%         detphi(2).detphi('830') = reshape(tempSPfwd(:,2),[size(cfg.phi0img,3) size(cfg.phi0img,3)]);
%         clear tempdetpatternTarg tempSPfwd
%     else % CW reconstruction
%         detphi('690') = reshape(tempSPfwd(:,1),[size(cfg.phi0img,3) size(cfg.phi0img,3)]);
%         detphi('830') = reshape(tempSPfwd(:,2),[size(cfg.phi0img,3) size(cfg.phi0img,3)]);
%         clear tempdetpatternTarg tempSPfwd
%     end
% 
%     % Edit 02/13/2024 to include image resid
%     if isfield(cfg,'phi0img') % Compute dref error
%         
%         cfg.residImg(1,iter) = nansum(abs(cfg.phi0img(:) - fwd_img(:)));
%         
%         fprintf(sprintf('Image residual = %d\n',cfg.residImg(1,iter)));
%         if strcmp(mode,'image')
%             cfg.fwd_img(:,:,:,:,iter)    = fwd_img;
%             cfg.residImg2D(:,:,:,:,iter) = cfg.phi0img - fwd_img;
%         end
%         
%         cbar  = [0 20e-7];
%         cbar2 = [0 90e-7];
%         figure,
%         subplot(6,2,1),imagesc(    fwd_img(:,:,10,1)),colorbar,axis equal, axis tight,caxis(cbar)
%         title(num2str(cfg.residImg(1,iter)))
%         subplot(6,2,2),imagesc(cfg.phi0img(:,:,10,1)),colorbar,axis equal, axis tight,caxis(cbar)
%         subplot(6,2,3),imagesc(    fwd_img(:,:,19,1)),colorbar,axis equal, axis tight,caxis(cbar)
%         subplot(6,2,4),imagesc(cfg.phi0img(:,:,19,1)),colorbar,axis equal, axis tight,caxis(cbar)
%         subplot(6,2,5),imagesc(    fwd_img(:,:,27,1)),colorbar,axis equal, axis tight,caxis(cbar)
%         subplot(6,2,6),imagesc(cfg.phi0img(:,:,27,1)),colorbar,axis equal, axis tight,caxis(cbar)
%         subplot(6,2,7),imagesc(    fwd_img(:,:,10,1)),colorbar,axis equal, axis tight,caxis(cbar2)
%         subplot(6,2,8),imagesc(cfg.phi0img(:,:,10,2)),colorbar,axis equal, axis tight,caxis(cbar2)
%         subplot(6,2,9),imagesc(    fwd_img(:,:,19,2)),colorbar,axis equal, axis tight,caxis(cbar2)
%         subplot(6,2,10),imagesc(cfg.phi0img(:,:,19,2)),colorbar,axis equal, axis tight,caxis(cbar2)
%         subplot(6,2,11),imagesc(    fwd_img(:,:,27,2)),colorbar,axis equal, axis tight,caxis(cbar2)
%         subplot(6,2,12),imagesc(cfg.phi0img(:,:,27,2)),colorbar,axis equal, axis tight,caxis(cbar2)        
%         drawnow
%         clear fwd_img
% 
%     end
   
    if isfield(cfg,'sdscale')
        for ii = 1:length(phi)
            for wv = phi(ii).phi.keys
                if ii == 1
                    phi(ii).phi(wv{1}) = phi(ii).phi(wv{1}).*complex(cfg.sdscale(:,1),cfg.sdscale(:,2))';
                elseif ii == 2
                    phi(ii).phi(wv{1}) = phi(ii).phi(wv{1}).*cfg.sdscale(:,1)';
                end
            end
        end
    end

    % build Jacobians on forward mesh
    if isa(cfg.omega,'containers.Map')
        omegas = cell2mat(cfg.omega.values);
    elseif (isfield(cfg,'omega') )
        omegas = cfg.omega;
      else
        omegas = 0;
    end
    if(isfield(cfg,'omega') && (any(omegas>0)) && ismember(1,rfcw)) % if RF data
%     if(isfield(cfg,'omega') && (cfg.omega>0 || any(omegas>0)) && ismember(1,rfcw)) % if RF data
        % currently, only support node-based values; rbjac supports
        % multiple wavelengths, in such case, it returns a containers.Map
        if((isfield(cfg,'seg') && length(cfg.seg)==size(cfg.elem,1)) || size(cfg.prop,1)==size(cfg.elem,1))
            % element based properties
            [Jmua,~,Jd]=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol, 1,'rfcw', rfcw);
        else
            % node based properties
            if(ismexjac) % use mex to build
                if(ismexjac>=2 && isfield(recon,'node'))
                    [Jmua, Jd]=rbjacmex(cfg, sd, phi, cfg.deldotdel, 1, ...
                        recon.mapid, recon.mapweight, size(recon.node,1), recon.elem);
                else
                    [Jmua, Jd]=rbjacmex(cfg, sd, phi);
                end
            else
                [Jmua,~,Jd]=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol,'rfcw',rfcw);
%                 [Jmua,Jd] = rbjtestnode(sd,phi,cfg.deldotdel,cfg.elem,cfg.evol,cfg);
            end
        end
    else % CW only
        if((isfield(cfg,'seg') && length(cfg.seg)==size(cfg.elem,1)) || size(cfg.prop,1)==size(cfg.elem,1))
            Jmua=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol, 1);
        else
            if(ismexjac)
                if(ismexjac>=2 && isfield(recon,'node'))
                    Jmua=rbjacmex(cfg, sd, phi, cfg.deldotdel, 1, ...
                        recon.mapid, recon.mapweight, size(recon.node,1), recon.elem);
                else
                    Jmua=rbjacmex(cfg, sd, phi);
                end
            else
                Jmua=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol,'rfcw',rfcw);
            end
        end
    end
    % Jmua/Jd are either containers.Map(wavelength) or single matrix
    
    % build Jacobians for chromophores in the form of a struct
    % TODO: need to handle Jmua is a map but cfg.param is not defined
    if(isa(cfg.prop,'containers.Map') && isfield(cfg,'param') && isa(cfg.param,'struct'))
        if(exist('Jd','var'))
            [Jmua,detphi0iter,detphi]=rbmultispectral(sd, cfg, Jmua, detphi0, detphi, cfg.param, rfcw, Jd, cfg.prop);
        else
            [Jmua,detphi0iter,detphi]=rbmultispectral(sd, cfg, Jmua, detphi0, detphi, cfg.param, rfcw);
        end
    else  % recon mua/d per wavelengths
        if(exist('Jd','var'))
            Jmua=struct('mua',Jmua,'dcoeff',Jd);
            sdkeep = find(sd(:,3) == 1);
            detphi0iter = detphi0(:);
            detphi = detphi(:);
            detphi0iter = detphi0iter(sdkeep);
            detphi = detphi(sdkeep);
        else
            Jmua=struct('mua',Jmua);
            detphi0iter = detphi0;
            detphi = detphi;
        end
    end

    % here, Jmua is a struct of unknown species; wavelengths
    % are vertically/row-wise concatenate; detphi and detphi0 are the
    % concatenated model and measurement RHS vectors

    % mapping jacobians from forward mesh to reconstruction mesh
    if(ismexjac<2 && isfield(recon,'elem') && isfield(recon,'node') && isfield(recon,'mapid') && isfield(recon,'mapweight')) % dual-mesh reconstruction
        if isstruct(Jmua)
            for ii = 1:length(Jmua)
                Jmua(ii)=structfun(@(x) transpose(meshremap(x.',recon.mapid, recon.mapweight,recon.elem,size(recon.node,1))), Jmua(ii),'UniformOutput',false); 
            end
        else
            Jmua=structfun(@(x) transpose(meshremap(x.',recon.mapid, recon.mapweight,recon.elem,size(recon.node,1))), Jmua,'UniformOutput',false); 
        end
    end

    if(isfield(recon,'seg') && isvector(recon.seg)) % reconstruction of segmented domains
        if (isstruct(Jmua) && (length(Jmua) > 1))
            for ii = 1:length(Jmua)
                Jmua(ii) = structfun(@(x) rbmasksum(x,recon.seg(:)'), Jmua(ii),'UniformOutput',false);
            end
        else
            Jmua=structfun(@(x) rbmasksum(x,recon.seg(:)'), Jmua,'UniformOutput',false);
        end
        isreduced=1;
    elseif(isfield(cfg,'seg') && isvector(cfg.seg)) % single-mesh bulk/seg recon
        if (isstruct(Jmua) && (length(Jmua) > 1))
            for ii = 1:length(Jmua)
                Jmua(ii)=structfun(@(x) rbmasksum(x,cfg.seg(:)'), Jmua(ii),'UniformOutput',false);
            end
        else
            Jmua=structfun(@(x) rbmasksum(x,cfg.seg(:)'), Jmua,'UniformOutput',false);
        end
        isreduced=1;
    end

    % blocks contains unknown names and Jacob size, should be Nsd*Nw rows
    % and Nn columns, Nsd is src/det pairs, Nw is number of wavelengths,
    % and Nn is the recon mesh, if present, or forward mesh node size.
    if (isstruct(Jmua) && (length(Jmua) > 1))
        for ii = 1:length(Jmua)
            blocks(ii) = structfun(@size, Jmua(ii), 'UniformOutput', false);
        end
    else
        blocks=structfun(@size, Jmua, 'UniformOutput', false);
    end
    
    % flatten Jmua into a horizontally/column contatenated matrix
    if(strcmp(reform,'complex')==0)
        [Jflat,misfit]=rbmatreform(rbmatflat(Jmua), detphi0iter(:), detphi(:), reform);
    else
        Jflat=rbmatflat(Jmua);
        misfit=detphi0iter(:)-detphi(:);
    end

    % store the residual
    resid(1,iter)         = sum(abs(misfit));
    updates(iter).detphi  = detphi;
        
    %cfg.absresid(1,iter)  = sum(abs(misfit)); 
    %cfg.CWSPmodel(:,iter) = detphi;
    %cfg.CWSPmeasured(:,1) = detphi0iter;
    
    % if(isreduced==0 && ~isempty(prior))
    if(isreduced == 0 && ~isempty(fieldnames(Aregu)) && iter == 1)
        if(size(Jflat,1)>=size(Jflat,2) && ~isfield(Aregu,'ltl'))
            Aregu.ltl=Aregu.lmat'*Aregu.lmat;
        end
        if(size(Jflat,1)<size(Jflat,2) && ~isfield(Aregu,'lir'))
            Lr=qr(Aregu.lmat);
            Aregu.lir=inv(triu(Lr));
        end
    end
    
    % Local recon test EXu 05262023
    if isfield(recon,'Jmask')
        if (size(Jflat,2) > size(recon.node,1))
            blockfields = fieldnames(blocks);
            for yy = 1:length(blockfields)
                Jflat(:,(yy-1).*size(recon.node,1) + setdiff(1:size(recon.node,1),recon.Jmask)) = 0;
            end
        else
            Jflat(:,setdiff(1:size(recon.node,1),recon.Jmask)) = 0;
        end
    end
    
    if ((blockscale == 1) && isstruct(Jmua))% && strcmp(mode,'image'))
        blockfields = fieldnames(blocks);
        cn = blocks(1).(blockfields{1})(2);
        scalefact = [];
        for zz = 1:length(fieldnames(Jmua))
            if (zz<=length(intersect(fieldnames(Jmua),{'hbo','hbr','mua'})))
                scalefact(zz) = 1/sqrt(sum(sum(Jflat(:,(zz-1)*cn+1:zz*cn).^2)));
            elseif (exist('Jd','var') &&  zz>length(intersect(fieldnames(Jmua),{'scatamp','scatpow','dcoeff'})))
                scalefact(zz) = 1/sqrt(sum(sum(Jflat(:,(zz-1)*cn+1:zz*cn).^2))).*musscale;
            end
            Jflat(:,(zz-1)*cn+1:zz*cn) = Jflat(:,(zz-1)*cn+1:zz*cn).*scalefact(zz);
        end
    end
             
    %% Adding L-curve method for finding regularization parameter (lambda)
       
%     if iter == 1% && (cfg.isbulkLcurve == 1 || cfg.isLcurve == 1)
%         [U,s,~] = csvd(Jflat); % Decomposition of Adjoint Jacobian @830nm
%         figure,
%         cfg.lambda = l_curve(U,s,misfit,'Tikh');
%         drawnow
%     end
        
    %%
    
    % solver the inversion (J*delta_x=delta_y) using regularized
    % Gauss-Newton normal equation
    if isfield(recon,'Jmask')
        dmu_recon = rbreginvover(Jflat,misfit,lambda,Aregu.ltl,blocks);
    else
        dmu_recon = rbreginv(Jflat, misfit, lambda, Aregu, blocks, solverflag{:});  % solve the update on the recon mesh
    end
    
    if (blockscale == 1 && exist('scalefact','var'))
        for zz = 1:length(scalefact)
            dmu_recon((zz-1)*cn+1:zz*cn) = dmu_recon((zz-1)*cn+1:zz*cn).*scalefact(zz);
        end
    end
   
    % obtain linear index of each output species
    len=cumsum([1; structfun(@(x) x(2), blocks(1))]);
    output=fieldnames(blocks);
    for i=1:length(output)
        dx=dmu_recon(len(i):len(i+1)-1);
        updates(iter).(output{i})=dx;
        switch output{i}
            case {'mua','dcoeff'}
                propidx=strcmp(output{i},'dcoeff')+1;
                if(length(dx)==size(recon.prop,1)-1) % label based prop
                    if(strcmp(output{i},'dcoeff')) % converting from d to musp
                        dcoeff=1./(3*recon.prop(2:end,propidx));
                        dcoeff=dcoeff+dx;
                        recon.prop(2:end,propidx)=1./(3*dcoeff);
                    else
                        recon.prop(2:end,propidx)=recon.prop(2:end, propidx)+dx;
                    end
                    cfg.prop=recon.prop;
                else % nodal or element based prop
                    if(strcmp(output{i},'dcoeff')) % converting from d to musp
                        dcoeff=1./(3*recon.prop(:, propidx));
                        dcoeff=dcoeff+dx;
                        recon.prop(:,propidx)=1./(3*dcoeff);
                    else
                        recon.prop(:,propidx)=recon.prop(:, propidx)+dx;
                    end
                    if(isfield(recon,'node'))
                        cfg.prop(:,propidx)=...
                            meshinterp(recon.prop(:,propidx),recon.mapid, recon.mapweight,recon.elem,cfg.prop(:,propidx)); % interpolate the update to the forward mesh
                    else
                        cfg.prop=recon.prop;
                    end
                end
            case {'hbo','hbr','water','lipid','scatamp','scatpow'}
                if(isfield(recon,'node')) % update recon mesh prop
                    recon.param.(output{i})=recon.param.(output{i})+dx;
                else
                    cfg.param.(output{i})=cfg.param.(output{i})+dx;
                end
            otherwise
                error('unknown type %s is not supported',output{i});
        end
    end
    if(doreport)
        fprintf(1,'iter [%4d]: residual=%e, relres=%e lambda=%e (time=%f s)\n',iter, resid(iter), resid(iter)/resid(1), lambda, toc);
    end
    if(iter>1 && abs((resid(iter)-resid(iter-1))/resid(1))<convergetol)
        resid=resid(1:iter);
        updates=updates(1:iter);
        break;
    end
    
    
%% EXTRA STUFF ADDED BY MM TO VISUALIZE FITTING PROCESS    

% if isstruct(detphi0iter) % RFCW bulk fitting
%     measDataRF  = detphi0iter(1).detphi;
%     modelDataRF = detphi(1).detphi;
%     if size(detphi0iter,2) > 1
%         measDataCW  = detphi0iter(2).detphi;
%         modelDataCW = detphi(2).detphi;
%     end
%     figure,
%     subplot(2,2,1),hold on,plot(abs(measDataRF) ,'DisplayName','Measured'),title('RF Amplitude')
%     subplot(2,2,1),hold on,plot(abs(modelDataRF),'DisplayName','Model')
%     legend
%     subplot(2,2,2),hold on,plot(angle(measDataRF) ,'DisplayName','Measured'),title('RF Phase')
%     subplot(2,2,2),hold on,plot(angle(modelDataRF),'DisplayName','Model')
%     legend
%     if size(detphi0iter,2) > 1
%         subplot(2,2,3:4),hold on,plot(measDataCW ,'DisplayName','Measured'),title('SP-CW')
%         subplot(2,2,3:4),hold on,plot(modelDataCW,'DisplayName','Model')
%     end
%     legend
%     drawnow
% else % CW only reconstruction
%     
%     measDataRF  = detphi0iter(:);
%     modelDataRF = detphi(:);
%     figure,
%     subplot(3,3,1:3),hold on,plot((measDataRF) ,'DisplayName','Measured'),title(char(strcat('Iter:',num2str(iter),' SP-CW lambda:',num2str(lambda),' relresid:',num2str(resid(1,iter)./resid(1,1)))))
%     subplot(3,3,1:3),hold on,plot((modelDataRF),'DisplayName','Model')
%     legend
%     miguel = cfg.prop('690');
%     adrian = cfg.prop('830');
%     subplot(3,3,4),plotmesh([cfg.node miguel(:,1)],cfg.elem, char(strcat('z=',num2str(round(max(cfg.node(:,3))./2)))));title('mua690');colorbar;shading interp,view(2)
%     subplot(3,3,5),plotmesh([cfg.node adrian(:,1)],cfg.elem, char(strcat('z=',num2str(round(max(cfg.node(:,3))./2)))));title('mua830');colorbar;shading interp,view(2)
%     subplot(3,3,6),plotmesh([recon.node recon.param.hbo+recon.param.hbr],recon.elem, char(strcat('z=',num2str(round(max(recon.node(:,3))./2)))));title('HbT');colorbar;shading interp,view(2)
%     
%     clear xi yi cutpos cutvalue vphi
%     yaxis_ref = 190;
%     [xi,yi] = meshgrid(min(cfg.node(:,1),[],1):1:max(cfg.node(:,1),[],1),min(cfg.node(:,2),[],1):1:max(cfg.node(:,2),[],1));
%     [cutpos,cutvalue,~] = qmeshcut(cfg.elem,cfg.node,miguel(:,1),char(strcat('z =',num2str(max(cfg.node(:,3))/2))));
%     vphi(:,:,1) = griddata(cutpos(:,1),cutpos(:,2),cutvalue,xi,yi);
%     clear cutpos cutvalue
%     [cutpos,cutvalue,~] = qmeshcut(cfg.elem,cfg.node,adrian(:,1),char(strcat('z =',num2str(max(cfg.node(:,3))/2))));
%     vphi(:,:,2) = griddata(cutpos(:,1),cutpos(:,2),cutvalue,xi,yi);
%     
%     [~,b] = min(abs(yi(:,1) - yaxis_ref),[],1);
%     subplot(3,3,7),
%     hold on,imagesc(xi(:,1),yi(1,:),vphi(:,:,1)),axis equal,axis tight
%     hold on,plot(xi(1,:),yi(b,:),'r-.','linewidth',2)
%     subplot(3,3,8),
%     hold on,imagesc(xi(:,1),yi(1,:),vphi(:,:,2)),axis equal,axis tight
%     hold on,plot(xi(1,:),yi(b,:),'b-.','linewidth',2)
%     subplot(3,3,9),
%     hold on,plot(xi(1,:),vphi(b,:,1),'r','linewidth',4)
%     hold on,plot(xi(1,:),vphi(b,:,2),'b','linewidth',4)
%     title(['Cross section @y =' num2str(yaxis_ref), ' mm'])
%     axis tight
%     drawnow
%     
%     clear miguel adrian
%         
% end
     
end

if(isfield(recon,'node') && isfield(recon,'elem'))
    [cfg,recon]=rbsyncprop(cfg,recon);
    cfg.prop=rbupdateprop(cfg);
end
