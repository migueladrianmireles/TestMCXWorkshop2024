function [detval, phi, Amat, rhs, sflag]=rbrunforward(cfg,varargin)
%
% [detval, phi]=rbrunforward(cfg)
%    or
% [detval, phi, Amat, rhs]=rbrunforward(cfg,'param1',value1,...)
%
% Perform forward simulations at all sources and all wavelengths based on the input structure
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the redbird data structure
%
% output:
%     detval: the values at the detector locations
%     phi: the full volumetric forward solution computed at all wavelengths
%     Amat: the left-hand-side matrices (a containers.Map object) at specified wavelengths 
%     rhs: the right-hand-side vectors for all sources (independent of wavelengths)
%     param/value pairs: (optional) additional parameters
%          'solverflag': a cell array to be used as the optional parameters
%               for rbfemsolve (starting from parameter 'method'), for
%               example  rbrunforward(...,'solverflag',{'pcg',1e-10,200})
%               calls rbfemsolve(A,rhs,'pcg',1e-10,200) to solve forward
%               solutions
% 
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

opt=varargin2struct(varargin{:});
rfcw = jsonopt('rfcw',[1],opt);

if(~isfield(cfg,'deldotdel'))
    cfg.deldotdel=rbdeldotdel(cfg);
end

wavelengths={''};
if(isa(cfg.prop,'containers.Map'))
   wavelengths=cfg.prop.keys;
end
sd = jsonopt('sd',containers.Map(wavelengths,cell(1,length(wavelengths))),opt);
if (isempty(sd(wavelengths{1})))
    sd = rbsdmap(cfg);
    if (~isa(sd,'containers.Map'))
        sd = containers.Map(wavelengths,{sd});
    end
end

Amat=containers.Map();
opt=varargin2struct(varargin{:});
solverflag=jsonopt('solverflag',{},opt);

for md = rfcw
    detval(md).detphi = containers.Map();
    phi(md).phi = containers.Map();
end

sdtmp = cell2mat(sd.values');
srcnum = length(unique(sdtmp(:,1)));

for waveid=wavelengths
	wv=waveid{1};
    
    sdwv = sd(wv);
    for md = rfcw
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%   Build RHS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [rhs,loc,bary]=rbfemrhs(cfg,sd,wv,md);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%   Build LHS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Amat(wv)=rbfemlhs(cfg,cfg.deldotdel,wv,md); % use native matlab code, 1 sec for 50k nodes

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%   Solve for solutions at all nodes: Amat*res=rhs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %solverflag={'pcg',1e-12,200}; % if iterative pcg method is used
        [phi(md).phi(wv),sflag]=rbfemsolve(Amat(wv),rhs,solverflag{:});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%   Extract detector readings from the solutions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tempdetval = rbfemgetdet(phi(md).phi(wv), cfg, rhs); % or detval=rbfemgetdet(phi(wv), cfg, rhs);
        if size(sdwv,2) < 4
            sdwv(:,4) = ones(size(sdwv,1),1).*md;
        end
        sdmd = sdwv(sdwv(:,4) == md | sdwv(:,4) == 3,:);
        detval(md).detphi(wv) = tempdetval(unique(sdmd(:,2))-srcnum,unique(sdmd(:,1)));
    end
end

% if only a single wavelength is required, return regular arrays instead of a map
if(length(wavelengths)==1)
    Amat=Amat(wavelengths{1});
    for md = rfcw
        phi(md).phi = phi(md).phi(wavelengths{1});
        detval(md).detphi = detval(md).detphi(wavelengths{1});
    end
end

if (length(rfcw) == 1)
    phi = phi(rfcw).phi;
    detval = detval(rfcw).detphi;
end