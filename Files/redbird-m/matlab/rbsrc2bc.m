function cfg=rbsrc2bc(cfg,isdet)
%
% srcbc=rbsrc2bc(cfg)
%
% Converting wide-field source forms into a boundary condition by defining
% in-ward flux on the mesh surface
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the simulation data structure, with srctype, srcpos, srcdir, 
%          srcparam1, srcparam2, srcpattern fields
%     isdet: default is 0; if set to 1, rbsrc2bc process widefield
%          detectors, the relevant fields are dettype, detpos, detdir,
%          detparam1, detparam2, detpattern
%
% output:
%     srcbc: an array of Ns x Nt, where Nt is size(cfg.face,1) and Ns is
%            the number of sources (if isdet=1, the detector counts)
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

srcbc=[];

if(nargin<2)
    isdet=0;
end

if (~isdet)
    if((~isfield(cfg,'srctype') && ~isfield(cfg,'widesrcid')) || (isfield(cfg,'srctype') && (strcmp(cfg.srctype,'pencil') || strcmp(cfg.srctype,'isotropic'))))
        return;
    end
    srcdir = cfg.srcdir;
    sources = cfg.srcpos;
    if isfield(cfg,'widesrcid')
        widesrcid = cfg.widesrcid;
        if ~isa(widesrcid,'containers.Map')
            widesrcid = containers.Map({''},{widesrcid});
        end
    else
        tempwf.srctype = {cfg.srctype};
        if isfield(cfg,'srcid')
            tempwf.srcid = {cfg.srcid};
        else
            tempwf.srcid = {1};
        end
        try
            tempwf.srcparam1 = {cfg.srcparam1};tempwf.srcparam2 = {cfg.srcparam2};
        catch
            error('Widefield sources must have srcparam1 and srcparam2');
        end
        if isfield(cfg,'srcpattern')
            tempwf.srcpattern = {cfg.srcpattern};
        end
        if isfield(cfg,'srcweight')
            tempwf.srcweight = {cfg.srcweight};
        end
        widesrcid = containers.Map({''},{tempwf});clear tempwf
    end
else
    if((~isfield(cfg,'dettype') && ~isfield(cfg,'widedetid')) || (isfield(cfg,'dettype') && (strcmp(cfg.dettype,'pencil') || strcmp(cfg.dettype,'isotropic'))))
        return;
    end
    srcdir = cfg.detdir;
    sources = cfg.detpos;
    if isfield(cfg,'widedetid')
        widesrcid = cfg.widedetid;
        if ~isa(widesrcid,'containers.Map')
            widesrcid = containers.Map({''},{widesrcid});
        end
    else
        tempwf.srctype = {cfg.dettype};
        if isfield(cfg,'detid')
            tempwf.srcid = {cfg.detid};
        else
            tempwf.srcid = {1};
        end
        try
            tempwf.srcparam1 = {cfg.detparam1};tempwf.srcparam2 = {cfg.detparam2};
        catch
            error('Widefield detectors must have detparam1 and detparam2');
        end
        if isfield(cfg,'detpattern')
            tempwf.srcpattern = {cfg.detpattern};
        end
        if isfield(cfg,'detweight')
            tempwf.srcweight = {cfg.detweight};
        end
        widesrcid = containers.Map({''},{tempwf});clear tempwf
    end
end

if isa(cfg.prop,'containers.Map')
    prop = cfg.prop;
else
    prop = containers.Map({''},{cfg.prop});
end

if (~isequal(cell2mat(widesrcid.keys),cell2mat(prop.keys)) && (~isempty(cell2mat(widesrcid.keys))))
    error('Keys for prop and widesrcid must be compatible')
elseif (~isequal(cell2mat(widesrcid.keys),cell2mat(prop.keys)) && (~isempty(cell2mat(prop.keys))))
    temp = widesrcid('');
    clear widesrcid
    widesrcid = containers.Map();
    wv = prop.keys;
    for ii = 1:length(prop.keys)
        widesrcid(wv{ii}) = temp;
    end            
end

widesrc = [];
wavelengths = prop.keys;
wfsrcmapping = containers.Map();
allidwf = [];

for wv = wavelengths
    wideparam = widesrcid(wv{1});
    op = prop(wv{1});
    
    z0 = 1/(op(2,1)+op(2,2)*(1-op(2,3)));
    
    srcmapping = [];    
    
    for wideidx = 1:length(wideparam.srcid)
        clear srcbc rhs
        
        srcid = wideparam.srcid{wideidx};
        allidwf = [allidwf;srcid];
        srctype = wideparam.srctype{wideidx};
        srcparam1 = wideparam.srcparam1{wideidx};
        srcparam2 = wideparam.srcparam2{wideidx};
        
        if(strcmp(srctype,'pattern'))
            srcpattern = wideparam.srcpattern{wideidx};
        end
        if isfield(wideparam,'srcweight')
            srcweight = wideparam.srcweight{wideidx};
        end
        
        srcpos = sources(srcid,:);

        switch srctype
            case {'planar','pattern','fourier'}
                    ps=[srcpos; srcpos+srcparam1(1:3); ...
                        srcpos+srcparam1(1:3)+srcparam2(1:3); srcpos+srcparam2(1:3); srcpos];

                    pnode=cfg.node;
                    pface=cfg.face;

                    % if src is colimated (default), sink it by 1/mus'
                    if(~isfield(cfg,'iscolimated') || cfg.iscolimated) 
                        sinkplane=srcdir; % the plane where sinked planar source is located as [A,B,C,D] where A*x+B*y+C*z+D=0
                        sinkplane(4)=-sum(srcdir.*(srcpos+srcdir*z0));
                        [cutpos,cutvalue,facedata,elemid,nodeid]=qmeshcut(cfg.elem,cfg.node,zeros(size(cfg.node,1),1),sinkplane);
                        pnode=cutpos;
                        idx=find(facedata(:,3)~=facedata(:,4));
                        pface=facedata(facedata(:,3)==facedata(:,4),1:3);
                        pface=[pface;facedata(idx,[1 2 3]);facedata(idx,[1 3 4])];
                    end
                    c0=meshcentroid(pnode,pface);
                    newnode=rotatevec3d([c0; ps],srcdir(1:3));
                    srcpoly=newnode(end-4:end, 1:2);
                    [isin,ison]=inpolygon(newnode(1:end-5,1),newnode(1:end-5,2),srcpoly(:,1),srcpoly(:,2));
                    isin=isin | ison;
                    idx=find(isin);
                    if(~isempty(idx)) % the below test only works for convex shapes
                            AB=pnode(pface(idx,2),1:3)-pnode(pface(idx,1),1:3);
                            AC=pnode(pface(idx,3),1:3)-pnode(pface(idx,1),1:3);
                            N=cross(AB',AC')';
                            dir=sum(N.*repmat(srcdir(:)',size(N,1),1),2);
                            if(exist('sinkplane','var'))
                                dir(dir>0)=-dir(dir>0);
                            end
                            if(all(dir>=0))
                                    error('please reorient the surface triangles');
                            end
                            srcbc=zeros(1,size(pface,1));
                            srcbc(idx(dir<0))=1;
                            pbc=newnode(idx(dir<0),1:2);

                            dp=pbc-repmat(srcpoly(1,:),size(pbc,1),1);
                            dx=srcpoly(2,:)-srcpoly(1,:);
                            dy=srcpoly(4,:)-srcpoly(1,:);
                            nx=dx/norm(dx);
                            ny=dy/norm(dy);

                            bary=[sum(dp.*repmat(nx/norm(dx),size(dp,1),1),2), sum(dp.*repmat(ny/norm(dy),size(dp,1),1),2)];
                            bary(bary<0)=0;
                            bary(bary>=1)=1-1e-6;

                            if(exist('srcpattern','var'))
                                    srcpattern=permute(srcpattern,[3 1 2]);
                                    pdim=size(srcpattern);
                                    patsize=pdim(1);
                                    srcbc=repmat(srcbc,patsize,1);
                                    for i=1:patsize
                                            srcbc(i,idx(dir<0))=srcpattern(sub2ind(pdim, i*ones(size(bary,1),1), floor(bary(:,1)*pdim(2))+1, floor(bary(:,2)*pdim(3))+1));
                                    end
                            elseif(strcmp(srctype,'fourier'))
                                    kx=floor(srcparam1(4));
                                    ky=floor(srcparam2(4));
                                    phi0=(srcparam1(4)-kx)*2*pi;
                                    M=1-(srcparam2(4)-ky);
                                    srcbc=repmat(srcbc,kx*ky,1);
                                    for i=1:kx
                                            for j=1:ky
                                                srcbc((i-1)*ky+j,idx(dir<0))=0.5*(1+M*cos((i*bary(:,1)+j*bary(:,2))*2*pi+phi0))';
                                            end
                                    end
                            end
                    else
                            error('source direction does not intersect with the domain');
                    end
            otherwise
                    error('this source type is not supported');
        end

        % at this point, srcbc stores the J- at each surface triangle (or sinked triangles)

        if isa(cfg.reff,'containers.Map')
            Reff=cfg.reff(wv{1});
        else
            Reff = cfg.reff;
        end
        maxbcnode=max(pface(:));
        if(exist('nodeid','var'))
            nodeweight=nodeid(:,3);
            nodeid=nodeid(:,1:2);
            maxbcnode=max(max(nodeid(pface,:)));
            parea=elemvolume(pnode,pface);
        else
            parea=cfg.area;
        end

        % 1/18 = 1/2*1/9, where 2 comes from the 1/2 in ls=(1+Reff)/(1-Reff)/2*D,
        % and 1/9 = (1/6+1/12+1/12)/3, where A/6 is <phi_i,phi_j> when i=j, and
        % A/12 is i!=j
        Adiagbc= parea(:)*((1-Reff)/(18*(1+Reff)));%ones(size(parea(:),1),1)*((1-Reff)/(18*(1+Reff)));
        Adiagbc=repmat(Adiagbc,1,size(srcbc,1)).*(srcbc');
        rhs=sparse(size(cfg.node,1),size(srcbc,1));

        for i=1:size(srcbc,1)
            if(exist('nodeid','var'))
                allnodes=nodeid(pface,:);
                rhs(1:maxbcnode,i)=sparse(allnodes(:), 1, [repmat(Adiagbc(:,i),3,1).*nodeweight(pface(:));repmat(Adiagbc(:,i),3,1).*(1-nodeweight(pface(:)))]);
            else
                rhs(1:maxbcnode,i)=sparse(cfg.face(:), 1, repmat(Adiagbc(:,i),1,3));
            end
            wsrc=1;
            if(exist('srcweight','var') && numel(srcweight)==size(srcbc,1))
                wsrc=srcweight(i);
            end
%             rhs(rhs(:,i)<1e-3*max(rhs(:,i)),i) = 0;
%             rhs(:,i)=rhs(:,i)*(wsrc/abs(sum(rhs(:,i))));
            rhs(:,i)=rhs(:,i)*(wsrc/sum(abs(rhs(:,i))));
%             rhs(:,i)=rhs(:,i)*(wsrc/sum(rhs(:,i)));
        end
        srcbc=rhs.';
        if exist('srcpattern','var')
            patsize = size(srcpattern,1);
        else
            patsize = 1;
        end
        
        indices = [(size(widesrc,1) + 1) (size(widesrc,1) + patsize)];
        widesrc = [widesrc; full(srcbc)];
        srcmapping = [srcmapping;srcid indices];
    end
    wfsrcmapping(wv{1}) = srcmapping;
end
sources(unique(allidwf),:) = [];

if (length(wfsrcmapping) == 1)
    wfsrcmapping = wfsrcmapping(wavelengths{1});
end

if (~isdet)
    cfg.widesrc = widesrc;
    cfg.wfsrcmapping = wfsrcmapping;
    cfg.srcpos = sources;
else
    cfg.widedet = widesrc;
    cfg.wfdetmapping = wfsrcmapping;
    cfg.detpos = sources;
end
    