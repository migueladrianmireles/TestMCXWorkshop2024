function [Lmat,reconprior] = rbmakeL (cfg,recon,prior,alpha,beta)

if (nargin < 4)
    alpha = 0.2;
end
if (nargin < 5)
    beta = 1.1;
end

tic
if (size(prior,1) == size(cfg.node,1))
%     [f2rid,f2rweight] = tsearchn(recon.node,recon.elem,cfg.node);
    reconprior = meshremap(prior,recon.mapid,recon.mapweight,recon.elem,size(recon.node,1));    
    reconprior = reconprior./sum(reconprior,2);
    reconprior(isnan(reconprior) | isinf(reconprior)) = 0;
else
    reconprior = prior;
end


%Lmat = zeros(size(reconprior,1),size(reconprior,1));

%for ii = 1:size(reconprior,1)
%    for jj = 1:size(reconprior,1)
%        Lmat(ii,jj) = sum(abs(reconprior(ii,:) - reconprior(jj,:)));        % L1 Norm? Shouldn't this be L2?
%         Lmat(ii,jj) = sqrt(sum((reconprior(ii,:) - reconprior(jj,:)).^2)); 
%        Lmat(jj,ii) = Lmat(ii,jj);
        
%    end
%     fprintf(['looping ' num2str(ii) ' - ' num2str(toc) '\n']);       %   Debugging
%end

%Lmat = zeros(size(reconprior,1),size(reconprior,1));
tic
Lmat = squeeze(sum(abs(permute(reconprior,[3 2 1]) - reconprior),2));

Lmat(Lmat.^2 < ((alpha*size(reconprior,2)).^2)) = -alpha - Lmat(Lmat.^2 < ((alpha*size(reconprior,2)).^2))./size(reconprior,2);
Lmat(Lmat >= alpha*size(reconprior,2))          = 0;
Lmat(1:size(Lmat,2)+1:end)                      = 0;

% dmat = sqrt(abs(sum(Lmat,1)));
% dmat(dmat < 1e-5) = 1e-5;
% dmat = 1./dmat;node
% Lmat = Lmat .* (1/beta) .* repmat(dmat,size(reconprior,1),1) .*repmat(dmat',1,size(reconprior,1));

dd          = abs(sum(Lmat));
dd(dd == 0) = 1;
Lmat        = Lmat./(beta.*sqrt(dd'*dd));

Lmat(1:size(Lmat,2)+1:end) = 1;
toc

fprintf('Computing L matrix - %f seconds\n',toc);
end