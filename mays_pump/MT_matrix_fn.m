function [A] = MT_matrix_fn(A,ndim,nspec,xmin,xmax,nxmt,DMT,dt,ds)
%%%%%
%   A comes in as % x,(y,z),A,(B,C),mobile,alive ...  
%%% I.e., size(A)=  (np,ndim+nspec+2)
dofig=false;
%gpu=gpuDeviceCount;      %  See if there's a gpu
gpu=0;                   %  Turn off gpu
%define numerical and physical parameters
np=size(A,1);
eps=1e-10;
particlecell=ones(np,3);
betaval=1.0;   % The SPH beta value
            
if(gpu>0) 
    A=gpuArray(A);
    xmin=gpuArray(xmin); xmax=gpuArray(xmax);
end

sigma = sqrt(4*max(DMT)*dt);  % Char. diffusion dist.
    % Do the mass-transfer
    % First chop up space: ideally each chunk is 3*sigma
if(gpu>0); nxmt=gpuArray(nxmt); end;
for i=1:ndim
    pad(i)=5*sigma; 
    if(gpu>0); pad=gpuArray(pad); end;
    %nxmt(i)=floor((xmax(i)-xmin(i))/dxmt(i)); 
    dxmt(i)=(xmax(i)-xmin(i))/nxmt(i);  % for constant dx
    for m=1:nxmt(i)    % This can allow variable dx's later (i.e. kdtree)
        xlo(i,m)=xmin(i)+(m-1)*dxmt(i);
        xhi(i,m)=xmin(i)+m*dxmt(i);
    end
end
% pad(ndim+1:3)=1

  % Figure out the mass-transfer i,j,k (cell) of each particle
  for i=1:ndim
      particlecell(1:np,i)=ceil((A(:,i)-xmin(i))./dxmt(i));
  end
  
  % Build matrix of separations including ghosts from adjoining cells
  %for i=1:ndim
  %    sz(i)=nxmt(i);
  %end

  massnew=A(:,ndim+1:ndim+nspec);  % Pull masses from A; push to massnew!
  masstemp=massnew;     % This holds temporary masses (incl. ghosts)
  lidx=1:prod(nxmt);   
  [is js ks]=ind2sub(nxmt,lidx);

  for m=1:prod(nxmt)            % go to every cell  (cell loop)
      ijk(1)=is(m); ijk(2)=js(m); ijk(3)=ks(m);   
      lims=zeros(3,2);  if(gpu>0); lims=gpuArray(lims); end;
      for i=1:ndim
         lims(i,1)=xlo(i,ijk(i))-pad(i); lims(i,2)=xhi(i,ijk(i))+pad(i);  %xlo and xhi
      end

     if(ndim==1) idxlocal=particlecell(:,1)==ijk(1); 
     elseif(ndim==2) idxlocal=particlecell(:,1)==ijk(1) & particlecell(:,2)==ijk(2);
     elseif(ndim==3) idxlocal=particlecell(:,1)==ijk(1) & ...
                     particlecell(:,2)==ijk(2) & particlecell(:,3)==ijk(3);
     end  

     if(ndim==1)      idxghost=A(:,1)>lims(1,1) & A(:,1)<=lims(1,2);
     elseif(ndim==2)  idxghost=A(:,1)>lims(1,1) & A(:,1)<=lims(1,2) & ...
                               A(:,2)>lims(2,1) & A(:,2)<=lims(2,2); 
     elseif(ndim==3)  idxghost=A(:,1)>lims(1,1) & A(:,1)<=lims(1,2) & ...
               A(:,2)>lims(2,1) & A(:,2)<=lims(2,2) & ... 
               A(:,3)>lims(3,1) & A(:,3)<=lims(3,2);
     end

if(dofig)
    figure(44)
    plot(A(idxghost,1),A(idxghost,2),'+')
    xlim([xmin(1), xmax(1)]);
    ylim([xmin(2), xmax(2)]); 
    xticks(xmin(1):dxmt(1):xmax(1)); 
    yticks([xmin(2):dxmt(2):xmax(2)]); 
    grid on; hold on;
    plot(A(idxlocal,1),A(idxlocal,2),'o')
    drawnow; hold off; 
    %pause;
end

  if(sum(idxghost)>1 & sum(idxlocal)>0)
      X=A(idxghost,:);
      npad=sum(idxghost);
      Dterm=DMT(idxghost,1); unity=ones(npad,1);
      Dave=0.5*(Dterm*unity' + unity*Dterm');  % D average matrix (symmetric) between particles

      P=squareform( pdist(X(:,1:ndim))  );
      P = exp((P.*P)./(-4*Dave*dt));

%%%% Traditional row and column sum normalized (run iii >= 1 as desired for horn-sinkhoff)
     for iii=1:1
%        colsum=sum(P); 
        rowsum=sum(P,2);
%        P = P ./ (0.5*(colsum+rowsum));
        P = P ./ (0.5*(rowsum'+rowsum));
        rowsum=sum(P,2);
        P=P+diag(1-rowsum);
     end
%      rowsum=sum(P,2);
%      P=P+diag(1-rowsum);
%%%%  end traditional method

%%%% Diagonal adjust only
     % prefactor=ds./((4*pi*Dave*dt).^(ndim/2));
     % P = prefactor.*P;
     % rowsum=sum(P,2);
     % max(rowsum)
     % min(rowsum)
     % P = P + diag(1-rowsum);
%%%% end diagonal adjust only

     masstemp(idxghost,:)= P * X(:,ndim+1:ndim+nspec);  % keep new mass on all ghosts
     massnew(idxlocal,:)=masstemp(idxlocal,:);       % push only local masses up


%   if(length(idxlocal)>0)   % any local particles?
%       [idx,r]=rangesearch(A(idxghost,1:ndim),A(idxlocal,1:ndim),6*sigma);
%       for i=1:len(idxlocal)
%           idxlocal
%           partnow=idxlocal(i)
%           DMTnow=DMT(partnow)
%           Dave=0.5*(DMT(partnow)+(DMT(idxnow))');
%           size(Dave)
%           prefactor=ds./((4*pi*Dave*dt).^(ndim/2));
%           size(prefactor)
%           size(r{i})
%           pause
%           numnear=length(idxnow);
%           if numnear>0
%             Prow=prefactor.*exp( (r{i}.*r{i}) ./(-4*DMT*dt) );
% %            if(abs(Prow(1)-prefactor)>1e-10)
% %                Prow(1)
% %            end
%             rowsum=sum(Prow);
%             Prow(1)=Prow(1)-diag(rowsum-1);
%             dM=beta*(Prow*A(idxghost(idxnow),7:end)-A(partnow,7:end));
%             massnew(partnow,:)=massnew(partnow,:)+dM;  % keep new mass on all ghosts
% 
%           end
%       end
%   end



%  Assumes beta = 1, which is correct
%      masstemp(idxghost,:)= P * massnew(idxghost,:);  % keep new mass on all ghosts
    
  end  % more than 1 particle if statement

  end  % cell loop

  A(:,ndim+1:ndim+nspec)=massnew;

%  if(kappa<1)  % do simple random walks
%      if gpu>0
%         rvec=randn(np,ndim,'gpuArray');
%      else
%         rvec=randn(np,ndim);
%      end
%      A(:,1:ndim)=A(:,1:ndim)+sqrt(2*DRW*dt)*rvec;
%      A(A(:,1)<xmin(1),1)=-A(A(:,1)<xmin(1),1);
%      A(A(:,2)<xmin(2),2)=-A(A(:,2)<xmin(2),2);
%      A(A(:,1)>xmax(1),1)=2*xmax(1)-A(A(:,1)>xmax(1),1);
%      A(A(:,2)>xmax(2),2)=2*xmax(2)-A(A(:,2)>xmax(2),2);
%
%  end  simple random walks loop

