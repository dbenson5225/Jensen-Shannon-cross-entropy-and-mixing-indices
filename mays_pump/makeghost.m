function [ghost] = makeghost(N,X,dx,xmin,xmax)

    ndim=size(X,2);
    eps=1e-10;
    ghost=zeros(4*N,ndim);
    
    if ndim==2
    ghost(1:N,1)=X(:,1)-dx(1); ghost(1:N,2)=X(:,2);
    ghost(N+1:2*N,1)=X(:,1)+dx(1); ghost(N+1:2*N,2)=X(:,2);
    ghost(2*N+1:3*N,1)=X(:,1); ghost(2*N+1:3*N,2)=X(:,2)-dx(2);
    ghost(3*N+1:4*N,1)=X(:,1); ghost(3*N+1:4*N,2)=X(:,2)+dx(2);
    
    ghost(ghost(:,1)<xmin(1)+dx(1),1)=xmin(1)+dx(1)+eps;
    ghost(ghost(:,1)>xmax(1)-dx(1),1)=xmax(1)-dx(1)-eps;
    %apple1=X(:,2)<=xmin(2);
    %apple2=X(:,2)>=xmax(2);
    ghost(ghost(:,2)<=xmin(2),2)=xmin(2)+eps;
    ghost(ghost(:,2)>=xmax(2),2)=xmax(2)-eps;
    
   % dxvec=ghost(:,1)-X(:,1);
   % dyvec=ghost(:,2)-X(:,2);

end

