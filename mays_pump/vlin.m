function [vpart] = vlin(X,vx,vy,dx,xmin,xmax)
  
%  Only in 2-D right now!!

    eps=1e-20; N=size(X,1); ndim=size(X,2);
    nnx=size(vx,2)-1;   nny=size(vy,1)-1;
    
    X(X(:,1)<xmin(1),1)=xmin(1)+eps;
    X(X(:,1)>xmax(1),1)=xmax(1)-eps;
    apple1=X(:,2)<xmin(2);
    apple2=X(:,2)>xmax(2);
    X(apple1,2)=xmin(1)+eps;
    X(apple2,2)=xmax(2)-eps;
    
    blahx = X(:,1)/dx(1);
    blahy = X(:,2)/dx(2);
        
    boxx = floor(blahx) + 1;
    boxy = floor(blahy) + 1;

    boxx(boxx<1)=1;
    boxx(boxx>nnx)=nnx;
    
    xdist = blahx - boxx + 1;
    ydist = blahy - boxy + 1;
        
    lindexcell=nny*(boxx-1)+boxy;
    lindexvx=lindexcell; 
    lindexvy=(nny+1)*(boxx-1)+boxy;
    
    vxpart = (1.0-xdist).*vx(lindexvx) + xdist.*vx(lindexvx+nny);
    vypart = (1.0-ydist).*vy(lindexvy) + ydist.*vy(lindexvy+1); 
    
    vxpart(apple1)=0; vxpart(apple2)=0;
    vypart(apple1)=0; vypart(apple2)=0;
 
    vpart(:,1)=vxpart;
    vpart(:,2)=vypart;
    
end

