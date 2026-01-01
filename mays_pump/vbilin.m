function [vpart] = vbilin(X,vx,vy,dx,xmin,xmax)
    ndim=size(X,2);
    eps=1e-10; N=size(X,1); vpart=zeros(N,2);
    nnx=size(vx,2)-1;   nny=size(vy,1)-1;
 
    fracdistx = X(:,1)/dx(1);
    fracdisty = X(:,2)/dx(2);
        
    boxx = floor(fracdistx) + 1;
    boxy = floor(fracdisty) + 1;

    boxx(boxx<1)=1;
    boxx(boxx>nnx)=nnx;
    
    xdist = fracdistx - boxx + 1;
    ydist = fracdisty - boxy + 1;
        
    lindexcell=nny*(boxx-1)+boxy;
    lindexvx=lindexcell; 
    lindexvy=(nny+1)*(boxx-1)+boxy;
    
    whichy=ones(N,1);  whichx=ones(N,1);
    whichy(ydist<0.5)=-1;  
    
    whichy(fracdisty<0.5)=0; whichy(fracdisty>nny-0.5)=0;  % These two are near no-flow
    
    whichx(xdist<0.5)=-1;  whichx(fracdistx<0.5)=0;
    
    vxpart = (1.0-xdist).*vx(lindexvx) + xdist.*vx(lindexvx+nny);
    othervx = (1.0-xdist).*vx(lindexvx+whichy) + xdist.*vx(lindexvx+nny+whichy);
    pear1=fracdisty<0.5;   pear2=fracdisty>nny-0.5;   % reflect at BC?
    othervx(pear1)=0; othervx(pear2)=0;
    
    % blah=ydist(pear2)
    % yposition=X(pear2,2)
    
    ydistweight=1.0-abs(ydist-0.5);
    ydistweight(pear1)=2.0*ydist(pear1);   %Adjust for zero-veloc BC at edges
    ydistweight(pear2)=2.0*(1.0-ydist(pear2));
    vxpart = ydistweight.*vxpart + (1.0-ydistweight).*othervx;   
    
    vypart = (1.0-ydist).*vy(lindexvy) + ydist.*vy(lindexvy+1); 
    othervy= (1.0-ydist).*vy(lindexvy+nny*whichx) + ydist.*vy(lindexvy+nny*whichx+1);
    
    %vypart=  (1.0-(xdist-0.5)).*vypart + abs(xdist-0.5).*othervy;
    vypart=  (1.0-(xdist-0.5)).*vypart + (xdist-0.5).*othervy;
    %vypart(apple1)=0.0;  vypart(apple2)=0.0;
     
    vpart(:,1)=vxpart;
    vpart(:,2)=vypart;
%     figure(5)
%     plot(vpart(:,1),X(:,2),'.')
    
end

