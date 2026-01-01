function [A] = diff_1d_spdiags(pore,A,D,dt,dx,methodd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% pore delineates active transport areas (with ones for active and zeros
% for inactive
nn=length(A);
r=D*dt/dx^2;
%ip=1-pore;
    
if(strcmpi(methodd,'explicit'))
A=(1-2*r)*A+r*(circshift(A,[-1])+circshift(A,[1]));
end
if(strcmpi(methodd,'implicit'))
%   Amat=diag((1+2*r)*ones(nn-2,1)) + diag(-r*ones(nn-3,1),-1) +  diag(-r*ones(nn-3,1),1);
  Amat = spdiags([-r*ones(nn,1) 1+2*r*ones(nn,1) -r*ones(nn,1)], -1 : 1, nn, nn);
  Amat(1)=Amat(1)-r; Amat(nn*nn)=Amat(nn*nn)-r;   % No-flow boundary
  %A(2)=A(2)+r*A(1);  A(nn-1)=A(nn-1)+r*A(nn);  % Constant values on bdys
  A=Amat\A;
% A=inv(Amat)*A
end
        
A=A.*pore;  
%B=B.*pore;
end

