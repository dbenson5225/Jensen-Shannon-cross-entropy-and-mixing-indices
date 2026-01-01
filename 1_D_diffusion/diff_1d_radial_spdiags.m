function [A] = diff_1d_radial_spdiags(pore,A,D,dt,dr,methodd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% pore delineates active transport areas (with ones for active and zeros
% for inactive
nn=length(A);
%nn=10; D=1e-3; dt=1; dr=0.1; methodd='implicit'
k=D*dt/dr^2;
%ip=1-pore;

rad=(dr*(0:nn))';  % vector of radius for each intercell bndry
rcell=(-dr/2)+(dr*(1:nn))';

if(strcmpi(methodd,'explicit'))
    A=(1-2*r)*A+r*(circshift(A,[-1])+circshift(A,[1]));
end
if(strcmpi(methodd,'implicit'))
%   Amat=diag((1+2*r)*ones(nn-2,1)) + diag(-r*ones(nn-3,1),-1) +  diag(-r*ones(nn-3,1),1);
%   Amat = spdiags([-r*ones(nn,1) 1+2*r*ones(nn,1) -r*ones(nn,1)], -1 : 1, nn, nn);

  Bin=[ [-k*rad(2:end-1)./rcell(2:end);0]  ...
        1+k*(rad(1:end-1)+rad(2:end))./rcell  ...
        [0;-k*rad(2:end-1)./rcell(1:end-1)] ];
  
  Amat=spdiags(Bin, -1 : 1, nn, nn);
  Amat(1)=Amat(1)-k*rad(1)/rcell(1); Amat(nn*nn)=Amat(nn*nn)-k*rad(nn+1)/rcell(nn);   % No-flow boundary
  %A(2)=A(2)+r*A(1);  A(nn-1)=A(nn-1)+r*A(nn);  % Constant values on bdys
 
%full (Amat) 

  A=Amat\A;
% A=inv(Amat)*A
end
        
A=A.*pore;  
%B=B.*pore;
end

