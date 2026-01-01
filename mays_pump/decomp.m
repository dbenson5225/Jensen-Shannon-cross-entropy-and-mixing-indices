function [B] = decomp(D)
N=size(D,1);
si=zeros(size(D));
%  Assumes D = [a b; c d] is a,b,c,d in columns

det=D(:,1).*D(:,4)-D(:,2).*D(:,3);
s=sqrt(det);
si(:,1)=s; si(:,4)=s;
tr=D(:,1)+D(:,4);
t=sqrt(2*s+tr);
B=(1./t).*(D+si);
B=real(B);
end

