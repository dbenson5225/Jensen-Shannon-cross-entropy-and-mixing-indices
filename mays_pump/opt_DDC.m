function [dxfactors] = opt_DDC(A,ndim,numcoresmin,numcoresmax)

%%%%% Only works in 2-D right now!
%xrange=zeros(ndim,1); aspect=zeros(ndim,1);
for i=1:ndim
    xrange(i)=max(A(:,i))-min(A(:,i));
end

if     ndim==1 
    aspect(1)=1.;
elseif ndim==2
    aspect(1)=xrange(1)/xrange(2);
elseif ndim==3
    aspect(1)=xrange(1)/xrange(2);
end

difference=1000;
numcores=numcoresmax;
while difference>1.1 & numcores>=numcoresmin
K = 1:ceil(sqrt(numcores));
D1 = K(rem(numcores,K)==0);
D1 = [D1 sort(numcores./D1)];
D2 = numcores./D1;

ratio=D1./D2;
[best,pos]=min(abs(ratio-aspect(1)));
if aspect(1)<1
f1=min(D1(pos),D2(pos));
f2=max(D1(pos),D2(pos));
else
f2=min(D1(pos),D2(pos));
f1=max(D1(pos),D2(pos));
end
bestratio = f1/f2;
difference=abs( (1-(bestratio-aspect(1))/aspect(1)));

numcores=numcores-1;
end
numcores=numcores+1;
dxfactors=[f1, f2];
disp(['Number of subdivisions = ',num2str(numcores),';  in [x,y] = ',num2str(dxfactors)]);
%pause

%figure(99)
%plot(D1,abs((D1./D2)-aspect),'r-o')
%hold on
%plot(D1,abs((D2./D1)-aspect),'b-o')
%axis([0 100 0 10])
%leg1=append('|f_1/f_2 - ',num2str(aspect),'|')
%leg2=append('|f_2/f_1 - ',num2str(aspect),'|')
%legend(leg1,leg2)
%titstr=append('numcores = ',num2str(numcores),'  aspect ratio = ',num2str(aspect));
%title(titstr)
%xlabel('factors');ylabel('Error')
%hold off