clear all
omega=10;
ntsteps =200; 
D = 1e-4;
dr = 1e-3;
nnodes=ceil(omega/dr);
methodd = 'implicit';
num = 1e2;
A=zeros(nnodes,1); B=zeros(nnodes,1);
KLAB=zeros(ntsteps,2); KLBA=zeros(ntsteps,2);
SDR_1=zeros(ntsteps,3); SDR_2=zeros(ntsteps,2);
SYMM=zeros(ntsteps,2); SYMMnew=zeros(ntsteps,6);
shannon=zeros(ntsteps,1); react=zeros(ntsteps,1);

r=(-dr/2+dr*(1:nnodes))';  % vector of radial distance to cell centers.
half_gap=dr;

%%%%%%%%%%
% delta IC at origin
%A(1:1)=1;  A=A/(pi*dr*dr*sum(A)); 
%B(2:end)=1; 
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%
% Mays-like IC
A(1:round(nnodes/2))=1;
B(round(nnodes/2)+1:end)=1; 
%%%%%%%%%%%%%%%%%%%%%

A=A/(2*pi*dr*sum(r.*A)); 
B=B/(2*pi*dr*sum(r.*B));
min(A)
min(B)
C=A-B;
pore = ones(nnodes, 1);
eps=1e-200;
dt=.01; tott=0.0;

for n=1:ntsteps
   m=0.5*(A+B);
   KLAB(n,1)=tott; KLBA(n,1)=tott;
   SDR_1(n,1)=tott; SDR_2(n,1)=tott;
   appleA=find(A>eps);
   appleB=find(B>eps);
   apple=find(A>eps & B>eps);
   Anorm=2*pi*dr*sum(r(apple).*A(apple)) 
   Bnorm=2*pi*dr*sum(r(apple).*B(apple))
%   apple=find
   KLAB(n,2) = (dr/Anorm)*sum(r(apple).*A(apple).*log(Bnorm*A(apple)./(Anorm*B(apple))));   % Assym KLD #1
   KLBA(n,2) = (dr/Bnorm)*sum(r(apple).*B(apple).*log(Anorm*B(apple)./(Bnorm*A(apple))));   % Assym KLD #2
   SYMM(n,2) = KLAB(n,2)+KLBA(n,2);                         % Symmetric KLD
   SYMMnew(n,2) = -log(dr)-dr*sum(r(appleB).*A(appleB).*log(B(appleB)));
   SYMMnew(n,3) = -log(dr)-dr*sum(r(appleA).*B(appleA).*log(A(appleA)));
   SYMMnew(n,4) = -log(dr)-dr*sum(r(A>0).*A(A>0).*log(A(A>0)));
   SYMMnew(n,5) = -log(dr)-dr*sum(r(B>0).*B(B>0).*log(B(B>0)));
   SYMMnew(n,6) = -(dr/Anorm)*sum(r(apple).*A(apple).*log(B(apple)/Bnorm)) - ...
                   (dr/Bnorm)*sum(r(apple).*B(apple).*log(A(apple)/Anorm)) + ...
                   (dr/Anorm)*sum(r(apple).*A(apple).*log(A(apple)/Anorm)) + ...
                   (dr/Bnorm)*sum(r(apple).*B(apple).*log(B(apple)/Bnorm));

   %shannon = 0.5*dx*sum(f(fgood).*log(f(fgood)./m(fgood))) + ...
   %          0.5*dx*sum(g(ggood).*log(g(ggood)./m(ggood)))

   Anorm=2*pi*dr*sum(r(appleA).*A(appleA)) 
   Bnorm=2*pi*dr*sum(r(appleB).*B(appleB))
   mnormA=2*pi*dr*sum(r(appleA).*m(appleA))
   mnormB=2*pi*dr*sum(r(appleB).*m(appleB))

   shannon(n) = 2*pi*(0.5*dr*sum(r(appleA).*A(appleA).*log(A(appleA)./m(appleA)) ) + ...
                      0.5*dr*sum(r(appleB).*B(appleB).*log(B(appleB)./m(appleB)) ) );


   %shannon(n) = 2*pi*((0.5/Anorm)*dr*sum(r(appleA).*A(appleA).*log((mnormA/Anorm)*A(appleA)./m(appleA))) + ...
   %             (0.5/Bnorm)*dr*sum(r(appleB).*B(appleB).*log((mnormB/Bnorm)*B(appleB)./m(appleB))) );
   
   
   % SDR_1 is int(grad C)D grad(c)
   % SDR_2 is int c^2  (later taking -0.5 d(SDR_2)/dt)

 %  SDR_1(n,2) = (D/dr)*(r.*diff(A-B))'*diff(A-B);
 %  SDR_1(n,3) = (D/dr)*(r.*diff(A))'*diff(B);
 %  SDR_2(n,2) = dx*sum((A-B).^2);
   %SDR_2(n,2) = dr*sum(r.*(A-B).^2);
   %SDR_2(n,2) = 2*pi*dr*sum(r.*A.^2);
   SDR_2(n,2) = 2*pi*dr*sum(r.*(A-B).^2);
   

   figure(1)
   plot(r,A,'r'); hold on 
   plot(r,B,'b');
   %plot(rad,C,'k--');
   %plot(r,exp(r.^2/(-4*D*tott))./(4*pi*D*tott), 'k' )
 %  axis([0 omega 0 1.5])
   hold off
   drawnow
%pause
% Analytic solutions
%    x1=omega/2 - dx*half_gap;
%    x2=omega/2 + dx*half_gap;
%    A=exp((dx*(1:nnodes)'-x1).^2/(-4*D*tott))/sqrt(4*pi*D*tott);
%    B=exp((dx*(1:nnodes)'-x2).^2/(-4*D*tott))/sqrt(4*pi*D*tott);
%    A=A/(dx*sum(A)); B=B/(dx*sum(B));
%    C=A-B;
% Numerical solutions
   
   [A] = diff_1d_radial_spdiags(pore,A,D,dt,dr,methodd);
   [B] = diff_1d_radial_spdiags(pore,B,D,dt,dr,methodd);
   [C] = diff_1d_radial_spdiags(pore,C,D,dt,dr,methodd);

   massA=(2*pi*dr*sum(r.*A))
   massB=(2*pi*dr*sum(r.*B))

   tott=tott+dt  
   dt=1.1*dt;

    %pause
end
hold off

shannon=min(log(2),shannon);
SYMM(:,1)=KLAB(:,1);  SYMMnew(:,1)=KLAB(:,1); time=SYMM(:,1);
factor=2;  % factor of 1 means scale SKLD down: 8 means scale JS up.
mixtime=((dr*2*half_gap)^2)/(8*2*D/factor)
SDR_2(1:end-1,1)=0.5*(SDR_2(1:end-1,1)+SDR_2(2:end,1)); 
SDR_2(1:end-1,2)=-0.5*(diff(SDR_2(:,2))./diff(SDR_2(:,1)) );
SDR_2(end,:)=[];

figure(2)
loglog(time,(factor/8)*(SYMMnew(:,2)+SYMMnew(:,3)),'r-'); hold on
axis([1e-3 1e4 1e-6 tott])
plot(time,(factor/8)*(SYMMnew(:,4)+SYMMnew(:,5)),'b-');
plot(time,(factor/8)*SYMMnew(:,6),'ko');
plot(time,factor*shannon,'g+')
%plot(time,mixtime./SYMMnew(:,1),'k--')
plot(SDR_2(:,1),SDR_2(:,2),'pentagram')
plot(SDR_1(:,1),(2/sqrt(128*pi*D))*SDR_1(:,1).^-1.5,'b-.')
legend('Cross','Self','SKLD/4','2*JS','Eq. (21)','SDR')
hold off

figure(22)
%loglog(SYMMnew(:,1),(factor/8)*(SYMMnew(:,2)+SYMMnew(:,3)),'r-'); hold on
%plot(SYMMnew(:,1),(factor/8)*(SYMMnew(:,4)+SYMMnew(:,5)),'b-');
loglog(SYMMnew(:,1),(factor/8)*SYMMnew(:,6),'ko');
hold on
plot(time,factor*shannon,'g+')
%plot(time,mixtime./time,'k--')
plot(time,SDR_1(:,2),'bd')
plot(SDR_1(:,1),SDR_1(:,3),'rd')
plot(SDR_2(:,1),SDR_2(:,2),'pentagram')
plot(SDR_1(:,1),(2/sqrt(128*pi*D))*SDR_1(:,1).^-1.5,'b-.')
axis([1e-3 1e5 1e-6 tott])
legend('SKLD/4','2*JS','Eq. (21)','SDR')
hold off


figure(5)
dtimes=0.5*(SDR_1(1:end-1,1)+SDR_1(2:end,1));
loglog(SDR_2(:,1),SDR_2(:,2),'-o'); hold on
axis([1e-3 1e5 1e-6 tott])
plot(SDR_1(:,1),SDR_1(:,2),'bd'); 
%loglog(dtimes,1*dtimes.^-1.5,'k--')
plot(dtimes,(1/sqrt(128*pi*D))*dtimes.^-1.5,'k--')
plot(dtimes,-diff(KLAB(:,2))./diff(KLAB(:,1)),'-+');
legend('SDR-2','SDR-1','t^-3/2','SKLD')
hold off
%loglog(SDR_2(:,1),SDR_2(:,2))


mix_JS=(1/(1-exp(-factor*log(2))))*(exp(-factor*shannon) - exp(-factor*log(2)));

figure(6)
semilogx(time,exp(-(factor/8)*SYMMnew(:,6)),'ok')
hold on
plot(time, mix_JS, 'r+')
plot(time,exp(-mixtime./time),'k--')
axis([1e-3 1e4 0 1])
xlabel('Time'); ylabel('Mixing Index M(t) (dimensionless)')
legend('F-D reflecting (SKLD)','F-D reflecting (JS)','Analytic infinite domain','Location','NW')
hold off

dmixdt=diff(exp(-(factor/8)*SYMMnew(:,6)))./diff(SYMMnew(:,1));
dmixdt_S=diff(mix_JS)./diff(SYMMnew(:,1));
ave_t = 0.5*(SYMMnew(1:(end-1),1)+SYMMnew(2:end,1))
dSDRdt=diff(exp(-SDR_1(:,2)))./diff(SDR_1(:,1));

figure(7)
%plot(ave_t/(separation^2/(16*D)),dmixdt,'ko',ave_t/(separation^2/(16*D)),dmixdt_S,'r+')
plot(ave_t,dmixdt,'ko',ave_t,dmixdt_S,'r+')
hold on
plot(SYMMnew(:,1),(mixtime*time.^-2).*exp(-mixtime./time),'k--')
%plot(ave_t/(separation^2/(16*D)),dSDRdt,'bd')
xlabel('Time'); ylabel('Mixing Rate dM(t)/dt (1/T)')
legend('F-D reflecting (SKLD)','F-D reflecting (JS)','Analytic SKLD infinite domain','d(SDR)/dt','Location','NE')
axis([0 1e2 0 5e-2])
hold off

figure(8)
loglog(ave_t,dmixdt,'co',ave_t,dmixdt_S,'k-')
hold on
plot(SDR_2(:,1),SDR_2(:,2),'pentagram')

%loglog(time,(mixtime*time.^-2).*exp(-mixtime./time),'k--')
%plot(ave_t/(separation^2/(16*D)),dSDRdt/mixtime,'bd')
xlabel('Time'); ylabel('Mixing Rate dM(t)/dt (1/T)')
legend('F-D reflecting (SKLD)','F-D reflecting (JS)','SDR','d(SDR)/dt','Location','NE')
axis([1e-3 1e6 1e-9 1e0])
%axis([0 2e2 0 2e-3])
hold off