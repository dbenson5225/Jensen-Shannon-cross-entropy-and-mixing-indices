clear all
omega=10;
ntsteps =200; 
D = 1e-2;
dx = 1e-3;
nnodes=ceil(omega/dx);
methodd = 'implicit';
num = 1e2;
A=zeros(nnodes,1); B=zeros(nnodes,1);
KLAB=zeros(ntsteps,2); KLBA=zeros(ntsteps,2);
SDR_1=zeros(ntsteps,3); SDR_2=zeros(ntsteps,2);
SYMM=zeros(ntsteps,2); SYMMnew=zeros(ntsteps,6);
shannon=zeros(ntsteps,1); react=zeros(ntsteps,1);
separation=1;
half_gap=(separation/2)*(1/dx);
max_mix=(2*half_gap*dx)^2/(4*D);

%A(2*nnodes/10:4*nnodes/10)=1;
%B(6*nnodes/10:8*nnodes/10)=1;
%%%%%%%%%%%%%%%%%%%%%%% IC #1: two deltas
 A(nnodes/2-half_gap)=1;
 B(nnodes/2+half_gap)=1;
%%%%%%%%%%%%%%%%%%%%%%% IC #2: two Heavisides
% A(1:nnodes/2-half_gap)=1;
% B(nnodes/2+half_gap:end)=1;
%%%%%%%%%%%%%%%%%%%%%% IC #3: All random, uncorrelated
 % frac_IC=0.05   % fraction of nodes filled with A or B
 % buffer=floor(nnodes/10000);
 % rand_nodes=buffer+randperm(nnodes-2*buffer,floor(frac_IC*nnodes));
 % A(rand_nodes) = randi(11, floor(frac_IC*nnodes), 1);
 % A(rand_nodes) = 1 + 0.5*(rand(size(rand_nodes))-0.5);
 % rand_nodes=buffer+randperm(nnodes-2*buffer,floor(frac_IC*nnodes));
 % B(rand_nodes) = randi(11, floor(frac_IC*nnodes), 1);
 % B(rand_nodes) = 1 + 0.5*(rand(size(rand_nodes))-0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=A/(dx*sum(A)); B=B/(dx*sum(B));
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
   Anorm=dx*sum(A(apple)); Bnorm=dx*sum(B(apple));
%   apple=find
   KLAB(n,2) = (dx/Anorm)*sum(A(apple).*log(Bnorm*A(apple)./(Anorm*B(apple))));   % Assym KLD #1
   KLBA(n,2) = (dx/Bnorm)*sum(B(apple).*log(Anorm*B(apple)./(Bnorm*A(apple))));   % Assym KLD #2
   SYMM(n,2) = KLAB(n,2)+KLBA(n,2);                         % Symmetric KLD
   SYMMnew(n,2) = -log(dx)-dx*sum(A(appleB).*log(B(appleB)));
   SYMMnew(n,3) = -log(dx)-dx*sum(B(appleA).*log(A(appleA)));
   SYMMnew(n,4) = -log(dx)-dx*sum(A(A>0).*log(A(A>0)));
   SYMMnew(n,5) = -log(dx)-dx*sum(B(B>0).*log(B(B>0)));
   SYMMnew(n,6) = -(dx/Anorm)*sum(A(apple).*log(B(apple)/Bnorm)) - ...
                   (dx/Bnorm)*sum(B(apple).*log(A(apple)/Anorm)) + ...
                   (dx/Anorm)*sum(A(apple).*log(A(apple)/Anorm)) + ...
                   (dx/Bnorm)*sum(B(apple).*log(B(apple)/Bnorm));

   %shannon = 0.5*dx*sum(f(fgood).*log(f(fgood)./m(fgood))) + ...
   %          0.5*dx*sum(g(ggood).*log(g(ggood)./m(ggood)))

   shannon(n) = 0.5*dx*sum(A(appleA).*log(A(appleA)./m(appleA))) + ...
                0.5*dx*sum(B(appleB).*log(B(appleB)./m(appleB)));

   % SDR_1 is int(grad C)D grad(c)
   % SDR_2 is int c^2  (later taking -0.5 d(SDR_2)/dt)

   SDR_1(n,2) = (D/dx)*(diff(A-B))'*diff(A-B);
   SDR_1(n,3) = (D/dx)*(diff(A))'*diff(B);
 %  SDR_2(n,2) = dx*sum((A-B).^2);
   SDR_2(n,2) = dx*sum((A-B).^2);
   
   figure(1)
   plot(dx*(1:nnodes),A,'r'); hold on 
   plot(dx*(1:nnodes),B,'b');
   plot(dx*(1:nnodes),C,'k--');
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
   [A] = diff_1d_spdiags(pore,A,D,dt,dx,methodd);
   [B] = diff_1d_spdiags(pore,B,D,dt,dx,methodd);
   [C] = diff_1d_spdiags(pore,C,D,dt,dx,methodd);

   k_f=0.0;
   A(A<0)=0; B(B<0)=0;
   reaction=k_f*A.*B;

   A=A-reaction; B=B-reaction;
   A(A<0)=0; B(B<0)=0; A=A/(dx*sum(A)) ; B=B/(dx*sum(B));
   react(n)=dx*sum(reaction);

    tott=tott+dt;  dt=1.1*dt;

   % if(tott>max_mix)
   %     pause
   % end
end
hold off

shannon=min(log(2),shannon);
SYMM(:,1)=KLAB(:,1);  SYMMnew(:,1)=KLAB(:,1); time=SYMM(:,1);
factor=2;  % factor of 1 means scale SKLD down: 8 means scale JS up.
mixtime=((dx*2*half_gap)^2)/(8*2*D/factor)
SDR_2(1:end-1,1)=0.5*(SDR_2(1:end-1,1)+SDR_2(2:end,1)); 
SDR_2(1:end-1,2)=-0.5*(diff(SDR_2(:,2))./diff(SDR_2(:,1)) );
SDR_2(end,:)=[];

figure(2)
loglog(time,(factor/8)*(SYMMnew(:,2)+SYMMnew(:,3)),'r-'); hold on
axis([1e-3 1e4 1e-6 tott])
plot(time,(factor/8)*(SYMMnew(:,4)+SYMMnew(:,5)),'b-');
plot(time,(factor/8)*SYMMnew(:,6),'ko');
plot(time,factor*shannon,'g+')
plot(time,mixtime./SYMMnew(:,1),'k--')
plot(SDR_2(:,1),SDR_2(:,2),'pentagram')
plot(SDR_1(:,1),(2/sqrt(128*pi*D))*SDR_1(:,1).^-1.5,'b-.')
legend('Cross','Self','SKLD/4','2*S','Eq. (21)','SDR')
hold off

figure(22)
%loglog(SYMMnew(:,1),(factor/8)*(SYMMnew(:,2)+SYMMnew(:,3)),'r-'); hold on
%plot(SYMMnew(:,1),(factor/8)*(SYMMnew(:,4)+SYMMnew(:,5)),'b-');
loglog(SYMMnew(:,1),(factor/8)*SYMMnew(:,6),'ko');
hold on
plot(time,factor*shannon,'g+')
plot(time,mixtime./time,'k--')
plot(time,SDR_1(:,2),'bd')
plot(SDR_1(:,1),SDR_1(:,3),'rd')
plot(SDR_2(:,1),SDR_2(:,2),'pentagram')
plot(SDR_1(:,1),(2/sqrt(128*pi*D))*SDR_1(:,1).^-1.5,'b-.')
axis([1e-3 1e5 1e-6 tott])
legend('SKLD/4','2*S','Eq. (21)','SDR')
hold off


figure(5)
dtimes=0.5*(SDR_1(1:end-1,1)+SDR_1(2:end,1));
loglog(SDR_2(:,1),SDR_2(:,2),'-o'); hold on
axis([1e-3 1e5 1e-6 tott])
plot(SDR_1(:,1),SDR_1(:,2),'bd'); 
%loglog(dtimes,1*dtimes.^-1.5,'k--')
plot(dtimes,(1/sqrt(128*pi*D))*dtimes.^-1.5,'k--')
plot(dtimes,-diff(KLAB(:,2))./diff(KLAB(:,1)),'-+');
hold off
%loglog(SDR_2(:,1),SDR_2(:,2))


mix_JS=(1/(1-exp(-factor*log(2))))*(exp(-factor*shannon) - exp(-factor*log(2)));

figure(6)
semilogx(time,exp(-(factor/8)*SYMMnew(:,6)),'ok')
hold on
plot(time, mix_JS, 'r+')
plot(time,exp(-mixtime./time),'k--')

%plot(SYMMnew(:,1), exp(-SDR_1(:,2)),'bd')
%plot(SYMMnew(:,1),omega*react/k_f,'d')

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
axis([1e-3 1e4 1e-6 1e2])
%axis([0 2e2 0 2e-3])
hold off