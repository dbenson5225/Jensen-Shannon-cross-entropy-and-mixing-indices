clear all
gpu=gpuDeviceCount;      %  See if there's a gpu
gpu=0;                   %  Turn off gpu

% Running mean interface velocity was determined to be 0.0607
average_V=0.0607;
numcycles=1;
numphases=12;
load('V_Mays_test1.mat','K','h','vx','vy','L','H','I','por');
npreg=20000; npinter=1000; ndim=2, nspec=3;  np=npreg+npinter;
rng('shuffle');  % actually randomize
subs_min = min(np/100,10); subs_max=max(np/1000,100)

figures=true;
along = 0.0e-3; %longitudinal dispersivity
atrans= 0.0e-3; %transverse dispersivity
Dmol = 0.0e-4; %molecular diffusion (3.33e-6 m^2/hr)
Dmol = 2.0e-3*average_V; %molecular diffusion (3.33e-6 m^2/hr)
graphfmt='b-o'  
kappa = 1;  % If less than one, does RW

beta = 1;    % SPH beta
small=1e-200;

nnKx=size(K,1);  nnKy=size(K,2);  %nnKz=1
ndim=2; nspec=3;
for i=1:ndim
   xmin(i)=0; 
end
xmax(1) = L; xmax(2) = H;
ttot=(1/36)*por*L*L; dt=10; ntsteps=ttot/dt;
contamR=L/6/4; bufferR=L/6/4;  areaR=pi*contamR^2; areabuf=pi*((bufferR+contamR)^2-contamR^2);
% Define ICs %%%%%%%%%%%
A=zeros(np,ndim+nspec+2);   % x,(y,z),A,(B,C),mobile,alive ...
count=0
while count<(areaR/areabuf)*npreg
    xval=2*contamR*(rand-0.5); yval=2*contamR*(rand-0.5);
    if (xval^2+yval^2)<contamR^2
        count=count+1;
        A(count,1)=L/2+xval;
        A(count,2)=L/2+yval;
    end
end
count_contam=count;
while count<npreg
    xval=2*(contamR+bufferR)*(rand-0.5); yval=2*(contamR+bufferR)*(rand-0.5);
    if ( ((xval^2+yval^2)>contamR^2) & ((xval^2+yval^2)<(contamR+bufferR)^2) ) 
        count=count+1;
        A(count,1)=L/2+xval;
        A(count,2)=L/2+yval;
    end
end

%r=contamR*rand(npreg/2,1); theta=2*pi()*rand(npreg/2,1);
%A(1:npreg/2,1)=L/2+r.*cos(theta); A(1:npreg/2,2)=L/2+r.*sin(theta); 
A(1:count_contam,ndim+1)=1;
%r=contamR+bufferR*rand(npreg/2,1); theta=2*pi()*rand(npreg/2,1);
%A(npreg/2+1:npreg,1)=L/2+r.*cos(theta); A(npreg/2+1:npreg,2)=L/2+r.*sin(theta); 
A(count_contam+1:npreg,ndim+2)=1;

%%% put in the interface for tracking:
r=contamR;
theta=linspace(0,2*pi,npinter)';
A(npreg+1:end,1)=L/2+r.*cos(theta); A(npreg+1:end,2)=L/2+r.*sin(theta); 
% randomly assign A or B to interface
blah=rand(np-npreg,1);
A(npreg+1:end,ndim+1)=round(blah);   A(npreg+1:end,ndim+2)=1-round(blah);  A(npreg+1:end,ndim+3)=1;
%%%%
massA=sum(A(:,ndim+1)); massB=sum(A(:,ndim+2)); massC=sum(A(:,ndim+3));
A(:,ndim+1)=A(:,ndim+1)/massA;  A(:,ndim+2)=A(:,ndim+2)/massB;  A(:,ndim+3)=A(:,ndim+3)/massC;  

ds=pi*(contamR+bufferR)^2/size(A,1);
KLAB=zeros(100,3); SDR=zeros(100,3); JS=zeros(100,2); v_ave=zeros(100,1)

n=0; elapsed_total=0; v_run_mean=0.0;
for cycles=1:numcycles

%  Run one Mays cycle:
for phase=1:numphases
infile=['V_Mays_test',num2str(phase),'.mat'];
load(infile,'K','h','vx','vy','L','H','I','por');
vy=-vy;
%vx=0.01*ones(size(vx)); vy=0.005*ones(size(vy));  % constant velocity test

dKx(2)=(xmax(2)-xmin(2))/nnKy;  dKx(1)=(xmax(1)-xmin(1))/nnKx; 
vxave=(vx(:,1:nnKx)+vx(:,2:nnKx+1))./2;
vyave=(vy(1:nnKy,:)+vy(2:nnKy+1,:))./2;
vmag=sqrt(vxave.^2+vyave.^2);
%figure(10), imagesc((log(vmag))); colormap('gray'); title('ln(|v|)'); 
%axis equal; axis tight; hold on; colorbar;

elapsed_one_phase=0;

DMT = zeros(np,1);        

if(gpu>0) 
    A=gpuArray(A);
    xmin=gpuArray(xmin); xmax=gpuArray(xmax); DMT=gpuArray(DMT)
end

dt=0.1*(dKx(1)/max(max(vmag))); 
dxvec=dKx(1)*ones(np,1); dyvec=dKx(2)*ones(np,1);
ghost1=A(1:ndim); ghost2=A(1:ndim); ghost3=A(1:ndim); ghost4=A(1:ndim);  
% Dead=zeros(np,ndim+nspec);  ndead=0;
while elapsed_one_phase<ttot
    n=n+1;     
    if (figures)
      figure(111)
      mixcolor=zeros(size(A,1),3);
      mixcolor(:,1)=(np/2)*A(:,ndim+1); mixcolor(:,3)=(np/2)*A(:,ndim+2); mixcolor(:,2)=npinter*A(:,ndim+3);
      pointsize=5;
      scatter(A(:,1),L-A(:,2),pointsize,mixcolor,'filled');
      hold on
      scatter([L/3 2*L/3 L/2 L/2],[L/2 L/2 2*L/3 L/3],'k+')
%      axis([0 nnKx*dKx(1) 0 nnKy*dKx(2)]); axis('square')
      axis([L/3-5 2*L/3+5 L/3-5 2*L/3+5]); axis('square')
      %scatter(A(:,1),A(:,2),'o','MarkerFaceColor',mixcolor,'MarkerEdgeColor',mixcolor,'MarkerSize',1.5);
      drawnow;
      hold off

      figure(112)
      plot(A(npreg+2:end,1),L-A(npreg+2:end,2),'b-'); hold on
      plot(A(npreg+2:end,1),L-A(npreg+2:end,2),'g.','MarkerSize',1)
      axis([L/3-5 2*L/3+5 L/3-5 2*L/3+5]); 
      axis('square')
      hold off
      drawnow
    end        
%pause
    % Get the linear interps of real particle velocity
   [vpart]= vlin(A(:,1:ndim),vx,vy,dKx,xmin,xmax);

  %%%% Do mass-transfer first.  Figure out best domain decomposition:
   vmag=sqrt(vpart(:,1).^2+vpart(:,2).^2);
   DMT=Dmol+atrans*vmag;
   for i=1:ndim
     xminpart(i)=min(A(:,i));
     xmaxpart(i)=max(A(:,i));
   end

   [nxmt] = opt_DDC(A,ndim,subs_min,subs_max);
%    [A] = MT_matrix_fn(A,ndim,nspec,xmin,xmax,nxmt,DMT,dt);
   [A] = MT_matrix_fn(A,ndim,nspec,xminpart,xmaxpart,nxmt,DMT,dt,ds);

%   [A] = MT_matrix_fn_2(A,ndim,nspec,xminpart,xmaxpart,nxmt,DMT,dt,ds);

%%%%% only do random walks if necessary ...
   dDdx=zeros(np,1); dDdy=zeros(np,1);
   if kappa<1 
    
       % Make ghpst particles for random walk correction
        [ghost] = makeghost(np,A(:,1:ndim),dKx,xmin,xmax);  % Make ghost particles
        
     % Get bilinear interps of real particle velocity
        [vpartbi]=vbilin(A(:,1:ndim),vx,vy,dKx,xmin,xmax);
    
    % Get bilinear velocity interp of all ghost particles
        [ghostv]=vbilin(ghost,vx,vy,dKx,xmin,xmax);
    
    %  Adjust (numerically) for Dm at no-flow boundaries
    % [Dmolvecpart]=Dmadjust(Dmol,X,dKx,dKy,xmin,xmax,ymin,ymax,Dfactor);    
    
    % Dmolvecghost=Dmol*ones(size(ghost,1),1);
    % can comment this out when D adjustment is right...
    % [Dmolvecghost]=Dmadjust(Dmol,ghost,dx,dy,xmin,xmax,ymin,ymax,Dfactor);
    
    %  Make the various D matrices and their corrections ...
        [DRW] = Dmatrix(vpartbi,along,atrans,Dmol);
        [Dghost]=Dmatrix(ghostv,along,atrans,Dmol);

        dDdx=sum((Dghost(np+1:2*np,1:2) - Dghost(1:np,1:2))./ ...
             (ghost(np+1:2*np,1)-ghost(1:np,1)) ,2);
        dDdy=sum((Dghost(3*np+1:4*np,3:4) - Dghost(2*np+1:3*np,3:4))./ ...
            (ghost(3*np+1:4*np,2)-ghost(2*np+1:3*np,2)), 2);

        [B]=decomp(DRW);   % B*B = D_RW
        rando=randn(np,2);
        jumpx=sqrt(2*dt)*(B(:,1).*rando(:,1)+B(:,2).*rando(:,2));
        jumpy=sqrt(2*dt)*(B(:,3).*rando(:,1)+B(:,4).*rando(:,2));

        A(:,1) = A(:,1) + jumpx;                   % Random walk Particles
        A(:,2) = A(:,2) + jumpy;
   end

    A(:,1) = A(:,1) + dt*(vpart(:,1)+dDdx);    % Move particles by advection and pseudoadvection
    A(:,2) = A(:,2) + dt*(vpart(:,2)+dDdy); 
  
  % Keep track of average velocity of interface particles...
  v_ave(n)=mean(abs(vmag));
  v_run_mean=(elapsed_total*v_run_mean+dt*v_ave(n))/(elapsed_total+dt)
  
  % Increment total time and store mixing measures
  elapsed_one_phase=elapsed_one_phase+dt;
  elapsed_total=elapsed_total+dt;

  interface(n,1)=elapsed_total;
  interface(n,2)=sum(  sqrt(  sum(  (A(npreg+2:end,1:2)-A(npreg+1:end-1,1:2)).^2, 2))  );
  SDR(n,1)=elapsed_total;
  SDR(n,2)=(0.5*(npreg+npinter)/(areaR+areabuf))^2*sum(A(:,ndim+1).^2);
  KLAB(n,1)=elapsed_total;  KLAB(n,3)=cycles;
  apple=find((A(:,ndim+1)>small) & (A(:,ndim+2)>small));
  appleA=find(A(:,ndim+1)>small);  appleB=find(A(:,ndim+2)>small);

  Anorm=sum(A(apple,ndim+1)); Bnorm=sum(A(apple,ndim+2));
  KLAB(n,2)=(1/Anorm)*sum( A(apple,ndim+1).*log((Anorm/Bnorm)*A(apple,ndim+1)./A(apple,ndim+2)) ); 
  KLAB(n,3)=(1/Bnorm)*sum( A(apple,ndim+2).*log((Bnorm/Anorm)*A(apple,ndim+2)./A(apple,ndim+1)) );
  dt=min([ 0.5*( dKx(1)/max(max( vmag )))    0.5*dKx(1)*dKx(1)/(along*max(max(vpart))) ]);
  disp([' ']); disp(['Cyle ',num2str(cycles),';  Pump phase ',num2str(phase),';  Elapsed time = ', ...
             num2str(elapsed_total), ';  dt = ', num2str(dt)]); disp([' ']);
  m=0.5*(A(:,ndim+1)+A(:,ndim+2));
  JS(n,1) = elapsed_total;
  JS(n,2) = 0.5*sum( A(appleA,ndim+1).*log(A(appleA,ndim+1)./m(appleA)) ) + ...
            0.5*sum( A(appleB,ndim+2).*log(A(appleB,ndim+2)./m(appleB)) );

    %  Remove and keep track of particles that make it to the end
%        temp=A(A(:,1)>(nnKx-1.5)*dKx(1),:);
%        A(A(:,1)>(nnKx-1.5)*dKx(1),:)=[];
%        temp(:,1)=elapsed;   % replace x-coordinate with elapsed time
%        ndnow=size(temp,1);
%        Dead(ndead+1:ndead+ndnow,:)=temp;
%        ndead=ndead+ndnow;
       
  end  % Time loop (for one pumping phase
  if(figures)
      print("-f111", strcat('plumes\all_parts_RW',num2str(phase)), "-dpdf" );
      print("-f112", strcat('plumes\interface_RW',num2str(phase)), "-dpdf" );
  end


end  % phase loop (for Mays pumping scheme)

end  % multiple cycles loop

figure(1)
plot(interface(:,1),v_ave,graphfmt); hold on
xlabel('Time'); ylabel('Mean interface velocity magnitude')

figure(2)
%loglog(KLAB(:,1),KLAB(:,2),'-'); 
semilogy(interface(:,1),interface(:,2),graphfmt);
%loglog(KLAB(:,1),1./KLAB(:,1))
%axis([10 1e4 1e-6 1e0]);
xlabel('Time [T]');ylabel('Interface Length' )
hold on

dLdt=diff(interface(:,2))./diff(interface(:,1));
figure(22)
plot(0.5*(interface(1:end-1,1)+interface(2:end,1)),dLdt,graphfmt); hold on
xlabel('Time'); ylabel('Interface Growth Rate')

KLAB(n+1:end,:)=[]
dtimes=0.5*(KLAB(1:end-1,1)+KLAB(2:end,1));
figure(3)
loglog(dtimes,-diff(KLAB(:,2))./diff(KLAB(:,1)),'-'); hold on
loglog(dtimes,dtimes.^-2)
plot(interface(:,1),1e-6*interface(:,2),'-+');
axis([10 1e4 1e-12 1e-2]);
xlabel('Time [T]');ylabel('Mixing rate = -dI/dt [1/T]' )

SDR(n+1:end,:)=[]
dtimes=0.5*(SDR(1:end-1,1)+SDR(2:end,1));
figure(4)
semilogy(dtimes,-0.5*diff(SDR(:,2))./diff(SDR(:,1)),graphfmt); hold on
plot(dtimes,.02*dtimes.^-2,'-')
%plot(interface(:,1),1e-6*interface(:,2),'-+');
axis([10 1400 1e-8 1e-3]);
xlabel('Time [T]');ylabel('SDR' )

factor=2;
mix_JS=(1/(1-exp(-factor*log(2))))*(exp(-factor*JS(1:n,2)) - exp(-factor*log(2)));

figure(5)
plot(JS(1:n,1),mix_JS(1:n),graphfmt)
xlabel('Time'); ylabel('Mixing Index')
hold on 
%plot(KLAB(:,1), exp(-(factor/4)*KLAB(:,2)),'d' )
%plot(KLAB(:,1), exp(-(factor/4)*KLAB(:,3)),'sq' )
%plot(KLAB(:,1), exp(-(factor/8)*(KLAB(:,2)+KLAB(:,3))),'ko' )
axis([0 1400 0 0.4 ])


rate_JS=diff(mix_JS)./diff(JS(:,1));
rate_KL=diff(exp(-(factor/8)*(KLAB(:,2)+KLAB(:,3)))) ./ diff(KLAB(:,1));

figure(6)
semilogy(0.5*(JS(2:end,1)+JS(1:end-1,1)),rate_JS,graphfmt)
hold on 
%semilogy(dtimes,-0.5*diff(SDR(:,2))./diff(SDR(:,1)),'*');
%axis([0 elapsed_total 0 max(rate_JS)])
%plot(KLAB(:,1), exp(-(factor/4)*KLAB(:,2)) )
%plot(KLAB(:,1), exp(-(factor/4)*KLAB(:,3)) )
%plot(0.5*(KLAB(2:end,1)+KLAB(1:end-1,1)), rate_KL,'ko' )
axis([0 1400 1e-8 1e-3 ])

figure(66)
plot(0.5*(JS(2:end,1)+JS(1:end-1,1)),rate_JS,graphfmt)
hold on 
xlabel('Time'); ylabel('Mixing Rate')
%axis([0 elapsed_total 0 max(rate_JS)])
%plot(KLAB(:,1), exp(-(factor/4)*KLAB(:,2)) )
%plot(KLAB(:,1), exp(-(factor/4)*KLAB(:,3)) )
%plot(0.5*(KLAB(2:end,1)+KLAB(1:end-1,1)), rate_KL,'ko' )
axis([0 1400 0 6e-4 ])


% construct a BTC
%Dead(Dead(:,1)<.1e-10,:)=[];
%prob=Nstart-1:-1:size(X,1);
%prob=prob'/Nstart;
%figure(5)
%loglog(Dead(:,1),prob,'o')
%loglog(Dead(:,1),1-cumsum(Dead(:,3)),'o')
%hold on


%    figure(2)
%    histogram(X(:,2),10)
%    title('Position Histogram')
    %drawnow
%plot(A,B);
%grid on;

