function [D] = Dmatrix(vpart,al,at,Dmolvec)
% Right now 2-D only!

% Make Dxx=1,Dxy=2, Dyx=3, Dyy=4
eps=1e-30;
D=zeros(size(vpart,1),4);
vmag=sqrt( vpart(:,1).^2 + vpart(:,2).^2);

% Spreading part:
D(:,1)=(al-at).*(vpart(:,1).*vpart(:,1))./vmag;
D(:,2)=(al-at).*(vpart(:,1).*vpart(:,2))./vmag;
D(:,3)=D(:,2);
D(:,4)=(al-at).*(vpart(:,2).*vpart(:,2))./vmag;

D(vmag<eps,:)=0;

% Mixing part: (add if no mass transfer performed)
%D(:,1)=D(:,1)+at*vmag + Dmolvec;
%D(:,4)=D(:,4)+at*vmag + Dmolvec;

end

