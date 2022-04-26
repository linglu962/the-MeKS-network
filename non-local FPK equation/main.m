clear;
clc;

%%% parameters
xmin=0;
xmax=3;
ymin=2;
ymax=7;
J=100;
v=linspace(-2,2,4*J+1);
v(1)=[];
v(end)=[];
hx=v(2)-v(1);
w=linspace(-2,2,4*J+1);
w(1)=[];
w(end)=[];
hy=w(2)-w(1);
deltat=1e-2;
T=100;
nT=T/deltat;
t=linspace(0,100,nT);
ak=0.004;bk=0.14;bs=0.68;k0=0.2;k1=0.222;n=2;p=5;
alpha=0.1;
% sigma=0.001;
% sigmax=(epsilong*sigma)^(1/alpha);
% sigmay=sigma^(1/alpha);
sigmax=0.25;
sigmay=0.25;
Calpha=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*sqrt(pi)*gamma(1-alpha/2));
sigmaB=0.0001;
Chx=sigmaB/2*(2/(xmax-xmin))^2-(2*sigmax/(xmax-xmin))^alpha*Calpha*zeta(alpha-1)*hx^(2-alpha);
Chy=sigmaB/2*(2/(xmax-xmin))^2-(2*sigmay/(ymax-ymin))^alpha*Calpha*zeta(alpha-1)*hy^(2-alpha);
C1x=Chx/hx^2;
C2x=Calpha/alpha*(2*sigmax/(xmax-xmin))^alpha;
C3x=Calpha*(2*sigmax/(xmax-xmin))^alpha*hx;
C1y=Chy/hy^2;
C2y=Calpha/alpha*(2*sigmay/(ymax-ymin))^alpha;
C3y=Calpha*(2*sigmay/(ymax-ymin))^alpha*hy;

%%% initialization
x=(xmax-xmin)/2*v(J+1:3*J-1)+(xmax+xmin)/2;
y=(ymax-ymin)/2*w(J+1:3*J-1)+(ymax+ymin)/2;
[X,Y]=meshgrid(x,y);
X=X';
Y=Y';
P0=40/pi*exp(-40*((X-0.15262).^2+(Y-4.3148).^2));
P=reshape(P0,(2*J-1)*(2*J-1),1);
PP0=40/pi*exp(-40*((X-1.5732).^2+(Y-3.1562).^2));
PP=PP0;


A1x=zeros(2*J-1,2*J-1);
for i=1:2*J-1
    A1x(i,:)=1./abs(v(2*J+1-i:4*J-1-i)).^(1+alpha);
    A1x(i,i)=0;
end
v1=0.5./abs(v(1:2*J-1)).^(1+alpha);
v2=0.5./abs(v(2*J+1:4*J-1)).^(1+alpha);
v1=flipud(v1);
v2=flipud(v2);
A0=[v1' A1x v2'];
for i=1:2*J-1
    A1x(i,i)=-sum(A0(i,:));
end
A1x=C3x*A1x;
v0=C2x./(1+v(J+1:3*J-1)).^alpha+C2x./(1-v(J+1:3*J-1)).^alpha;
A1x=A1x-2*C1x*eye(2*J-1)-diag(v0);
A1x(1:2*J-2,2:2*J-1)=A1x(1:2*J-2,2:2*J-1)+C1x*eye(2*J-2);
A1x(2:2*J-1,1:2*J-2)=A1x(2:2*J-1,1:2*J-2)+C1x*eye(2*J-2);

A1y=zeros(2*J-1,2*J-1);
for i=1:2*J-1
    A1y(i,:)=1./abs(w(2*J+1-i:4*J-1-i)).^(1+alpha);
    A1y(i,i)=0;
end
w1=0.5./abs(w(1:2*J-1)).^(1+alpha);
w2=0.5./abs(w(2*J+1:4*J-1)).^(1+alpha);
w1=flipud(w1);
w2=flipud(w2);
A0=[w1' A1y w2'];
for i=1:2*J-1
    A1y(i,i)=-sum(A0(i,:));
end
A1y=C3y*A1y;
w0=C2y./(1+w(J+1:3*J-1)).^alpha+C2y./(1-w(J+1:3*J-1)).^alpha;
A1y=A1y-2*C1y*eye(2*J-1)-diag(w0);
A1y(1:2*J-2,2:2*J-1)=A1y(1:2*J-2,2:2*J-1)+C1y*eye(2*J-2);
A1y(2:2*J-1,1:2*J-2)=A1y(2:2*J-1,1:2*J-2)+C1y*eye(2*J-2);

% L1=1:(2*J-1);
% L2=0:(2*J-1):(2*J-2)*(2*J-1);
% A1=sparse((2*J-1)*(2*J-1),(2*J-1)*(2*J-1));
% for i=1:2*J-1
%     A1((i-1)*(2*J-1)+L1,(i-1)*(2*J-1)+L1)=A1((i-1)*(2*J-1)+L1,(i-1)*(2*J-1)+L1)+A1x;
%     A1(L2+i,L2+i)=A1(L2+i,L2+i)+A1y;
% end

%%% drift
Kx=0.3;
Ky=20;
fx=10*(ak+bk*(X/10).^n./(k0^n+(X/10).^n)-X/10./(1+X/10+Y/2));
Mfx1=(fx+max(max(abs(fx))))/2;
Mfx2=(fx-max(max(abs(fx))))/2;
fy=2*(bs./(1+(X/10/k1).^p)-Y/2./(1+X/10+Y/2));
Mfy1=(fy+max(max(abs(fy))))/2;
Mfy2=(fy-max(max(abs(fy))))/2;

L=nT/10;
Ptotal=[];

%%% integral
PA=zeros(2*J-1,2*J-1,nT);PB=zeros(2*J-1,2*J-1,nT);
for i=1:nT
    PA(:,:,i)=reshape(P,2*J-1,2*J-1);
    U=P;
    
    Ux=reshape(U,2*J-1,2*J-1);
    phi=Mfx1.*Ux;
    phixp=2/(xmax-xmin)*leftbiased(phi,hx);
    phixp=reshape(phixp,(2*J-1)*(2*J-1),1);
    phi=Mfx2.*Ux;
    phixn=2/(xmax-xmin)*rightbiased(phi,hx);
    phixn=reshape(phixn,(2*J-1)*(2*J-1),1);
    
    phi=Mfy1.*Ux;
    phi=phi';
    phiyp=2/(ymax-ymin)*leftbiased(phi,hy);
    phiyp=reshape(phiyp',(2*J-1)*(2*J-1),1);
    phi=Mfy2.*Ux;
    phi=phi';
    phiyn=2/(ymax-ymin)*rightbiased(phi,hy);
    phiyn=reshape(phiyn',(2*J-1)*(2*J-1),1);
    U01=A1x*Ux+(A1y*Ux')';
    U01=reshape(U01,(2*J-1)*(2*J-1),1);
    U1=U+deltat*(U01-phixp-phixn-phiyp-phiyn);
    
    
    Ux=reshape(U1,2*J-1,2*J-1);
    phi=Mfx1.*Ux;
    phixp=2/(xmax-xmin)*leftbiased(phi,hx);
    phixp=reshape(phixp,(2*J-1)*(2*J-1),1);
    phi=Mfx2.*Ux;
    phixn=2/(xmax-xmin)*rightbiased(phi,hx);
    phixn=reshape(phixn,(2*J-1)*(2*J-1),1);
    
    phi=Mfy1.*Ux;
    phi=phi';
    phiyp=2/(ymax-ymin)*leftbiased(phi,hy);
    phiyp=reshape(phiyp',(2*J-1)*(2*J-1),1);
    phi=Mfy2.*Ux;
    phi=phi';
    phiyn=2/(ymax-ymin)*rightbiased(phi,hy);
    phiyn=reshape(phiyn',(2*J-1)*(2*J-1),1);
    U02=A1x*Ux+(A1y*Ux')';
    U02=reshape(U02,(2*J-1)*(2*J-1),1);
    U2=3/4*U+1/4*U1+deltat/4*(U02-phixp-phixn-phiyp-phiyn);
    
    
    Ux=reshape(U2,2*J-1,2*J-1);
    phi=Mfx1.*Ux;
    phixp=2/(xmax-xmin)*leftbiased(phi,hx);
    phixp=reshape(phixp,(2*J-1)*(2*J-1),1);
    phi=Mfx2.*Ux;
    phixn=2/(xmax-xmin)*rightbiased(phi,hx);
    phixn=reshape(phixn,(2*J-1)*(2*J-1),1);
    
    phi=Mfy1.*Ux;
    phi=phi';
    phiyp=2/(ymax-ymin)*leftbiased(phi,hy);
    phiyp=reshape(phiyp',(2*J-1)*(2*J-1),1);
    phi=Mfy2.*Ux;
    phi=phi';
    phiyn=2/(ymax-ymin)*rightbiased(phi,hy);
    phiyn=reshape(phiyn',(2*J-1)*(2*J-1),1);
    U03=A1x*Ux+(A1y*Ux')';
    U03=reshape(U03,(2*J-1)*(2*J-1),1);
    P=1/3*U+2/3*U2+deltat*2/3*(U03-phixp-phixn-phiyp-phiyn);
    
    if rem(i,L)==0
        Ptotal=[Ptotal P];
    end
    
    PB(:,:,nT-i+1)=PP;
    Pvs1=[zeros(1,2*J-1);PP];Pvs1(end,:)=[];
    Pvp1=[PP;zeros(1,2*J-1)];Pvp1(1,:)=[];
    Pws1=[zeros(2*J-1,1) PP];Pws1(:,end)=[];
    Pwp1=[PP zeros(2*J-1,1)];Pwp1(:,1)=[];
    P1=(Pvp1-Pvs1)/hx/2.*fx*2/(xmax-xmin)+(Pwp1-Pws1)/hy/2.*fy*2/(ymax-ymin);
    U01=A1x*PP+(A1y*PP')';
    U1=PP+deltat*(U01+P1);
    
    Pvs1=[zeros(1,2*J-1);U1];Pvs1(end,:)=[];
    Pvp1=[U1;zeros(1,2*J-1)];Pvp1(1,:)=[];
    Pws1=[zeros(2*J-1,1) U1];Pws1(:,end)=[];
    Pwp1=[U1 zeros(2*J-1,1)];Pwp1(:,1)=[];
    P1=(Pvp1-Pvs1)/hx/2.*fx*2/(xmax-xmin)+(Pwp1-Pws1)/hy/2.*fy*2/(ymax-ymin);
    U02=A1x*U1+(A1y*U1')';
    U2=3/4*PP+1/4*U1+deltat/4*(U02+P1);
    
    Pvs1=[zeros(1,2*J-1);U2];Pvs1(end,:)=[];
    Pvp1=[U2;zeros(1,2*J-1)];Pvp1(1,:)=[];
    Pws1=[zeros(2*J-1,1) U2];Pws1(:,end)=[];
    Pwp1=[U2 zeros(2*J-1,1)];Pwp1(:,1)=[];
    P1=(Pvp1-Pvs1)/hx/2.*fx*2/(xmax-xmin)+(Pwp1-Pws1)/hy/2.*fy*2/(ymax-ymin);
    U03=A1x*U2+(A1y*U2')';
    PP=1/3*PP+2/3*U2+2/3*deltat*(U03+P1);
end

k=zeros(1,nT);s=zeros(1,nT);
PPP=PA.*PB;
for i=1:nT
    p=PPP(:,:,i);
    [M,I]=max(p(:));
    k(i)=X(I);s(i)=Y(I);
end
figure;
plot(k,s);

% figure;
% plot(t,k);
% figure;
% plot(t,s);


% for i=1:10
%     P3=reshape(Ptotal(:,i),2*J-1,2*J-1);
%     figure;
%     mesh(X,Y,P3);
% end

