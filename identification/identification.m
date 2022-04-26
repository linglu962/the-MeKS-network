clear;clc;
% mesh parameters
xmin=0;xmax=3;Mx=2e4;
ymin=2;ymax=7;My=2e4;
M=Mx*My;
x0=linspace(xmin,xmax,Mx);
y0=linspace(ymin,ymax,My);
[X0,Y0]=meshgrid(x0,y0);
Y0=flipud(Y0);
delta_t=0.001;
% system parameters
ak=0.004;bk=0.14;bs=0.68;k0=0.2;k1=0.222;n=2;p=5;
fx=10*(ak+bk*(X0/10).^n./(k0^n+(X0/10).^n)-X0/10./(1+X0/10+Y0/2));
fy=2*(bs./(1+(X0/10/k1).^p)-Y0/2./(1+X0/10+Y0/2));
alphax=1.8;alphay=1.8;betax=0;betay=0;epsilonk=0.25;epsilons=0.25;d=0;


% generating data and necessary parameters in learning
X=X0+fx*delta_t+delta_t^(1/2)*sqrt(d)*randn(Mx,My)+epsilonk*delta_t^(1/alphax)*stblrnd(alphax,betax,1,0,Mx,My);
Y=Y0+fy*delta_t+delta_t^(1/2)*sqrt(d)*randn(Mx,My)+epsilons*delta_t^(1/alphay)*stblrnd(alphay,betay,1,0,Mx,My);
n_basis_function=10;
Qx=50;Qy=50;
Q=Qx*Qy;% Number of bins

% coordinate transformation
xxmin=-1;xxmax=1;
yymin=-1;yymax=1;
xx0=linspace(xxmin,xxmax,Mx);
yy0=linspace(yymin,yymax,My);
[XX0,YY0]=meshgrid(xx0,yy0);
YY0=flipud(YY0);
XX=(xxmax-xxmin)/(xmax-xmin)*X+(xmax*xxmin-xmin*xxmax)/(xmax-xmin);
YY=(yymax-yymin)/(ymax-ymin)*Y+(ymax*yymin-ymin*yymax)/(ymax-ymin);
Rx=XX-XX0;
Ry=YY-YY0;
% identification of stable parameter ��,levy noise intensity ��
% x-direction
N=2;e=1;m=5;
n_pos=zeros(1,N+1);n_neg=zeros(1,N+1);
for i=1:N+1
    n_pos(i)=sum((sum(Rx>=e*m^(i-1)&Rx<e*m^i)));
    n_neg(i)=sum((sum(Rx>=-e*m^i&Rx<-e*m^(i-1))));
end
k=1:N;
alphax1=(k*log(m)).^(-1).*log((n_pos(1)+n_neg(1))./(n_pos(2:end)+n_neg(2:end)));
alphax_test=sum(alphax1)/N;

rho=sum(n_neg)/sum(n_pos);
betax=(1-rho)/(1+rho);

k=0:N;
if alphax_test<0.98||alphax_test>1.02
    k_alphax=alphax_test*(1-alphax_test)/gamma(2-alphax_test)/cos(pi*alphax_test/2);
else
    k_alphax=2/pi;
end
epsilonk1=(alphax_test*e^alphax_test*m.^(k*alphax_test).*(n_pos+n_neg)/k_alphax/delta_t/Mx/My/(1-m^(-alphax_test))).^(1/alphax_test);
epsilonk_test_p=sum(epsilonk1)/(N+1);
epsilonk_test=(xmax-xmin)/(xxmax-xxmin)*epsilonk_test_p;

% y-direction
n_pos=zeros(1,N+1);n_neg=zeros(1,N+1);
for i=1:N+1
    n_pos(i)=sum((sum(Ry>=e*m^(i-1)&Ry<e*m^i)));
    n_neg(i)=sum((sum(Ry>=-e*m^i&Ry<-e*m^(i-1))));
end
k=1:N;
alphay1=(k*log(m)).^(-1).*log((n_pos(1)+n_neg(1))./(n_pos(2:end)+n_neg(2:end)));
alphay_test=sum(alphay1)/N;

rho=sum(n_neg)/sum(n_pos);
betay=(1-rho)/(1+rho);

k=0:N;
if alphay_test<0.98||alphay_test>1.02
    k_alphay=alphay_test*(1-alphay_test)/gamma(2-alphay_test)/cos(pi*alphay_test/2);
else
    k_alphay=2/pi;
end
epsilons1=(alphay_test*e^alphay_test*m.^(k*alphay_test).*(n_pos+n_neg)/k_alphay/delta_t/Mx/My/(1-m^(-alphay_test))).^(1/alphay_test);
epsilons_test_p=sum(epsilons1)/(N+1);
epsilons_test=(ymax-ymin)/(yymax-yymin)*epsilons_test_p;

% binning
p=find(abs(Rx)<e&abs(Ry)<e);
XX01=XX0(p);XX1=XX(p);Rx1=XX1-XX01;
YY01=YY0(p);YY1=YY(p);Ry1=YY1-YY01;
M1=size(p,1);
x1Q=zeros(1,Q);y1Q=zeros(1,Q);Rx1Q=zeros(1,Q);Ry1Q=zeros(1,Q);Weight=zeros(1,Q);RR1Q=zeros(1,Q);
hx=(xxmax-xxmin)/Qx;hy=(yymax-yymin)/Qy;
for i=1:Qx
    for j=1:Qy
        aa=find(XX01>=xxmin+(i-1)*hx&XX01<xxmin+i*hx&YY01>=yymin+(j-1)*hy&YY01<yymin+j*hy);
        aa=aa';
        x1Q(Qy*(i-1)+j)=sum(XX01(aa))/size(aa,2);
        y1Q(Qy*(i-1)+j)=sum(YY01(aa))/size(aa,2);
        Rx1Q(Qy*(i-1)+j)=sum(Rx1(aa))/size(aa,2);
        Ry1Q(Qy*(i-1)+j)=sum(Ry1(aa))/size(aa,2);
        Weight(Qy*(i-1)+j)=size(aa,2)/M1;
    end
end

A=diag(Weight)*basis(x1Q,y1Q);
if alphax_test<0.98||alphax_test>1.02
    Rx_alpha_beta=k_alphax*betax*e^(1-alphax)/epsilonk_test^(-alphax_test)/(1-alphax_test);
else
    Rx_alpha_beta=epsilonk_test*2/pi*betax*log(e);
end
if alphay_test<0.98||alphay_test>1.02
    Ry_alpha_beta=k_alphay*betay*e^(1-alphay)/epsilons_test^(-alphay_test)/(1-alphay_test);
else
    Ry_alpha_beta=epsilons_test*2/pi*betay*log(e);
end
Bx=diag(Weight)*(M1/M/delta_t*Rx1Q'-Rx_alpha_beta);
By=diag(Weight)*(M1/M/delta_t*Ry1Q'-Ry_alpha_beta);

% % identification of the drift terms
% x-direction
X00=A;Y00=Bx;
% initial solution
c0=(X00'*X00)\(X00'*Y00);
X=X00;Y=Y00;

% intermediate solution and record matrix
c=c0;
cc=zeros(n_basis_function,n_basis_function);
cc(:,1)=c0;

% pointer vector
p=1:n_basis_function;

% sparsity enforcement
for i=1:n_basis_function-1
    [MM,I]=min(abs(c));
    X(:,I)=[];
    c=(X'*X)\(X'*Y);
    p(abs(c0)==MM)=0;
    c0(abs(c0)==MM)=0; 
    c0(find(p))=c;  %#ok<FNDSB>
    cc(:,i+1)=c0;
end

% % Cross Validation
% run times and numbers of folds
Nk=100;fold=5;
% vector of CV scores
delta=zeros(1,n_basis_function);
%  calculation of CV scores
for i=1:Nk
    %  data grouping
    s = randperm(Q);
    ngroup=zeros(1,fold);
    ngroup(1) = ceil(rand*(Q-fold+1));
    for j=2:fold-1
        ngroup(j) = ceil(rand*(Q-sum(ngroup(1:j-1))-fold+j));
    end
    ngroup(end)=Q-sum(ngroup(1:fold-1));
    for j=1:n_basis_function
        part=zeros(1,fold);
        p1=s(1:ngroup(1));
        p2=s;p2(1:ngroup(1))=[];
        part(1)=sum((Y00(p1)-X00(p1,:)*SSR(X00(p2,:),Y00(p2),j-1)).^2);
        for m=2:fold
            p1=s(sum(ngroup(1:m-1))+1:sum(ngroup(1:m)));
            p2=s;p2(sum(ngroup(1:m-1))+1:sum(ngroup(1:m)))=[];
            part(m)=sum((Y00(p1)-X00(p1,:)*SSR(X00(p2,:),Y00(p2),j-1)).^2);
        end
        delta0=1/fold*sum(part);
        delta(end-j+1)=delta(end-j+1)+delta0;
    end
end
delta=sqrt(delta/Nk);
ratio=delta(1:n_basis_function-1)./delta(2:n_basis_function);

% drawing of CV scores and ratios
figure
plot(1:n_basis_function-1,delta(1:end-1));
hold off
figure
plot(2:n_basis_function,ratio(1:end));
hold off

% y-direction
X00=A;Y00=By;
% initial solution
c0=(X00'*X00)\(X00'*Y00);
X=X00;Y=Y00;

% intermediate solution and record matrix
c=c0;
cc=zeros(n_basis_function,n_basis_function);
cc(:,1)=c0;

% pointer vector
p=1:n_basis_function;

% sparsity enforcement
for i=1:n_basis_function-1
    [MM,I]=min(abs(c));
    X(:,I)=[];
    c=(X'*X)\(X'*Y);
    p(abs(c0)==MM)=0;
    c0(abs(c0)==MM)=0; 
    c0(find(p))=c;  %#ok<FNDSB>
    cc(:,i+1)=c0;
end

% % Cross Validation
% run times and numbers of folds
Nk=100;fold=5;
% vector of CV scores
delta=zeros(1,n_basis_function);
%  calculation of CV scores
for i=1:Nk
    %  data grouping
    s = randperm(Q);
    ngroup=zeros(1,fold);
    ngroup(1) = ceil(rand*(Q-fold+1));
    for j=2:fold-1
        ngroup(j) = ceil(rand*(Q-sum(ngroup(1:j-1))-fold+j));
    end
    ngroup(end)=Q-sum(ngroup(1:fold-1));
    for j=1:n_basis_function
        part=zeros(1,fold);
        p1=s(1:ngroup(1));
        p2=s;p2(1:ngroup(1))=[];
        part(1)=sum((Y00(p1)-X00(p1,:)*SSR(X00(p2,:),Y00(p2),j-1)).^2);
        for m=2:fold
            p1=s(sum(ngroup(1:m-1))+1:sum(ngroup(1:m)));
            p2=s;p2(sum(ngroup(1:m-1))+1:sum(ngroup(1:m)))=[];
            part(m)=sum((Y00(p1)-X00(p1,:)*SSR(X00(p2,:),Y00(p2),j-1)).^2);
        end
        delta0=1/fold*sum(part);
        delta(end-j+1)=delta(end-j+1)+delta0;
    end
end
delta=sqrt(delta/Nk);
ratio=delta(1:n_basis_function-1)./delta(2:n_basis_function);

% drawing of CV scores and ratios
figure
plot(1:n_basis_function-1,delta(1:end-1));
hold off
figure
plot(2:n_basis_function,ratio(1:end));
hold off