clear;clc;
global ak bk bs k0 k1 n p;
% syms X Y ak bk bs k0 k1 p n;
ak=0.004;bk=0.14;bs=0.68;k0=0.2;k1=0.222;n=2;p=5;
X=0.8568;Y=4.4938;
fx=10*(ak+bk*(X/10).^n./(k0^n+(X/10).^n)-X/10./(1+X/10+Y/2));
fy=2*(bs./(1+(X/10/k1).^p)-Y/2./(1+X/10+Y/2));
% fx1=diff(fx,X);fx2=diff(fx,Y);
% fy1=diff(fy,X);fy2=diff(fy,Y);
% fx1 = X/(10*(X/10 + Y/2 + 1)^2) - 1/(X/10 + Y/2 + 1) + (bk*n*(X/10)^(n - 1))/(k0^n + (X/10)^n) - (bk*n*(X/10)^n*(X/10)^(n - 1))/(k0^n + (X/10)^n)^2;
% fx2 = X/(2*(X/10 + Y/2 + 1)^2);
% fy1 = Y/(10*(X/10 + Y/2 + 1)^2) - (bs*p*(X/(10*k1))^(p - 1))/(5*k1*((X/(10*k1))^p + 1)^2);
% fy2 = Y/(2*(X/10 + Y/2 + 1)^2) - 1/(X/10 + Y/2 + 1);
fx1 = 0.1358;fx2 = 0.0386; 
fy1 = -0.0263;fy2 = -0.0978;
A=[fx1 fx2;fy1 fy2];
[XX,B]=eig(A);
saddle=[X Y];
delta_x=0.01;a=XX(1,2);b=XX(2,2);
x_left=[X-delta_x Y-delta_x*b/a];
x_right=[X+delta_x Y+delta_x*b/a];
x=linspace(x_left(1),x_right(1),10);
y=linspace(x_left(2),x_right(2),10);
stable_left=zeros(2,3300,5);stable_right=zeros(2,3300,5);
stable_left(1,1,:)=x(1:5);stable_right(1,1,:)=x(6:10);
stable_left(2,1,:)=y(1:5);stable_right(2,1,:)=y(6:10);
for i=1:5
    for j=1:3299
    stable_left(:,j+1,i)=rk4(0,0.01,stable_left(:,j,i)');
    stable_right(:,j+1,i)=rk4(0,0.01,stable_right(:,j,i)');
    end
end
figure
i=5;
    plot(stable_left(1,:,i),stable_left(2,:,i));
    hold on
    plot(stable_right(1,:,i),stable_right(2,:,i));
    hold off
    stable_manifold=[stable_right(1,:,5);stable_right(2,:,5)];
    
    