clear;clc;
load('1.85ks.mat');
load('stable_manifold.mat');
d=zeros(1,size(k,2));
t=linspace(0,100,size(k,2));
for i=1:size(k,2)
    dd=zeros(1,3300);
    for j=1:3300
        dd(j)=sqrt((k(i)-stable_manifold(1,j))^2+(s(i)-stable_manifold(2,j))^2);
    end
    d(i)=min(dd);
end
[M,I]=min(d);
% tipping_time=zeros(1,37);
load('tipping time.mat');
% tipping_time(5)=t(I);
% save('tipping time', 'tipping_time');
figure
plot(k,s);
hold on
plot(stable_manifold(1,:),stable_manifold(2,:));
hold on
plot(k(I),s(I),'*');
hold off


figure
plot(0.15262,4.3148,'rp');
hold on
plot(1.5732,3.1562,'rp');
hold on
plot(0.8568,4.4938,'y*');
hold on
plot(stable_manifold(1,:),stable_manifold(2,:));
hold off

alpha=0.1:0.05:1.9;
figure
plot(alpha,tipping_time,'rp-');
hold off

figure
plot(alpha(8:end),tipping_time(8:end),'rp-');
hold off