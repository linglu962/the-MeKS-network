function A=fun(~,x)
global ak bk bs k0 k1 n p;
A(1)=10*(ak+bk*(x(1)/10).^n./(k0^n+(x(1)/10).^n)-x(1)/10./(1+x(1)/10+x(2)/2));
A(2)=2*(bs./(1+(x(1)/10/k1).^p)-x(2)/2./(1+x(1)/10+x(2)/2));
end