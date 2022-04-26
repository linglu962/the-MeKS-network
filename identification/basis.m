function A=basis(x,y)
A=[ones(size(x,2),1) x' y' x'.^2 x'.*y' y'.^2 x'.^3  x'.^2.*y' x'.*y'.^2 y'.^3];
