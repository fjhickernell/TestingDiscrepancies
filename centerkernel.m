function [gram,intker,intintker] = centerkernel(x)
[n,d] = size(x);
intintker = (13/12)^d;
intker = prod(1 + 0.5*abs(x-1/2) - 0.5*(x-1/2).^2,2);
x1 = reshape(x,[n 1 d]);
x2 = reshape(x,[1 n d]);
gram = prod(1 + 0.5*(abs(x1-1/2) + abs(x2-1/2) - abs(x1-x2)),3);

