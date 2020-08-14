function [gram,intker,intintker] = centerkernel(x,weights)
[n,d] = size(x);
intintker = prod(1 + weights/12);
intker = prod(1 + (0.5*weights).*(abs(x-1/2) - (x-1/2).^2),2);
x1 = reshape(x,[n 1 d]);
x2 = reshape(x,[1 n d]);
gram = prod(1 + (0.5*reshape(weights,[1 1 d])).*(abs(x1-1/2) + abs(x2-1/2) - abs(x1-x2)),3);

