function [gram,intker,intintker] = starkernel(x,weights)
[n,d] = size(x);
intintker = prod(1 + weights/3);
intker = prod(1 + weights.*(1 - x.^2)/2,2);
gram = prod(1 + reshape(weights,[1 1 d]).*(1 - max(reshape(x,[n 1 d]),reshape(x,[1 n d]))),3);

