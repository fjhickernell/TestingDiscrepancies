function [gram,intker,intintker] = starkernel(x)
[n,d] = size(x);
intintker = (4/3)^d;
intker = prod((3 - x.^2)/2,2);
gram = prod(2 - max(reshape(x,[n 1 d]),reshape(x,[1 n d])),3);

