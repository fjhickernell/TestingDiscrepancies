function D = discrepancy(kernel,x,weights)
n = size(x,1);
[gram,intker,intintker] = kernel(x,weights);
D2  = intintker - (2/n)*sum(intker) + (1/n^2)*sum(sum(gram));
D = sqrt(D2);
