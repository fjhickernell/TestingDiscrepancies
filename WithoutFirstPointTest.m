%% Testing discrepancies of sets of points with and without the first one
gail.InitializeWorkspaceDisplay
tic
format short e
d = 2;
nvec = 2.^(4:12)';
nlen = length(nvec);
nrep = 500;
D(nlen,nrep) = 0;
Dmiss(nlen,nrep) = 0;
kernelvec = {@starkernel,@centerkernel}; 
kernamevec = {'Star','Centered'};
nker = length(kernelvec);
summary(nlen,nrep,nker) = 0;
for kk = 1:nker
   kernel = kernelvec{kk};
   kername = kernamevec{kk};
   for jj = 1:nrep
      p = sobolset(d);
      if jj > 1
         p = scramble(p,'MatousekAffineOwen');
      end
      x = net(p,nvec(end)+1);
      for ii = 1:nlen
         n = nvec(ii);
         D(ii,jj) = discrepancy(kernel,x(1:n,:));
         Dmiss(ii,jj) = discrepancy(kernel,x(2:n+1,:));
      end
   end
   summmary = ...
      [nvec D(:,1) Dmiss(:,2) mean(D(:,2:end),2) mean(Dmiss(:,2:end),2)];

   figure
   hold on
   set(gca,'XScale','log','YScale','log')
   fudge = 1.1;
   h(4,1) = 0;
   for ii = 1:nlen
      n = nvec(ii);
      h(2) = loglog([n n]/fudge, quantile(D(ii,2:end),[0.1 0.9]),'s-','color',MATLABBlue, ...
         'markersize',6,'markerfacecolor',MATLABBlue);
      h(1) = loglog(n/fudge,D(ii,1),'.','color',MATLABBlue,'markersize',20);
      h(4) = loglog([n n]*fudge, quantile(Dmiss(ii,2:end),[0.1 0.9]),'s-','color',MATLABOrange, ...
         'markersize',6,'markerfacecolor',MATLABOrange);
      h(3) = loglog(n*fudge,Dmiss(ii,1),'.','color',MATLABOrange,'markersize',20);
   end
   xlabel('\(n\)')
   ylabel([kername ' Discrepancy'])
   legend(h,{'1st no scr.', '1st 10\%-90\% scr.', ...
      '2nd no scr.', '2nd 10\%-90\% scr.'},'box','off')
   print('-depsc',[kername 'DiscrepancyPlot.eps'])
end


toc
