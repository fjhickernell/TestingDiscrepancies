%% Testing discrepancies of sets of points with and without the first one
InitializeWorkspaceDisplay
tic
format short e
dvec = [2 5 10];
nd = length(dvec);
powvec = [0 1 2];
npow = length(powvec);
nvec = 2.^(4:12)';
nlen = length(nvec);
nstart = 2;
nrep = 499 + nstart;
D(nlen,nrep) = 0;
Dmiss(nlen,nrep) = 0;
kernelvec = {@starkernel,@centerkernel}; 
kernamevec = {'Star','Centered'};
nker = length(kernelvec);
summary(nlen,nrep,nker) = 0;
for mm = 1:npow
   pow = powvec(mm);
   weights = (1:dvec(nd)).^(-pow);
   for ll = 1:nd
      d = dvec(ll);
      for kk = 1:nker
         kernel = kernelvec{kk};
         kername = kernamevec{kk};
         for jj = 1:nrep
            p = sobolset(d);
            if jj > nstart
               p = scramble(p,'MatousekAffineOwen');
            end
            x = net(p,nvec(end)+1);
            for ii = 1:nlen
               n = nvec(ii);
               if jj == 2
                  D(ii,jj) = discrepancy(kernel,x(1:n,:)+1/(2*n),weights(1:d));
               else
                  D(ii,jj) = discrepancy(kernel,x(1:n,:),weights(1:d));
                  Dmiss(ii,jj) = discrepancy(kernel,x(2:n+1,:),weights(1:d));
               end
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
            h(3) = loglog([n n]/fudge, quantile(D(ii,3:end),[0.1 0.9]),'s-','color',MATLABBlue, ...
               'markersize',6,'markerfacecolor',MATLABBlue);
            h(1) = loglog(n/fudge,D(ii,1),'.','color',MATLABBlue,'markersize',20);
            h(2) = loglog(n/fudge,D(ii,2),'d','color',MATLABBlue, ...
               'markersize',6,'markerfacecolor',MATLABBlue);
            h(5) = loglog([n n]*fudge, quantile(Dmiss(ii,3:end),[0.1 0.9]),'s-','color',MATLABOrange, ...
               'markersize',6,'markerfacecolor',MATLABOrange);
            h(4) = loglog(n*fudge,Dmiss(ii,1),'.','color',MATLABOrange,'markersize',20);
         end
         % axis([10 1e4 1e-4 1])
         xlabel('\(n\)')
         ylabel([kername ' Discrepancy'])
         title(['\(d = ' int2str(d) '\), \quad weights \( = j^{-' int2str(pow) '}\)'])
         legend(h,{'1st no scr.', '1st no scr.\ ctr.', '1st 10\%-90\% scr.', ...
            '2nd no scr.', '2nd 10\%-90\% scr.'},'box','off','location','southwest')
         print('-depsc',[kername 'DiscrepancyPlot n=' int2str(n) 'd=' int2str(d) 'p=' int2str(pow) '.eps'])
      end
   end
end
toc
