function [bchkrms,fixflag,sucflag]=DoMultSRC(zhat,ztrue,bhat,L,D,na,ns2,Kcoef2,beflt)
ncands=2;
if ns2>0
    k2=na-ns2+1;
   [zpar2,~] = ssearch(zhat(k2:end),L(k2:end,k2:end),D(k2:end),ncands);
   sucflag=sum(~(ztrue(k2:end)==zpar2(:,1)))==0;
   fixflag=1;
   sucflag=fixflag&sucflag;

   bchk=bhat-Kcoef2*(zhat(k2:end)-zpar2(:,1));
   bchkrms=norm(bchk(1:3));
else
   bchkrms=beflt;
   fixflag=0;
   sucflag=0;
end
    