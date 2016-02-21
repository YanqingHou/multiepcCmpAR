function [bchkrms,fixflag,sucflag]=DoMultFAR(zhat,ztrue,bhat,L,D,mu,Kfar,beflt)
ncands=2;
[zfixed,sqnorm] = ssearch(zhat,L,D,ncands);
sucflag=sum(~(ztrue==zfixed(:,1)))==0;
  
if sqnorm(1)/sqnorm(2)<=mu % FFRT passed
   bchk=bhat-Kfar*(zhat-zfixed(:,1));
   fixflag=1;
   bchkrms=norm(bchk(1:3));
else  
   bchkrms=beflt;
   fixflag=0;
end
sucflag=sucflag&fixflag;

    