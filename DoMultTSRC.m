function [bchkrms,fixflag,sucflag,ns]=DoMultTSRC(zhat,ztrue,bhat,L,D,na,ns1,ns2,mu1,mu2,Kcoef1,Kcoef2,beflt)
ncands=2;
        if ns1>0
            k1=na-ns1+1;
            [zpar1,sqnorm1] = ssearch(zhat(k1:end),L(k1:end,k1:end),D(k1:end),ncands);
            sucflag1=sum(~(ztrue(k1:end)==zpar1(:,1)))==0;
            if sqnorm1(1)/sqnorm1(2)>mu1
                if ns2>0
                    k2=na-ns2+1;
                    [zpar2,sqnorm2] = ssearch(zhat(k2:end),L(k2:end,k2:end),D(k2:end),ncands);
                    sucflag2=sum(~(ztrue(k2:end)==zpar2(:,1)))==0;
                    %                     ratio test for the second step.
                    if sqnorm2(1)/sqnorm2(2)>mu2
                        %                             use the float solution, update the
                        %                             availability
                        fixflag=0; sucflag=0; bchkrms=beflt;ns=0;
                    else
                        fixflag=1;
                        sucflag=fixflag&sucflag2;
                        ns=ns2;
%                         failnum2=fixnum2-sucnum2;
                        bchk2=bhat-Kcoef2*(zhat(k2:end)-zpar2(:,1));
%                         availnum=availnum+double(norm(bchk2(1:3))<=errs);
                        bchkrms=norm(bchk2(1:3));
                        %                             update the baseline solution with ns2(ipf)
                    end
                else % ns2==0, use float solution.
                    fixflag=0; sucflag=0; bchkrms=beflt;ns=0;
                end
            else % subset in the first step accepted

                fixflag=1;
                sucflag=sucflag1&fixflag;
                ns=ns1;
                bchk1=bhat-Kcoef1*(zhat(k1:end)-zpar1(:,1));
                bchkrms=norm(bchk1(1:3));
                %                       use the optimal subset
                %                       collect the fix num, success num and fail num
                %                       calculate the availability
            end
            
        else % optsub(ipf)==0
            %                   use the float solution
            fixflag=0; sucflag=0; bchkrms=beflt;ns=0;
        end