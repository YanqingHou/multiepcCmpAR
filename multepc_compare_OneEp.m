function [resallARmethods]=multepc_compare_OneEp(ep0,res)
%% rewrite to multi epoch case
Pfttls=[0.001,0.01]';
Pfttlslen=length(Pfttls);
Pf1s=0.9*Pfttls;


resallARmethods=[];
noepc=size(res,2);
farbchkrmss=nan(Pfttlslen,noepc); srcbchkrmss=farbchkrmss; tsrcbchkrmss=farbchkrmss;
farfixflags=farbchkrmss; srcfixflags=farbchkrmss; tsrcfixflags=farbchkrmss;
farsucflags=farbchkrmss; srcsucflags=farbchkrmss; tsrcsucflags=farbchkrmss;tsrcfixns=zeros(Pfttlslen,noepc);
fltrmss=farbchkrmss;
for epc=1:noepc
    Qa=res(ep0,epc).Qa;
    Qb=res(ep0,epc).Qb;
    Qab=res(ep0,epc).Qab;
    %     Psb=res(ep0,epc).Psb;
    %% 0.Initialize the float ambiguities%%%%%%%%%%%%%%%%%%%%%%%%%%%
    na      = size(Qa,1);
    nb      = size(Qb,1);
    Qx=[Qa,Qab;
        Qab', Qb];
    Qx      = (tril(Qx,0)+tril(Qx,-1)');
    
    atrue=randi([-200 200],na,1);% generate random integers
    btrue=zeros(nb,1);
    xtrue=[atrue;btrue];
    
    [Qzhat,Z,L,D,ztrue,~] = decorrel(Qa,atrue);
    Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');
    
    Qbz=Qab'*Z;
    Kfar=Qbz/Qzhat;
    
    Psb = prod ( 2 * normcdf(1./(2*sqrt(D))) -1 );
    
    %
    Ps0=zeros(Pfttlslen,1);
    ns2=Ps0;
    g2=Ps0;
    optsub=zeros(Pfttlslen,1);
    FARmus=zeros(Pfttlslen,1);
    
    DoStep2flags=zeros(Pfttlslen,1);
    paras=struct('mu1',[],'mu2',[],'K1',[],'K2',[]);
    for ipf=1:Pfttlslen
        Pf1=Pfttls(ipf);%*0.9;
        Ps0(ipf)=1-Pf1;
        ns2(ipf)=Ps2ns(Psb,Ps0(ipf),na,D);
        
        k2=na-ns2(ipf)+1;
        Qzpar2 = Qzhat(k2:end,k2:end); Zpar2 = Z(:,k2:end); Qbzpar2=Qab'*Zpar2;
        Kcoef2=Qbzpar2/Qzpar2;
        Qbp2=Qb-Kcoef2*Qbzpar2';
        g2(ipf)=sqrt(det(Qb(1:3,1:3))/det(Qbp2(1:3,1:3)))^(1/3);%large than 1
        
        FARmus(ipf)=cpmu1(1-Psb,Pfttls(ipf),na);
        
        Pfreq1=Pf1s(ipf);
        %     find the optimal subset that maximizes Pfix*(g1-g2)
        EstPfixMethod=1;%use fitting function
        [optsub(ipf),DoStep2flags(ipf),paras(ipf)]=SearchOptsub(Pfttls(ipf),Pfreq1,Qzhat,Z,D,Qab,Qb,ns2(ipf),g2(ipf),0,0,EstPfixMethod);
        
        % calculate
        xhat=mvnrnd(xtrue,Qx,1)';
        ahat=xhat(1:na); bhat=xhat(na+1:na+nb);
        zhat=Z'*ahat;
        beflt=norm(bhat(1:3));
        fltrmss(ipf,epc)=beflt;
        [farbchkrmss(ipf,epc),farfixflags(ipf,epc),farsucflags(ipf,epc)]=DoMultFAR(zhat,ztrue,bhat,L,D,FARmus(ipf),Kfar,beflt);
        if ns2(ipf)==na
            srcbchkrmss(ipf,epc)=farbchkrmss(ipf,epc); tsrcbchkrmss(ipf,epc)=farbchkrmss(ipf,epc);
            srcfixflags(ipf,epc)=farfixflags(ipf,epc); tsrcfixflags(ipf,epc)=farfixflags(ipf,epc);
            srcsucflags(ipf,epc)=farsucflags(ipf,epc); tsrcsucflags(ipf,epc)=farsucflags(ipf,epc);
        else
            [srcbchkrmss(ipf,epc),srcfixflags(ipf,epc),srcsucflags(ipf,epc)]=DoMultSRC(zhat,ztrue,bhat,L,D,na,ns2(ipf),paras(ipf).K2,beflt);
            [tsrcbchkrmss(ipf,epc),tsrcfixflags(ipf,epc),tsrcsucflags(ipf,epc),tsrcfixns(ipf,epc)]=DoMultTSRC(zhat,ztrue,bhat,L,D,na,optsub(ipf),ns2(ipf),paras(ipf).mu1,paras(ipf).mu2,paras(ipf).K1,paras(ipf).K2,beflt);
        end
        resallARmethods=[resallARmethods;
            ep0,epc,Pf1,0,fltrmss(ipf,epc),0,0,0;
            ep0,epc,Pf1,1,farbchkrmss(ipf,epc),na,farfixflags(ipf,epc),farsucflags(ipf,epc);
            ep0,epc,Pf1,2,srcbchkrmss(ipf,epc),ns2(ipf),srcfixflags(ipf,epc),srcsucflags(ipf,epc);
            ep0,epc,Pf1,3,tsrcbchkrmss(ipf,epc),tsrcfixns(ipf,epc),tsrcfixflags(ipf,epc),tsrcsucflags(ipf,epc)];
        % do far
        % do tsrc
        % do src
    end
end

%method: 0--float, 1--FAR, 2--SRC, 3--TSRC
% display('finish');
% [ep0, epc, Pfttl, mtd, fltrms,ns, fixflag,sucflag]
% [ep0, epc, Pfttl, mtd, farrms,ns, fixflag,sucflag]
% [ep0, epc, Pfttl, mtd, srcrms,ns, fixflag,sucflag]
% [ep0, epc, Pfttl, mtd, tsrcrms,ns, fixflag,sucflag]



display(ep0);
%%End of the file-0222220**********************9
