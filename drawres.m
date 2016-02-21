% drawres
figure;
subplot(2,1,1);
plot(1:50,farbchkrmss(1,:),'r-*',1:50,srcbchkrmss(1,:),'b-o',1:50,tsrcbchkrmss(1,:),'m-p');
legend('FAR','SRC','TSRC');
subplot(2,1,2);
plot(1:50,farsucflags(1,:),'r-*',1:50,srcsucflags(1,:),'b-o',1:50,tsrcsucflags(1,:),'m-p');
legend('FAR','SRC','TSRC');
