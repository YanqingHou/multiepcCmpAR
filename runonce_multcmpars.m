function usdt=runonce_multcmpars(num)
%% Initialization
tic
S=load('options.mat');
opts=S.opts;
clear S;
filename=opts(num).filename;
Nsamp=opts(num).Nsamp;
Nsamp=1E6;%debug
loadfile=strcat('../Multsimumat/Mult_',filename,'.mat');
if ~exist(loadfile,'file')
    display(strcat('File doesnt exist: ',loadfile));
    usdt=0;
    return;
end
S=load(loadfile);
res=S.res; clear S;
lenres=size(res,1); %336
% lenres=48;%to shorten the running time.
savefilebench=strcat('../MultCmpARs/MultCmpARs_',filename,'.txt');

display('initial time');
toc


printformatstrs=cell(lenres,1);

epoch2simu=[];
for i=1:lenres
    tempfilebench=strcat('../MultCmpARs/MultCmpARsModel',num2str(num),'Epoch',num2str(i),'.txt');
    if ~exist(tempfilebench,'file')
        epoch2simu=[epoch2simu,i];
    end
end
epochsimulens=length(epoch2simu);


%% Parallel running
parfor epi=1:epochsimulens
%     for epi=1:lenres%:lenres%75%
    i=epoch2simu(epi);
    tempfilebench=strcat('../MultCmpARs/MultCmpARsModel',num2str(num),'Epoch',num2str(i),'.txt');
    [resbench]=multepc_compare_OneEp(i,res);
    printformatstr=makeSformat(size(resbench,2));
    printformatstrs{epi}=printformatstr;
    fid1=fopen(tempfilebench,'w');
    fprintf(fid1,printformatstr,resbench');
    fclose(fid1);
    
end
printformatstr=printformatstrs{1};
display('after parfor');
toc
%% Parallel saving result files
slicefiles1=cell(16,1);
% slicefiles2=cell(16,1);

lenslice=lenres/16;
parfor i=1:16
    % for i=1:16
    slicefiles1{i,1}=strcat('../MultCmpARs/MultCmpARsSliceModel',num2str(num),'Slice',num2str(i),'.txt');
    %     slicefiles2{i,1}=strcat('TSRCSliceModel',num2str(num),'Slice',num2str(i),'.txt');
    
    epcs=lenslice*(i-1)+1:lenslice*i;
    
    %     bench_slicewrite(slicefiles1{i,1},'CmpBench','',num,epcs);
    general_slicewrite(slicefiles1{i,1},'../MultCmpARs/MultCmpARs','',num,epcs,printformatstr);
    %     tsrcslicewrite(slicefiles2{i,1},'TSRC','',num,epcs);
    
    %     write to 16 slices
end
display('SliceFiles written');
toc
%% Combine all result files into one
% numslice=min(numslice1,numslice2,numslice3,numslice4);

fid1=fopen(savefilebench,'w');
% fid2=fopen(savefiletsrc,'w');

for i=1:16
    tempfilebench=slicefiles1{i,1};
    %     tempfiletsrc=slicefiles2{i,1};
    %     if ~exist(tempfileffrt)|| ~exist(tempfiletsrc)
    if ~exist(tempfilebench,'file')
        display(strcat('File doesnt exist: ',tempfilebench));
        continue;
    end
    reseprt01=load(tempfilebench);
    fprintf(fid1,printformatstr,reseprt01');
    delete(tempfilebench);
    
    
    %     write to one file
end

fclose(fid1);
% fclose(fid2);
display('all files combined');
toc
usdt=toc;
