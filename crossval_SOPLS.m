function model=crossval_SOPLS(model,cvi,wb)
%model=crossval_SOPLS(model,cvi,wb)

X=model.X;
Y=model.Y;
nSamp=size(X{1}.d,1);
nBlocks=length(X);
nY=size(Y.d,2);
options=model.options;

if nargin==2
    wb=1;
end

if ischar(cvi)
    %set cv segments
    switch lower(cvi)
        case 'loo'
            cvi=1:nSamp;
        case 'r10'
            nn=floor(nSamp/10);
            rest=nSamp-10*nn;
            cvi=[repmat(1:10,1,nn) 1:rest];
            cvi=cvi(randperm(nSamp));
        case 'r20'
            nn=floor(nSamp/20);
            rest=nSamp-20*nn;
            cvi=[repmat(1:20,1,nn) 1:rest];
            cvi=cvi(randperm(nSamp));
        otherwise
            if length(cvi)~=nSamp
                error('cvi segments does not match number of samples')
            end
    end
end


Yval=zeros(nSamp,nY);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preprocess outside the loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:nBlocks
    X{i}=mypreprocess(X{i},options.preprocX{i});
    X{i}=rmfield(X{i},'pp');
end

Y=mypreprocess(Y,options.preprocY);
Y=rmfield(Y,'pp');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if wb==1
    H=waitbar(0,['SO-PLS cross-validation with ' num2str(max(cvi)) ' segments...']);
end

Yhat=zeros(nSamp,nY,nBlocks+1);
Yval_blockwise = cell(nBlocks,1);
%cross-validation loop
for i=1:max(cvi)
    if wb==1
        waitbar(i/max(cvi),H)
    end
    idx=find(cvi==i);
    
    %keep out and keep in
    Xcal=cell(1,nBlocks);
    Xval=cell(1,nBlocks);
    for j=1:nBlocks
        Xcal{j}=deleterow(X{j},idx);
        Xval{j}=selectrow(X{j},idx);
    end
    Ycal=deleterow(Y,idx);
    
    
    %calibrate and predict
    modi=SOPLS(Xcal,Ycal,options);
    [Yval(idx,:),T,Yval_blockwise_tmp]=predict_SOPLS(Xval,modi);
    
    for k=1:nBlocks
       Yval_blockwise{k}(idx,:)=Yval_blockwise_tmp{k};
    end
    
    Yval0(idx,:)=repmat(mean(Ycal.d),length(idx),1);
    
end


if wb==1
    close(H)
end

if isequal(model.options.preprocY,'autoscale')
    stdY=std(model.Y.d);
else
    stdY=ones(1,size(model.Y.d,2));
end



%rescale predicted values
Yval = Yval.*repmat(stdY,nSamp,1)+repmat(mean(model.Y.d),nSamp,1);
Yval0 = Yval0.*repmat(stdY,nSamp,1)+repmat(mean(model.Y.d),nSamp,1);
for i=1:nBlocks
Yval_blockwise{i} = Yval_blockwise{i}.*repmat(stdY, nSamp,1)+repmat(mean(model.Y.d),nSamp,1);
SSQres_blockwise(i,:)=sum((model.Y.d-Yval_blockwise{i}).^2);
end


%Calculate explained variances
SSQres=(sum((model.Y.d - Yval).^2));
SSQ0=sum((model.Y.d-Yval0).^2);

ExpVarY=(1-SSQres./SSQ0)*100;
ExpVarY_blockwise=(1-SSQres_blockwise./repmat(SSQ0,nBlocks,1))*100;

ExpVarYtot=(1-sum(SSQres)./sum(SSQ0))*100;
ExpVarYtot_blockwise=(1-sum(SSQres_blockwise,2)./repmat(sum(SSQ0),nBlocks,1))*100;



model.cvres.Yval=Yval;
model.cvres.Yval_blockwise=Yval_blockwise;
model.cvres.ExpVarY=ExpVarY;
model.cvres.ExpVarY_blockwise=ExpVarY_blockwise;
model.cvres.ExpVarYtot=ExpVarYtot;
model.cvres.ExpVarYtot_blockwise=ExpVarYtot_blockwise;
model.cvres.cvi = cvi;


model.cvres.RMSEcv=sqrt(mean((model.Y.d-Yval).^2));



end
