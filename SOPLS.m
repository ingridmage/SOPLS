function model=SOPLS(X,Y,options)
% options=SOPLS(X,Y)
% model=SOPLS(X,Y,options)
%
% INPUT
% X           cell array with X data blocks (saisir)
% Y           Y data (saisir)
% options:
%   nComps      Number of components from each block, empty if not known
%   preprocX    cell array with preprocessing method for each X block. Either
%               'mean center' or 'autoscale'
%   preprocY    preprocessing method for Y
%   compsel     3 x 1 cell array. component selection is performed
%               if nComps is empty
%               {1} 'global' or 'seq'
%               {2} 'manual' or 'auto'
%               {1} cvi. Either 'loo', 'r20', 'r10', or a vector with segments


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 %set default options
    
    for i=1:length(X)
    model.maxComps(i)=min(rank(X{i}.d),20);
    model.nComps =[];
    model.nCompAlternatives = [];
    model.preprocX{i}='mean center';
    end
    model.preprocY='mean center';
    model.BlockNames = strcat('Block',num2str((1:length(X))'));
    model.compsel{1} = 'global';
    model.compsel{2} = 'manual';
    model.compsel{3} = 'r10';
    

    
else
    
    nBlocks=length(X);
    nSamp=size(Y.d,1);
    nY=size(Y.d,2);
    
    if length(options.preprocX)~=nBlocks
        error('One preprocessing method for each X-block must be specified.')
    end
    
    
    
      



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialise model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.type = 'SO-PLS';
model.X=X; %save raw data
model.Y=Y;
model.options=options;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select components automatically if not given
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(options.nComps)
    
    if strcmp(lower(options.compsel{1}),'global')
        
        if strcmp(lower(options.compsel{2}),'manual')
            plots=1;
        else
            plots = 0;
        end
        globalopt_results=SOPLS_globalopt(X,Y,options,options.maxComps,options.compsel{3},plots);
        
        if strcmp(lower(options.compsel{2}),'manual')
            options.nComps = input('Number of components in each block (vector) : ')
        else
            options.nComps = globalopt_results.Aselected;
        end
        
        
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preprocess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:nBlocks
    X{i}=mypreprocess(X{i},options.preprocX{i});
    X{i}=rmfield(X{i},'pp');
end

Y=mypreprocess(Y,options.preprocY);
Y=rmfield(Y,'pp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model blocks sequentially
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExpVarY=zeros(nBlocks,nY);
Yresid=Y;
Ycal_blockwise=cell(nBlocks,1);

for i=1:nBlocks
    
    
    %sequential component selection
    if length(options.nComps)<i & contains(lower(options.compsel{1}),'seq')
          dummy=mypls(X{i},Yresid,options.maxComps(i),'none','none');
         dummy = mycrossval(dummy,options.compsel{3},0);
         if strcmp(lower(options.compsel{2}),'manual')
             plsplots(dummy)
             options.nComps(i) = input(['Number of components from  ' options.BlockNames{i} ' : '])
         else
             options.nComps(i)=dummy.Aopt;
         end
         
    end
    
    %fit model
    plsmod{i}=mypls(X{i},Yresid,options.nComps(i),'none','none');
   
    T{i}=plsmod{i}.T(:,1:options.nComps(i));
    
    if i<nBlocks
    %orthogonolize the next block based on all previous
    Tcurrent=cell2mat(T);
    model.Borth{i+1}=inv(Tcurrent'*Tcurrent)*Tcurrent'*X{i+1}.d;
    X{i+1}.d=X{i+1}.d-Tcurrent*model.Borth{i+1};
    end
    
    model.mod{i}=plsmod{i};
    xx=cell2mat(T); %scores from all previous blocks
    Ycal=xx*inv(xx'*xx)*xx'*Y.d;
    Ycal_blockwise{i}=Ycal;
    SSQres_blockwise(i,:)=sum((Y.d-Ycal).^2);
    Yresid.d=Y.d-Ycal;
end

%final model
xx=cell2mat(T);
Beta=inv(xx'*xx)*xx'*Y.d;
Ycal=xx*Beta;

%rescale Ycal
if isequal(model.options.preprocY,'autoscale')
    stdY=std(model.Y.d);
else
    stdY=ones(1,size(model.Y.d,2));
end
Ycal=Ycal.*repmat(stdY,nSamp,1)+repmat(mean(model.Y.d),nSamp,1);
for i=1:nBlocks
Ycal_blockwise{i} = Ycal_blockwise{i}.*repmat(stdY, nSamp,1)+repmat(mean(model.Y.d),nSamp,1);
end

%calculate explained variance
SSQres=(sum((model.Y.d - Ycal).^2));
SSQ0=(sum((Y.d.^2)));
ExpVarY=(1-SSQres./SSQ0)*100;
ExpVarY_blockwise=(1-SSQres_blockwise./repmat(SSQ0,nBlocks,1))*100;

ExpVarYtot=(1-sum(SSQres)./sum(SSQ0))*100;
ExpVarYtot_blockwise=(1-sum(SSQres_blockwise,2)./repmat(sum(SSQ0),nBlocks,1))*100;



model.options=options;
model.T=T;
model.BetaT=Beta;
model.Ycal=Ycal;
%model.Ycal_blockwise=Ycal_blockwise;
model.ExpVarY=ExpVarY;
model.ExpVarY_blockwise=ExpVarY_blockwise;
model.ExpVarYtot=ExpVarYtot;
model.ExpVarYtot_blockwise=ExpVarYtot_blockwise;

model.RMSEcal=sqrt(mean((model.Y.d-model.Ycal).^2));

end