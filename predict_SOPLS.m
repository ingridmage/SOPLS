function [Yhat,T,Yhat_blockwise]=predict_SOPLS(X,model)
%Yhat=predict_SOPLS(X,model)

nBlocks=length(X);
if nBlocks ~= length(model.X)
    error('Number of blocks in X must be same as in model')
end

n=size(X{1}.d,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preprocess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:nBlocks
    if isequal(model.options.preprocX{i},'mean center')
        X{i}.d=X{i}.d-repmat(mean(model.X{i}.d),n,1);
    elseif isequal(model.options.preprocX{i},'autoscale')
        X{i}.d=(X{i}.d-repmat(mean(model.X{i}.d),n,1))./repmat(std(model.X{i}.d),n,1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% predict
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%estimate scores for the first block
xx=X{1}.d;
T{1}=[];
    for j=1:model.options.nComps(1)
    T{1}(:,j)=xx*model.mod{1}.W(:,j);
    xx=xx-T{1}(:,j)*model.mod{1}.P(:,j)';
    end

for i=2:nBlocks
    
    xx=X{i}.d;
    %orthogonalise
    if ~isempty(model.Borth{i})
    xx=xx-cell2mat(T)*model.Borth{i};
    end
    
    %estimate scores
    T{i}=[];
    for j=1:model.options.nComps(i)
    T{i}(:,j)=xx*model.mod{i}.W(:,j);
    xx=xx-T{i}(:,j)*model.mod{i}.P(:,j)';
    end
    
end

if isequal(model.options.preprocY,'autoscale')
    stdY=std(model.Y.d);
else
    stdY=ones(1,size(model.Y.d,2));
end

Yhat = (cell2mat(T)*model.BetaT).*repmat(stdY,n,1)+repmat(mean(model.Y.d),n,1);

for i=1:nBlocks
    tt=cell2mat(T(1:i));
    if isempty(tt) %no components from the first blocks
        Yhat_blockwise{i}=repmat(mean(model.Y.d),n,1);
    else
Yhat_blockwise{i} = (tt*model.BetaT(1:size(tt,2),:)).*repmat(stdY,n,1)+repmat(mean(model.Y.d),n,1);
    end
end
