
%% data

T1=randn(100,3); 
T2=[T1(:,1) randn(100,2)];

X1=mat2saisir(T1*randn(3,10)+randn(100,10)*0.1);
X2=mat2saisir(T2*randn(3,500)+randn(100,500)*0.1);

Y=mat2saisir([T1(:,1:2) T2(:,2:3)]*[1 1 1 1]'+randn(100,1)*0.1);




%% set options

options=SOPLS({X1 X2},Y);
options.preprocX{1}='autoscale'

%% select optimal number of components (global method)



Amax=[5 5];
cvi='r10'
plots=1;
globalopt_results=SOPLS_globalopt({X1 X2},Y,options,Amax,cvi,plots)

%% fit model

options.nComps=[2 2]

model=SOPLS({X1 X2},Y,options)
model=crossval_SOPLS(model,'loo')

%% predict new data

[Yhat,T,Yhat_blockwise]=predict_SOPLS({X1 X2},model)

%% plot

