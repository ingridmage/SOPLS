function plotoptions=plot_SOPLS(model,plotoptions)
% plotoptions=plot_SOPLS(model)
% plot_SOPLS(model,plotoptions)

%Modifications, aug 2013 by Ingunn Berget
%some changes for readability in making labels (use format in num2str)
%dealing with one component models (using additional args for diff)
%options to create different types of plots

%modifications march 2014 by Ingrid Måge
% collect options in plotoptions struct with fields:
% plots:            'all' (default),'expVarPie','Predplot' or 'ScoresLoadings'.
% spec:             =1 if data are spectra 
% PCs:              which two components to plot. Default is [1 2]
% samplegroups:     numeric vector or character array defining sample groups
% variablegroups:   cell array containing numeric array or character array
%                   defining variable groups.

if nargin<2
    plotoptions.plots='all';
    plotoptions.spec=0;
    plotoptions.PCs=[1 2];
    plotoptions.Yvar=[];
    plotoptions.samplegroups=[];
    plotoptions.variablegroups=[];
else
    
    
    
    
    if strcmp('all',plotoptions.plots) == 1 | strcmp('expVarPie',plotoptions.plots)
        plotExpVarY(model)
    end
    
    if strcmp('all',plotoptions.plots) == 1 | strcmp('PredPlot',plotoptions.plots)
        makePredPlot(model)
    end
%     
%     if strcmp('all',plotoptions.plots) == 1 | strcmp('ScoresLoadings',plotoptions.plots)
%         plot_scores_loadings(model,plotoptions.spec,plotoptions.samplegroups,plotoptions.variablegroups)
%     end
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INTERNAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plotExpVarY(model)
        %if less than 20 responses:
        %total variance for each response (bar)
        %and pie chart of explained variance for each block combinations  if less
        %than 20 responses
        %for models with only one comp: bar plot
        %
        
        nBlocks=length(model.X);
        nY=size(model.Y.d,2);
        
        for i=1:nY
        BlockEV = diff([0; model.cvres.ExpVarY_blockwise(:,i)]);
        BlockEV(BlockEV==0)=0.0001;

        if isfield(model.options,'BlockNames')
            labels = model.options.BlockNames;
        else
            labels = strcat('Block',num2str((1:nBlocks)'));
        end
if iscell(labels)
    labels{end+1} = 'Residuals';
else
labels = cellstr(strvcat(labels,'Residuals'));
end
BlockEV = [BlockEV; 100-sum(BlockEV)];
for j=1:length(labels)
labels{j} = [labels{j} ' (' num2str(round(BlockEV(j),1)) '%)'];
end
figure; h=pie(BlockEV/100,[],labels); 
%title(strvcat(model.Y.v(i,:),'Cross-validation'));
h(end-1).FaceColor = [1 1 1];
        
       
        
        end
    end %function end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predicted versus measured
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function makePredPlot(model)
        
        %model: From POSO_PLS
        %Yvar: Optional input in main function
        %if given, look at only selected y variables
        %if not given, look at all when number is less than 20
        
        nBlocks=length(model.X);
        nY=size(model.Y.d,2);
        
        for i=1:nY
            figure; 
            scatter(model.Y.d(:,i),model.Ycal(:,i),'filled');
            hold on;
            scatter(model.Y.d(:,i),model.cvres.Yval(:,i),'filled');
            legend(['Calibration, R^2 = ' num2str(round(model.ExpVarY(i))/100)],['Cross-validation, R^2 = ' num2str(round(model.cvres.ExpVarY(i))/100)],'Autoupdate','off','Location','northwest');
            xlabel('Reference'); ylabel('Predicted');
            ymin = min([model.Y.d(:,i);model.cvres.Yval(:,i)]);
            ymax = max([model.Y.d(:,i);model.cvres.Yval(:,i)]);
            line([ymin ymax],[ymin ymax],'color','k');
            axis tight
            title(strcat(model.Y.v(i,:)));
        end
        
    end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scores and loadings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_scores_loadings(model,spec,sample_groups,variable_groups)
        col=colormap(lines);%strvcat('b','m','c','g','k','y');
        col(3,:)=[]; %remove the red color, which is reserved for y loadings
        
        %make one figure for each contribution
        for i=1:length(model.options.blockCombinations) %each contribution
            
            figure('Name',['Scores & Loadings, ' model.ContrLabel{i}]);
            
            nComps=size(model.Scores{i},2);
            
            if nComps==1 %only one component
                
                %scores
                if spec==0
                    subplot(1,2,1)
                else
                    subplot(1,1+length(model.options.blockCombinations{i}),1)
                end
                
                if isempty(sample_groups) %no sample grouping
                    plot(zeros(size(model.Scores{i})),model.Scores{i},'.')
                    hold on
                    text(zeros(size(model.Scores{i}))+0.03,model.Scores{i},model.X{1}.i)
                else
                    sample_grouping(zeros(size(model.Scores{i})),model.Scores{i},sample_groups,model.X{1}.i);
                    colorbar off
                end
                ylabel(model.ContrLabel{i})
                title(['Scores ' model.ContrLabel{i}])
                
                %loadings
                if spec==0
                    subplot(1,2,2)
                    title(['Correlation Loadings ' model.ContrLabel{i}])
                    
                    idx=model.options.blockCombinations{i};
                    for j=1:length(idx)
                        
                        if isempty(variable_groups)
                            h=plot(zeros(size(model.X{idx(j)}.d,2),1),corr(model.Scores{i}(:,1),model.X{idx(j)}.d),'.');
                            set(h,'MarkerEdgeColor',col(idx(j),:))
                            hold on
                            H=text(zeros(size(model.X{idx(j)}.d,2),1),corr(model.Scores{i}(:,1),model.X{idx(j)}.d),model.X{idx(j)}.v,'Color',col(idx(j),:));
                        else
                            sample_grouping(zeros(size(model.X{idx(j)}.d,2),1),corr(model.Scores{i}(:,1),model.X{idx(j)}.d),variable_groups{idx(j)},model.X{idx(j)}.v)
                            colorbar off
                            hold on
                        end
                        
                    end
                    plot(zeros(size(model.Y.d,2),1),corr(model.Scores{i}(:,1),model.Y.d),strcat('.','r'))
                    text(zeros(size(model.Y.d,2),1),corr(model.Scores{i}(:,1),model.Y.d),strcat(upper(model.Y.v)),'Color','r','EdgeColor','r');
                    
                    H=line(get(gca,'xlim'),[0 0]);
                    set(H,'linestyle',':','Color','k')
                    H=line(get(gca,'xlim'),[-sqrt(0.5) -sqrt(0.5)]);
                    set(H,'linestyle',':','Color','k')
                    H=line(get(gca,'xlim'),[sqrt(0.5) sqrt(0.5)]);
                    set(H,'linestyle',':','Color','k')
                else %spec=1
                    
                    idx=model.options.blockCombinations{i};
                    
                    
                    for k=1:length(idx)
                        
                        subplot(1,1+length(idx),1+k)
                        
                        if length(idx)>1
                            plot(str2num(model.X{idx(k)}.v),model.Loadings{i}{idx(k)});
                            axis tight
                            if ~isempty(variable_groups)
                                sample_grouping(str2num(model.X{idx(k)}.v),model.Loadings{i}{idx(k)},variable_groups{idx(k)},repmat(' ',length(model.X{idx(k)}.v),1));
                            end
                        else
                            plot(str2num(model.X{idx(k)}.v),model.Loadings{i});
                            axis tight
                            if ~isempty(variable_groups)
                                sample_grouping(str2num(model.X{idx(k)}.v),model.Loadings{i},variable_groups{idx(k)},repmat(' ',length(model.X{idx(k)}.v),1));
                            end
                        end
                        colorbar off
                        legend off
                        title([model.ContrLabel{i} '-Loadings, Block ' num2str(idx(k))])
                        
                    end
                    
                end
                
            else % more than one component --> scatter plot
                nplot=1;
                
                nPlots=floor(nComps/2);
                compCombinations=[];
                dummy=1:nComps;
                if ~isequal(length(dummy)/2,floor(length(dummy)/2))
                    dummy=[dummy(1:end-1) 1 dummy(end)];
                end
                compCombinations=reshape(dummy,2,length(dummy)/2)';
                
                for j=1:nPlots
                    
                    %scores
                    if spec==0
                        subplot(nPlots,2,nplot)
                    else
                        subplot(nPlots,1+length(model.options.blockCombinations{i}),nplot)
                    end
                    
                    title(['Scores, ' model.ContrLabel{i}])
                    if isempty(sample_groups)
                        plot(model.Scores{i}(:,compCombinations(j,1)),model.Scores{i}(:,compCombinations(j,2)),'.')
                        hold on
                        text(model.Scores{i}(:,compCombinations(j,1)),model.Scores{i}(:,compCombinations(j,2)),model.Y.i)
                        
                    else
                        sample_grouping(model.Scores{i}(:,compCombinations(j,1)),model.Scores{i}(:,compCombinations(j,2)),sample_groups,model.Y.i);
                        colorbar off
                    end
                    
                    %loadings
                    if spec==0
                        subplot(nPlots,2,nplot+1)
                        idx=model.options.blockCombinations{i};
                        for k=1:length(idx)
                            
                            if isempty(variable_groups)
                                h=plot(corr(model.Scores{i}(:,compCombinations(j,1)),model.X{idx(k)}.d),corr(model.Scores{i}(:,compCombinations(j,2)),model.X{idx(k)}.d),'.');
                                set(h,'MarkerEdgeColor',col(idx(k),:))
                                hold on
                                text(corr(model.Scores{i}(:,compCombinations(j,1)),model.X{idx(k)}.d),corr(model.Scores{i}(:,compCombinations(j,2)),model.X{idx(k)}.d),model.X{idx(k)}.v,'Color',col(idx(k),:));
                            else
                                
                                sample_grouping(corr(model.Scores{i}(:,compCombinations(j,1)),model.X{idx(k)}.d),corr(model.Scores{i}(:,compCombinations(j,2)),model.X{idx(k)}.d),variable_groups{idx(k)},model.X{idx(k)}.v)
                                hold on
                                colorbar off
                            end
                        end
                        
                        plot(corr(model.Scores{i}(:,compCombinations(j,1)),model.Y.d),corr(model.Scores{i}(:,compCombinations(j,2)),model.Y.d),'r.')
                        h=text(corr(model.Scores{i}(:,compCombinations(j,1)),model.Y.d),corr(model.Scores{i}(:,compCombinations(j,2)),model.Y.d),strcat(upper(model.Y.v)),'Color','r','EdgeColor','r');
                        %circle 100%
                        [X,Y] = pol2cart(linspace(0,2*pi,100),ones(1,100));
                        plot(X,Y,'k');
                        %circle 50%
                        [X,Y] = pol2cart(linspace(0,2*pi,100),ones(1,100)*sqrt(0.5));
                        plot(X,Y,'k');
                        H=line([0 0],[-1 1]);
                        set(H,'linestyle',':','Color','k')
                        H=line([-1 1],[0 0]);
                        set(H,'linestyle',':','Color','k')
                        title(['Correlation loadings, ' model.ContrLabel{i}])
                        
                    else %spec=1
                        
                        idx=model.options.blockCombinations{i};
                        
                        
                        for k=1:length(idx)
                            subplot(nPlots,1+length(idx),nplot+k)
                            
                            if length(idx)>1
                                plot(str2num(model.X{idx(k)}.v),model.Loadings{i}{idx(k)});
                                axis tight
                                if ~isempty(variable_groups)
                                    for ii=1:size(model.Loadings{i}{idx(k)},2)
                                        sample_grouping(str2num(model.X{idx(k)}.v),model.Loadings{i}{idx(k)}(:,ii),variable_groups{idx(k)},repmat(' ',length(model.X{idx(k)}.v),1));
                                        hold on
                                    end
                                    
                                    
                                end
                                
                            else %one block
                                
                                if ~isempty(variable_groups)
                                    for ii=1:size(model.Loadings{i},2)
                                        sample_grouping(str2num(model.X{idx(k)}.v),model.Loadings{i}(:,ii),variable_groups{idx(k)},repmat(' ',length(model.X{idx(k)}.v),1));
                                        hold on
                                    end
                                    
                                    
                                end
                                plot(str2num(model.X{idx(k)}.v),model.Loadings{i});
                                axis tight
                            end
                            colorbar off
                            legend off
                            title([model.ContrLabel{i} '-Loadings, Block ' num2str(idx(k))])
                        end
                        
                        
                        
                    end
                end
                nplot=nplot+2;
                
            end
        end
        
    end
end %function

