% ImageryLatencyLR(varargin)

%% dynamic classification across conditions (Imag)

% Get data (change path to your directory)
load('/Users/tysonaflalo/Dropbox/GitCode/Tyson/ForDistribution/eLifeTactile/ExpImagData.mat')

DC=DynCLass('overwrite',0,'Type','classification');

prefix='BasicAnalysis';



for i=1:length(ExpImagData)
basename=sprintf('%d-%s',i,prefix);

[AnalysisResults(i)]=DC.fit_transform(basename,ExpImagData(i).NeuralData,ExpImagData(i).Labels,ExpImagData(i).Time);

DC.plotAcc({AnalysisResults(i).cvAccuracy},AnalysisResults(i).timeWindow,'date',prefix)

end



%% Collapse across sessions
% Note, this includes imagery and touch to the right and left body sides;
% results in paper are only for imagery and right side.
cvAll=cell(3,3);
for i=1:length(AnalysisResults)
    tmp=AnalysisResults(i).cvAccuracy;
    for j1=1:3
        for j2=1:3
            cvAll{j1,j2}(:,:,i)=tmp{j1,j2};
        end
    end
end

for j1=1:3
    for j2=1:3
        cvAllMu{j1,j2}=mean(cvAll{j1,j2},3);
    end
end

%% Cross condition classification
DC.plotAcc({cvAllMu},round(AnalysisResults(i).timeWindow*10),'basename',prefix);


    
%% Generalization
 plt.fig('units','inches','width',5,'height',5,'font','Arial','fontsize',16);
        
            clear V
              for i=1:size(cvAll{j1,j2},3)
             
             V{2}(i,:)=(mean(cvAll{3,3}(18:22,:,i),1));

             end
      
             
        Analyze.plotEventRelatedAverage(V,{'early','late'},'useBootStrap')
        ylim([-20 150])
        
        xlim([25 65])
        ylim([20 90])
        
%         ylim([-2 2])
plt.hline(33.3,{'k--'})
        
set(gca, 'XTick',(25:5:65)+3)
tmp=AnalysisResults(1).timeWindow((25:5:65));
tmp(1)=-1;
        set(gca, 'XTickLabel',tmp)
        ylabel(sprintf('Cross-Validated \n Accuracy'))
        xlabel('Time (S)')
plt.vline(38,{'k--'})
title(sprintf('Generalization between \n Cue/Delay and Active Imagery'))

        %%
 