% ImageryLatencyLR(varargin)
%%
% Dates = {'20180702','20180709','20180711','20180716','20180718'}; % new
Dates = {'20180709','20180711','20180716','20180718','20180723','20180730','20180806','20180815'}; % new, for imagery
Labels = {'CheekL','CheekR','ImagCheek','ShoulderL','ShoulderR','ImagShoulder','HandL','HandR','ImagHand','XX'};
Phase = 'Go';


%% DEfine params

Dates = {'20180709','20180711','20180716','20180718','20180723','20180730','20180806','20180815'}; % new, for imagery

DF=DataFilters.FaceScr_ImagExp_PopDist('winSize',.5);
DC=Analyze.FaceScratch.DynDistance('overwrite',1);

prefix='PopDist'; % classifiation .5 kernel
%% do actual processing of distances.
parfor i=1:length(Dates)

cDate=Dates{i};

basename=sprintf('%s-%s',cDate,prefix);

AllTrialData=Analyze.LoadConvertedData('FaceScratch',cDate);


[FRTrue,LabelTrue,timeWindow]=DF.fit_transform(AllTrialData);

[cvAccuracy]=DC.fit_transform(basename,FRTrue,LabelTrue,timeWindow);

DC.plotAcc({cvAccuracy.cvAccuracy},cvAccuracy.timeWindow,LabelTrue,'date',basename)


end

%% load saved per session distances and combine

Dates = {'20180709','20180711','20180716','20180718','20180723','20180730','20180806','20180815'}; % new, for imagery

for i=1:length(Dates)
cDate=Dates{i};
basename=sprintf('%s-%s',cDate,prefix)

AllTrialData=Analyze.LoadConvertedData('FaceScratch',cDate)

DF=DataFilters.FaceScr_ImagExp_PopDist('winSize',.5);
DC=Analyze.FaceScratch.DynDistance('overwrite',0);

[cvAccuracyAll{i}]=DC.fit_transform(basename,[],[],[]);
end

%
cvAll=cell(2,2); clear cvAllMu
 for i=1:length(cvAccuracyAll)     
     tmp=cvAccuracyAll{i}.cvAccuracy;
     for j1=1:2
         for j2=1:2
        cvAll{j1,j2}(:,:,i)=tmp{j1,j2};     
         end
     end
 end

for j1=1:2
         for j2=1:2
        cvAllMu{j1,j2}=mean(cvAll{j1,j2},3);
         end
     end    

%%
DC.plotAcc({cvAllMu},round(cvAccuracyAll{1}.timeWindow*10)-3,LabelTrue,'basename',prefix);

%%
 plt.fig('units','inches','width',5,'height',5,'font','Arial','fontsize',12);
            pnl = panel();  pnl.margin=10; pnl.pack(2,2); pnl.fontsize=12;pnl.fontname='arial';
            
for j1=1:2
         for j2=1:2
             pnl(j1,j2).select()
             clear V
             for i=1:size(cvAll{j1,j2},3)
             V(i,:)=diag(cvAll{j1,j2}(:,:,i));
             end
        Analyze.plotEventRelatedAverage({V},{''},'useBootStrap')
        ylim([-20 150])
%         plt.hline(30,{'k--'})
         end
end   
%%
     
%%
 plt.fig('units','inches','width',5,'height',5,'font','Arial','fontsize',16);
        
            clear V
              for i=1:size(cvAll{j1,j2},3)
%              V{1}(i,:)=mean(cvAll{1,1}(20:35,:,i),1);
%              V{2}(i,:)=mean(cvAll{2,2}(20:35,:,i),1);
             
             V{1}(i,:)=(mean(cvAll{1,1}(20,:,i),1));
             V{2}(i,:)=(mean(cvAll{2,2}(20,:,i),1));
%              V(i,:)=mean(cvAll{3,3}(:,28:37,i)',1);
             end
      
             
        Analyze.plotEventRelatedAverage(V,{'early','late'},'useBootStrap')
        ylim([-20 150])
        
        xlim([25 55])
        ylim([-10 90])
        
%         ylim([-2 2])

        
set(gca, 'XTick',(25:5:55)+3)
tmp=cvAccuracyAll{1}.timeWindow((25:5:55));
tmp(1)=-1;
        set(gca, 'XTickLabel',tmp)
        ylabel(sprintf('Cross-Validated \n Mahalanobis Distance'))
        xlabel('Time (S)')
plt.vline(38,{'k--'})
title(sprintf('Seperable Neural Processes between \n Cue/Delay and Active Imagery'))
        plt.SaveFigure(1,'/Users/tysonaflalo/Dropbox/Lab/Papers/NS tactile/NewFigs','MahDistPerEffect','PNG','PDF','SVGI')
        %%
 plt.fig('units','inches','width',5,'height',5,'font','Arial','fontsize',16);
        
            clear V
              for i=1:size(cvAll{j1,j2},3)
%            
             V{1}(i,:)=(mean(cvAll{1,1}(20,:,i),1));
             V{2}(i,:)=(mean(cvAll{2,2}(20,:,i),1));
                          
             end
      
             
        Analyze.plotEventRelatedAverage({cat(1,V{:})},{'early'},'useBootStrap')
        ylim([-20 150])
        
        xlim([25 65])
        set(gca, 'XTick',(25:5:65)+3)

        ylim([-5 70])
        
%         ylim([-2 2])

        tmp=cvAccuracyAll{1}.timeWindow((25:5:65));
tmp(1)=-1;
        set(gca, 'XTickLabel',tmp)
        ylabel(sprintf('Cross-Validated \n Mahalanobis Distance'))
        xlabel('Time (S)')
plt.vline(38,{'k--'})
title(sprintf('Seperable Neural Processes between \n Cue/Delay and Active Imagery'))
plt.SaveFigure(1,'/Users/tysonaflalo/Dropbox/Lab/Papers/NS tactile/NewFigs','MahDist','PNG','PDF','SVGI')

%%
figure
tmp=mean(cat(1,V{:}),1);
findchangepts(tmp(25:65),'MaxNumChanges',2,'statistic','mean')

%%
        %%
 plt.fig('units','inches','width',5,'height',5,'font','Arial','fontsize',16);
        
            clear V
              for i=1:size(cvAll{j1,j2},3)
%            
             V{1}(i,:)=(mean(cvAll{1,1}(20,:,i),1));
             V{2}(i,:)=(mean(cvAll{2,2}(20,:,i),1));
                          
             end
      
             
        Analyze.plotEventRelatedAverage({cat(1,V{:})},{'early'},'useBootStrap')
        ylim([-20 150])
        
        xlim([25 65])
        set(gca, 'XTick',(25:5:65)+3)

        ylim([-5 70])
        
%         ylim([-2 2])

        tmp=cvAccuracyAll{1}.timeWindow((25:5:65));
tmp(1)=-1;
        set(gca, 'XTickLabel',tmp)
        ylabel(sprintf('Cross-Validated \n Mahalanobis Distance'))
        xlabel('Time (S)')
plt.vline(38,{'k--'})
title(sprintf('Seperable Neural Processes between \n Cue/Delay and Active Imagery'))
% plt.SaveFigure(1,'/Users/tysonaflalo/Dropbox/Lab/Papers/NS tactile/NewFigs','MahDist_Gen','PNG','PDF','SVGI')