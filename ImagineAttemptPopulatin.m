% ImageryLatencyLR(varargin)
%%
% Dates = {'20180702','20180709','20180711','20180716','20180718'}; % new
Dates = {'20180709','20180711','20180716','20180718','20180723','20180730','20180806','20180815'}; % new, for imagery
Phase = 'Go';
Labels = {'CheekL','CheekR','ImagCheek','ShoulderL','ShoulderR','ImagShoulder','HandL','HandR','ImagHand','XX'};
Tag = 'XImagery2';
GlobalTrialCriteria = {'TaskLabel',Tag};
OutDir = fullfile(env.get('result'),'FaceScratch','SUAnal',[Tag '-Go']);
% don't include imagine conditions in baseline (because working memory/delay activity)
BaselinePhase = {'Phase','Delay','Condition',{'CheekL','CheekR','ShoulderL','ShoulderR','HandL','HandR','XX'}};
BaselineWindow = [0 1]; % relative to delay phase
PhaseStart = .5; % for Go phase
PhaseDur = 2; %[0 - 2]?
PhaseDur = 1; %[0 - 2]?


%% dynamics test

Dates = {'20180709','20180711','20180716','20180718','20180723','20180730','20180806','20180815'}; % new, for imagery

DF=DataFilters.FaceScr_ImagExp('winSize',.5);
DC=Analyze.FaceScratch.DynCLass('overwrite',1,'Type','correlation');

prefix='wide'; % classifiation .5 kernel
% prefix='corr'; % correlation

parfor i=1:length(Dates)
cDate=Dates{i};

basename=sprintf('%s-%s',cDate,prefix)

AllTrialData=Analyze.LoadConvertedData('FaceScratch',cDate)


[FRTrue,LabelTrue,timeWindow]=DF.fit_transform(AllTrialData);

[cvAccuracy]=DC.fit_transform(basename,FRTrue,LabelTrue,timeWindow);

DC.plotAcc({cvAccuracy.cvAccuracy},cvAccuracy.timeWindow,'date',basename)


end

%% Combine

Dates = {'20180709','20180711','20180716','20180718','20180723','20180730','20180806','20180815'}; % new, for imagery
prefix='wide'; % classifiation .5 kernel

for i=1:length(Dates)
cDate=Dates{i};
basename=sprintf('%s-%s',cDate,prefix)

AllTrialData=Analyze.LoadConvertedData('FaceScratch',cDate)

DF=DataFilters.FaceScr_ImagExp();
DC=Analyze.FaceScratch.DynCLass('overwrite',0);

[cvAccuracyAll{i}]=DC.fit_transform(basename,[],[],[]);
end

%
cvAll=cell(3,3)
 for i=1:length(cvAccuracyAll)     
     tmp=cvAccuracyAll{i}.cvAccuracy;
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

%%
DC.plotAcc({cvAllMu},round(cvAccuracyAll{1}.timeWindow*10)-3,'basename',prefix);

%%
 plt.fig('units','inches','width',5,'height',5,'font','Arial','fontsize',12);
            pnl = panel();  pnl.margin=10; pnl.pack(3,3); pnl.fontsize=12;pnl.fontname='arial';
            
for j1=1:3
         for j2=1:3
             pnl(j1,j2).select()
             clear V
             for i=1:size(cvAll{j1,j2},3)
             V(i,:)=diag(cvAll{j1,j2}(:,:,i));
             end
        Analyze.plotEventRelatedAverage({V},{''},'useBootStrap')
        ylim([20 100])
        plt.hline(30,{'k--'})
         end
end   
%%
     
%%
 plt.fig('units','inches','width',5,'height',5,'font','Arial','fontsize',12);
        
            clear V
              for i=1:size(cvAll{j1,j2},3)
             V{1}(i,:)=mean(cvAll{3,3}(18:35,:,i),1);
             V{2}(i,:)=mean(cvAll{3,3}(50:70,:,i),1);
%              V(i,:)=mean(cvAll{3,3}(:,28:37,i)',1);
             end
      
             
        Analyze.plotEventRelatedAverage(V,{'early','late'},'useBootStrap')
        ylim([20 100])
        plt.hline(30,{'k--'})

%%
 plt.fig('units','inches','width',5,'height',5,'font','Arial','fontsize',16);
        
            clear V
              for i=1:size(cvAll{j1,j2},3)
%              V{1}(i,:)=mean(cvAll{1,1}(20:35,:,i),1);
%              V{2}(i,:)=mean(cvAll{2,2}(20:35,:,i),1);
             
             V{2}(i,:)=(mean(cvAll{3,3}(18:22,:,i),1));
%              V(i,:)=mean(cvAll{3,3}(:,28:37,i)',1);
             end
      
             
        Analyze.plotEventRelatedAverage(V,{'early','late'},'useBootStrap')
        ylim([-20 150])
        
        xlim([25 65])
        ylim([20 90])
        
%         ylim([-2 2])
plt.hline(33.3,{'k--'})
        
set(gca, 'XTick',(25:5:65)+3)
tmp=cvAccuracyAll{1}.timeWindow((25:5:65));
tmp(1)=-1;
        set(gca, 'XTickLabel',tmp)
        ylabel(sprintf('Cross-Validated \n Accuracy'))
        xlabel('Time (S)')
plt.vline(38,{'k--'})
title(sprintf('Generalization between \n Cue/Delay and Active Imagery'))
        plt.SaveFigure(1,'/Users/tysonaflalo/Dropbox/Lab/Papers/NS tactile/NewFigs','EarlyLateGen','PNG','PDF','SVGI')
%%
 plt.fig('units','inches','width',5,'height',5,'font','Arial','fontsize',16);
        
            clear V
              for i=1:size(cvAll{j1,j2},3)
%    
             % Imagery_> Exp
             V{1}(i,:)=(mean(cvAll{3,1}(:,45:65,i),2));
             V{2}(i,:)=(mean(cvAll{3,2}(:,45:65,i),2));
             
             
             
             %Exp2Imag
%              V{1}(i,:)=(mean(cvAll{1,3}(:,35:65,i),2));
%              V{2}(i,:)=(mean(cvAll{2,3}(:,35:65,i),2));
             
%             V{1}(i,:)=(mean(cvAll{1,3}(44:66,:,i),1));
%              V{2}(i,:)=(mean(cvAll{2,3}(44:66,:,i),1));
             end
      
             
        Analyze.plotEventRelatedAverage(V,{'R','L'},'useBootStrap','Legend')
        ylim([-20 150])
        

        ylim([15 70])
        plt.hline(33.3,{'k--'})

        if 1==1
            bd=[25 65];
                xlim(bd)
set(gca, 'XTick',(bd(1):5:bd(2))+3)
tmp=cvAccuracyAll{1}.timeWindow((bd(1):5:bd(2)));
tmp(1)=-1;
        set(gca, 'XTickLabel',tmp)
        end
        ylabel(sprintf('Generalization Accuracy '))
        xlabel('Imagery Training Window')
plt.vline(38,{'k--'})
title(sprintf('Generalization from \n Imagery to Experience'))
        plt.SaveFigure(1,'/Users/tysonaflalo/Dropbox/Lab/Papers/NS tactile/NewFigs','GenImag2Exp','PNG','PDF','SVGI')
        
        %%
 plt.fig('units','inches','width',5,'height',5,'font','Arial','fontsize',16);
        
            clear V
              for i=1:size(cvAll{j1,j2},3)
%    
             % Imagery_> Exp
             V{1}(i,:)=(mean(cvAll{3,1}(:,45:65,i),2));
%              V{2}(i,:)=(mean(cvAll{3,2}(:,45:65,i),2));
             
             
             
             %Exp2Imag
%              V{1}(i,:)=(mean(cvAll{1,3}(:,35:65,i),2));
%              V{2}(i,:)=(mean(cvAll{2,3}(:,35:65,i),2));
             
%             V{1}(i,:)=(mean(cvAll{1,3}(44:66,:,i),1));
%              V{2}(i,:)=(mean(cvAll{2,3}(44:66,:,i),1));
             end
      
             
        Analyze.plotEventRelatedAverage(V,{'R'},'useBootStrap','Legend')
        ylim([-20 150])
        

        ylim([15 70])
        plt.hline(33.3,{'k--'})

        if 1==1
            bd=[25 65];
                xlim(bd)
set(gca, 'XTick',(bd(1):5:bd(2))+3)
tmp=cvAccuracyAll{1}.timeWindow((bd(1):5:bd(2)));
tmp(1)=-1;
        set(gca, 'XTickLabel',tmp)
        end
        ylabel(sprintf('Generalization Accuracy '))
        xlabel('Imagery Training Window')
plt.vline(38,{'k--'})
title(sprintf('Generalization from \n Imagery to Experience'))
        plt.SaveFigure(1,'/Users/tysonaflalo/Dropbox/Lab/Papers/NS tactile/NewFigs','GenImag2Exp','PNG','PDF','SVGI')