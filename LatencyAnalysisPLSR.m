% ImageryLatencyLR(varargin)

%% Latency, full population
%%
Dates = {'20180709','20180711','20180716','20180718','20180723','20180730','20180806','20180815'}; % new, for imagery
AllTrialData=Analyze.LoadConvertedData('FaceScratch',Dates)

%%
    timeWindow=[-.15:.002:.25];
%     timeWindow=[0:.002:.15];

    DF=DataFilters.FaceScr_Latency2('GoTimes',timeWindow,'kernelWidth',.002,'Effector',{'Shoulder','Cheek'},'SmoothType','BC');
DC=Analyze.FaceScratch.LatencyDecode('overwrite',1,'CondType', {'L','R'});

plsrt_=dimred.plsr_ae('SmoothParams',[.8 .002],'n_components',5,'shouldTrackResiduals',0);

LA=AU.LatencyAnalysis('dimred',plsrt_);

prefix='LatencyPop'; % classifiation .5 kernel
% prefix='corr'; % correlation

%
[FRTrue,LabelTrue,timeWindow]=DF.fit_transform(AllTrialData);
disp('done')
%
% LA.fit_DiscrimTransform('LA',FRTrue,LabelTrue,timeWindow)

%
LA.fit_transform('LA',FRTrue,timeWindow)
%
LA.PlotData()
LA.PlotPWFit()
% LA.BootStrapPWFit()
