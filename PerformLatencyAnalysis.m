
%%
load('/Users/tysonaflalo/Dropbox/GitCode/Tyson/ForDistribution/eLifeTactile/LatencyData')
plsrt_=dimred.plsr_ae('SmoothParams',[.8 .002],'n_components',5,'shouldTrackResiduals',0);
LA=AU.LatencyAnalysis('dimred',plsrt_);


LA.fit_transform('LA',LatencyData.NeuralData,LatencyData.Time)
%%
LA.PlotData()
LA.PlotPWFit()

