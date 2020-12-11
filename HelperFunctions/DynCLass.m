classdef DynCLass < handle  & Utilities.StructableHierarchy & Utilities.Structable
    
    % Transform nsX file into wavelet features sampled at the task
    % intervals.
    % this uses the standard {nspIDX}{FeatureIDX} format for the output.
    
    
    properties (GetAccess = 'public', SetAccess = 'public')
        % parameters
        dec=Predictor.FWClassifier(@Analyze.FaceScratch.ClassifierConfig);
        overwrite=0
        CondType = {'R','L','Imag'};
        Type='Classification'; % Correlatin
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        % parameters
        cvAccuracy % CAR pca coefs
        timewindow % CAR pca score for CARinds                
    end
    
    
    methods
        
        function obj=DynCLass(varargin)
            [varargin]=Utilities.argobjprop(obj,varargin);
            Utilities.argempty(varargin)
        end
        
        
        function [out,SaveID]=fit_transform(obj,basename,FRTrue,LabelTrue,timeWindow,varargin)
            % transform braodband data to features binned at the same time
            % intervals as the task file.

            
            SaveID=sprintf('%s_FaceScratch', basename);
            Basedir=fullfile(env.get('result'),'FaceScratch','ImagExp');
            if exist(fullfile(Basedir,[SaveID '.mat'])) && ~obj.overwrite
                disp('file exists - loading')
                tmp=load(fullfile(Basedir,SaveID))                
                out.cvAccuracy=tmp.cvAccuracy;
                out.timeWindow=tmp.timeWindow;
                return
            end
            
            switch lower(obj.Type)
                case 'classification'
                    cvAccuracy=DPA.DynamicClassification(FRTrue,LabelTrue,timeWindow,obj.dec);
                case 'correlation'
                    cvAccuracy=DPA.DynamicGroupCorrelation(FRTrue,LabelTrue,timeWindow);
                otherwise
                    error('Unsupported method')
            end
         
            out.cvAccuracy=cvAccuracy;
            out.timeWindow=timeWindow;
            
            
            save(fullfile(Basedir,SaveID),'-struct','out')
            
        end
        
        function plotAcc(obj,cvAccuracy,timeWindow,varargin)
            
                        [varargin,basename] = Utilities.ProcVarargin(varargin,'basename',[]);
            

                        
            args.BrainArea={''};%
            args.Conditions=obj.CondType;
            args.timeWindow=timeWindow;
            args.scale=[33.3 90];
            args.scale=[-1 1];
            args.scaleType='Values'; %'Values' 'Percentile'
            
            args.scale=[5 95];
            args.scaleType='Percentile'; %'Values' 'Percentile'
            args.commonScale=1;
            args.FigSize(1)=8;
            args.FigSize(2)=10;
            
            DPA.PlotDynamicAccuracy(cvAccuracy,args)
            colormap(cool)
            colormap(jet)
            if ~isempty(basename)
            SaveID=sprintf('%s_FaceScratch', basename);
            Basedir=fullfile(env.get('result'),'FaceScratch','ImagExp','Figs');
            colormap(jet)
            plt.SaveFigure(1,Basedir,SaveID,'PDF')
            end
        end
        
        
        
        function Z=transform(obj,Z)
        end
    end
    
    methods(Static)
        function Raw2Wavelet()
            
        end
    end
    
end