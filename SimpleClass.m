classdef SimpleClass < handle  & Utilities.StructableHierarchy & Utilities.Structable
    
    % Transform nsX file into wavelet features sampled at the task
    % intervals.
    % this uses the standard {nspIDX}{FeatureIDX} format for the output.
    
    
    properties (GetAccess = 'public', SetAccess = 'public')
        % parameters
        dec=Predictor.FWClassifier(@Analyze.FaceScratch.ClassifierConfig);
        overwrite=0
        CondType = {'R','Imag'};
        Type='Classification'; % Correlatin
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        % parameters
        cvAccuracy % CAR pca coefs
        timewindow % CAR pca score for CARinds
    end
    
    
    methods
        
        function obj=SimpleClass(varargin)
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
            
            
            opts.cvOptions.ValidationType='ClassicCrossValidation';
            opts.cvOptions.NReps=10;
            opts.cvOptions.NFolds=10;
            
            
            obj.dec.Train('TrainingData',...
                {FRTrue,[1:length(LabelTrue)],LabelTrue},'CrossValidate','cvOptions',opts.cvOptions);
            
            out.BootCI=obj.dec.getBootCI;
            
            out.cvAccuracy=obj.dec.Results.cvAccuracy;
            out.confusionMat=obj.dec.Results.ConfusionMat;
            
            save(fullfile(Basedir,SaveID),'-struct','out')
            
        end
        
        function plotAcc(obj,varargin)
            
            [varargin,basename] = Utilities.ProcVarargin(varargin,'basename',[]);
            
            plt.fig('units','inches','width',4,'height',4,'font','Helvetica','fontsize',16)
            obj.dec.plotConfusion
            axis image
            colormap(cool)
            if ~isempty(basename)
                SaveID=sprintf('%s_FaceScratch', basename);
                Basedir=fullfile(env.get('result'),'FaceScratch','ImagExp','Figs');
                
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