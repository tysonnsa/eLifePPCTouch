classdef plsr_ae < handle  & Utilities.StructableHierarchy & Utilities.Structable
    
    % PLSR_dimred. Linear autoencoder using PLSR regression.
    % For time series, consider the smooth option. Finds mapping that
    % reconstrcts the smoothed originial data. Reasonable if we expect
    % autocorrelation.
    %
    % INPUT
    % Z: Num Observations X num Features
    
    
    
    properties (GetAccess = 'public', SetAccess = 'public')
        % parameters
        n_components = 15;
        SmoothParams = [];
        resTrackingSmoothParam=.99;
        silenceOutliers = false;
        shouldTrackResiduals = true;
        plotResiduals=false;
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        % attributes
        mu; % mean firing
        components; % Principal axes in feature space.
        beta; % mapping to original feature space (for tracking residuals)
        
        resTrack
        FitInfo; % Information about the dim reduction.
    end
    
    
    methods
        
        function this=plsr_ae(varargin)
            [varargin]=Utilities.argobjprop(this,varargin);
            Utilities.argempty(varargin)  
        end
        
        
        function XS=fit_transform(this,Z,~)
            
            % smoothing option smooths the
            if ~isempty(this.SmoothParams)
                [sZ]=Smooth.SmoothPopulation(Z,'mj',this.SmoothParams(1),this.SmoothParams(2));
            else
                sZ=Z;
            end
            [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(Z,sZ,this.n_components);
            
            
            this.mu=mean(Z,1);
            this.components=stats.W;
            this.beta=BETA;
            
            L=(Z-this.mu)*this.components;
            sigma=std(L);
            this.components=stats.W./sigma; % normalize to unit variance
            
            Zfit = [ones(size(Z,1),1) Z]*BETA;
            
            this.FitInfo.XS=XS;
            this.FitInfo.BETA=BETA;
            this.FitInfo.Zfit=Zfit;
            
            [a,b]=corr(Z,Zfit);
            this.FitInfo.R2=diag(a).^2;
            
            this.FitInfo.PCTVAR=PCTVAR;
            this.FitInfo.residuals=Z-Zfit;
            
            this.FitInfo.residualsSmooth=Smooth.SmoothPopulation(this.FitInfo.residuals,'exp',[10 this.resTrackingSmoothParam] ,0.03);
            this.FitInfo.resMu=mean(this.FitInfo.residualsSmooth);
            this.FitInfo.resSigma=std(this.FitInfo.residualsSmooth);
            this.FitInfo.normResiduals=(this.FitInfo.residuals-this.FitInfo.resMu)./this.FitInfo.resSigma;
            %%
            %             figure; plot(cumsum(this.FitInfo.residuals))
            % over longer time horizons, if a channel shows consistent bias, we may
            % want to minimize its impact.
            if this.plotResiduals
                plt.fig('units','inches','width',15,'height',4,'font','Helvetica','fontsize',14);
                plot((1:size(Zfit,1))*0.03,Smooth.SmoothPopulation(this.FitInfo.residuals,'exp',[10 .999] ,0.03))
                axis tight
            end
            %%
        end
        
        function [L,res]=transform(this,Z,varargin)
            [varargin,plot_residuals] = Utilities.ProcVarargin(varargin, 'plot_residuals', false);
            
            L=(Z-this.mu)*this.components;
            
            res=this.trackResidual(Z);
            if this.shouldTrackResiduals
                if size(Z,1)==1
                    this.resTrack=this.resTrackingSmoothParam*this.resTrack+...
                        (1-this.resTrackingSmoothParam)*res;
                    
                else
                    this.resTrack(1,:)=zeros(size(res(1,:)));
                    for i=2:size(Z,1)
                        this.resTrack(i,:)=this.resTrackingSmoothParam*this.resTrack(i-1,:)+...
                            (1-this.resTrackingSmoothParam)*res(i,:);
                    end
                end
                
                if this.silenceOutliers
                    GoodVals=this.resTrack<this.FitInfo.resSigma*2;
                    Z=Z.*GoodVals;
                    L=(Z-this.mu)*this.components;
                end
            end
            
            if plot_residuals
                X=cumsum(res);
                figure; plot(X)
                for i=1:size(X,2)
                    text(size(X,1),X(end,i),num2str(i));
                end
            end
        end
        
        function res=trackResidual(this,Z)
            % Compute the residual between the reconstructed feature and the
            % actual feature
            Zfit = [ones(size(Z,1),1) Z]*this.beta;
            res=Z-Zfit;
        end
    end
end