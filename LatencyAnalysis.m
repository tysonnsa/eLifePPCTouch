classdef LatencyAnalysis < handle  & Utilities.StructableHierarchy & Utilities.Structable
    
    % Transform nsX file into wavelet features sampled at the task
    % intervals.
    % this uses the standard {nspIDX}{FeatureIDX} format for the output.
    
    
    properties (GetAccess = 'public', SetAccess = 'public')
        % parameters
        dimred=[]
        overwrite=1;
        MaxNumChanges=2;
        stati='linear';
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        % parameters
        Lat;
        Mu;
        LatS;
        MuLat1;
        MuLat2;
        CILat1
        dPrime
        timeWindow % CAR pca score for CARinds
    end
    
    
    methods
        
        function obj=LatencyAnalysis(varargin)
            [varargin]=Utilities.argobjprop(obj,varargin);
            Utilities.argempty(varargin)
        end
        
        
        function [out,SaveID]=fit_transform(obj,basename,FRTrue,timeWindow,varargin)
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
            
            
            for i=1:size(FRTrue,1)
                TMP= FRTrue(i,:);
                TMP=cellfun(@(x)cat(1,x{:}),TMP,'UniformOutput',false);
                missingVals=~any(isnan(TMP{1}),1);
                TMP=cellfun(@(x)(x(:,missingVals)),TMP,'UniformOutput',false);
                FROut{i}=TMP;
                
                
                %replace nans with means
                
                
                A=cat(1,TMP{:});
                
                obj.dimred.fit_transform(A);
                
                for j=1:length(TMP)
                    [L{j}]=transform(obj.dimred,TMP{j});
                end
                
                for jj=1:size(L{1},2)
                    disp(jj)
                    tmp=cell2mat(cellfun(@(x)nanmean(x(:,jj),1),L,'UniformOutput',false));                                        
                    
                    if mean(tmp(1:floor(length(tmp)/2))) > mean(tmp(ceil(length(tmp)/2):end))
                            for kk=1:length(L)
                                L{kk}(:,jj)=-L{kk}(:,jj);
                            end
%                         L=cellfun(@(x)x(:,jj)=-(x(:,jj)),L,'UniformOutput',false);
                    end                    
                end
                
                Lat{i}=L;
                Mu{i}=cell2mat(cellfun(@(x)nanmean(x(:)),FROut{i},'UniformOutput',false));
                
                LatS{i}=cell2mat(cellfun(@(x)(x(:,1)),L,'UniformOutput',false));
                MuLat1{i}=cell2mat(cellfun(@(x)nanmean(x(:,1),1),Lat{i},'UniformOutput',false));
                
                CILat1{i}=cell2mat(cellfun(@(x)bootci(2000,@nanmean,x(:,1),1),Lat{i},'UniformOutput',false));
                
                MuLat2{i}=cell2mat(cellfun(@(x)nanmean(x(:,2),1),Lat{i},'UniformOutput',false));
                MuLat3{i}=cell2mat(cellfun(@(x)nanmean(x(:,3),1),Lat{i},'UniformOutput',false));
            end
            
            obj.Lat=Lat;
            obj.Mu=Mu;
            obj.LatS=LatS;
            obj.MuLat1=MuLat1;
            obj.MuLat2=MuLat2;
            obj.CILat1=CILat1;
            obj.timeWindow=timeWindow;
            
            out.Lat=Lat;
            out.Mu=Mu;
            out.LatS=LatS;
            out.MuLat1=MuLat1;
            out.MuLat2=MuLat2;
            out.CILat1=CILat1;
            out.timeWindow=timeWindow;
            
            save(fullfile(Basedir,SaveID),'-struct','out')
            
        end
        
         function [out,SaveID]=fit_DiscrimTransform(obj,basename,FRTrue,Labls,timeWindow,varargin)
            % transform braodband data to features binned at the same time
            % intervals as the task file.
            CondIDX=[1 2]
            
            SaveID=sprintf('%s_FaceScratch', basename);
            Basedir=fullfile(env.get('result'),'FaceScratch','Discrim');
            if exist(fullfile(Basedir,[SaveID '.mat'])) && ~obj.overwrite
                disp('file exists - loading')
                tmp=load(fullfile(Basedir,SaveID))
                out.cvAccuracy=tmp.cvAccuracy;
                out.timeWindow=tmp.timeWindow;
                return
            end
            
            for i=1:size(FRTrue,1)
                TMP= FRTrue(i,:);
                
                TMP=cellfun(@(x)cat(1,x{:}),TMP,'UniformOutput',false);
                missingVals=~any(isnan(TMP{1}),1);
                TMP=cellfun(@(x)(x(:,missingVals)),TMP,'UniformOutput',false);
                FROut{i}=TMP;
                
                
                %replace nans with means
                 for j=1:length(TMP)
                     j
                     N=cellfun(@(x)size(x,1),FRTrue{i});
                     A=[zeros(N(1),1); ones(N(2),1)];
                     
%                      Results=Analyze.FitPLSLM(TMP{j},A,'Ncomp',5)
                     Results=Analyze.FitPLSLM([obj.Lat{i}{j}],A,'Ncomp',5);
                     dP(j)=Results.dP;
                     dPp(j)=Results.dPp;
                end
                
              DDPP{i}=dP;
            end
            
            
            plot(obj.timeWindow,DDPP{1})
            
            obj.dPrime=DDPP;
          
            
          
            out.dPrime=obj.dPrime;
            
%             save(fullfile(Basedir,SaveID),'-struct','out')
            
        end
        
        function PlotData(obj,cvAccuracy,timeWindow,varargin)
            
            plt.fig('units','inches','width',10,'height',7,'font','Helvetica','fontsize',16)
            clrs=lines(length(obj.CILat1));
            
            pnl = panel();  pnl.margin=20; pnl.pack(2,2); pnl.fontname='arial';pnl.fontsize=12;
            pnl(1,1).select(); hold on
            for i=1:length(obj.CILat1)
                plot(obj.timeWindow,obj.Mu{i},'color',clrs(i,:),'linewidth',2)
            end
            xlabel('Time')
            ylabel('Mean Rate (Hz)')
            
            
            pnl(1,2).select(); hold on
            % plot(timeWindow,-MuLat1{1},'r')
            % plot(timeWindow,MuLat1{2},'k')
            
            for i=1:length(obj.CILat1)
                CI=obj.CILat1{i};
                CI=CI-obj.MuLat1{i};
                Utilities.boundedline(obj.timeWindow,obj.MuLat1{i},mean(abs(CI'),2),'alpha','cmap',clrs(i,:),'linewidth',2);
                %     Utilities.boundedline(timeWindow,MuLat1{1},abs(CI'),'alpha','cmap',clrs(1,:),'linewidth',2);
                
                preIDX=find(obj.timeWindow<0);
                preVal{i}=obj.MuLat1{i}(:,preIDX);
                
                V95Mu(i)=prctile(preVal{i}(:),95)
                
                plt.hline(V95Mu(i),{'color',clrs(i,:)});
                
            end
            
            xlabel('Time')
            ylabel('Latent Proj')
            axis tight
            pnl(2,1).select(); hold on
            
            for i=1:length(obj.CILat1)
                CI=obj.CILat1{i};
                CI=CI-obj.MuLat1{i};
                Utilities.boundedline(obj.timeWindow,obj.MuLat1{i},mean(abs(CI'),2),'alpha','cmap',clrs(i,:),'linewidth',2);
                %     Utilities.boundedline(timeWindow,MuLat1{1},abs(CI'),'alpha','cmap',clrs(1,:),'linewidth',2);
                
            end
            
            xlim([0 .2])
            
            pnl(2,2).select(); hold on
            % plot(timeWindow,-MuLat1{1},'r')
            % plot(timeWindow,MuLat1{2},'k')
            
            for i=1:length(obj.CILat1)
                plot(obj.timeWindow,obj.LatS{i}','color',clrs(i,:))
                preIDX=find(obj.timeWindow<0);
                preVal{i}=obj.LatS{i}(:,preIDX);
                
                V95(i)=prctile(preVal{i}(:),95)
                
                plt.hline(V95(i),{'color',clrs(i,:)});
            end
            axis tight
            
            xlabel('Time')
            ylabel('Latent Proj')
            
            pnl(1,2).select(); hold on
            % plot(timeWindow,-MuLat1{1},'r')
            % plot(timeWindow,MuLat1{2},'k')
            
            
            
            %
            %             if ~isempty(basename)
            %                 SaveID=sprintf('%s_FaceScratch', basename);
            %                 Basedir=fullfile(env.get('result'),'FaceScratch','ImagExp','Figs');
            %                 colormap(jet)
            %                 plt.SaveFigure(1,Basedir,SaveID,'PDF')
            %             end
            
        end
        
        
        function BootStrapPWFit(obj,varargin)
            
%             plt.fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16)
            
            clrs=lines(length(obj.CILat1));
            
            %piecewise linear
            tmp=find(obj.timeWindow>.1);
            P0=[0     0     0    tmp(1)    30     0];
            lb = -inf*[1 1 1 1 1 1];
            ub  = inf*[1 1 1 1 1 1];
            plusfun = @(x) max(x,0);
            
            model = @(P,x) P(1) + P(6)*x + (P(6)-P(2))*plusfun(P(4)-x) + (P(3)-P(6))*plusfun(x-(P(4)+P(5)));
            options = optimoptions('lsqcurvefit','Display','off');
            
            
            for i=1:length(obj.CILat1)
                for j=1:500
                    RandSamp=obj.LatS{i};
                    RandSamp=mean(RandSamp(randi(size(RandSamp,1),size(RandSamp,1),1),:),1);
%                     tmp=my_findchangepts( RandSamp,'Statistic',obj.stati,'MaxNumChanges',obj.MaxNumChanges,'MinDistance',10);
%                     LatSamp{i}(j)=obj.timeWindow(tmp(1));
                    
                    
                    
                    Pfit = lsqcurvefit(model,P0,1:length(RandSamp),RandSamp,lb,ub,options);
                    
                    modelpred = model(Pfit,sort(1:length(RandSamp)));
%                     plot(1:length(RandSamp),RandSamp,'o',sort(1:length(RandSamp)),modelpred,'r-')
                    
                    
                    tmp=find(modelpred>prctile( RandSamp(obj.timeWindow<0),95));
%                     [obj.timeWindow(tmp(1)) obj.timeWindow(round(Pfit(4))) obj.timeWindow(tmp(1))]
                    
                    LatSampA{i}(j)=obj.timeWindow(tmp(1));
                    LatSampB{i}(j)=obj.timeWindow(round(Pfit(4)));
                end
                
                LatCIA{i}=prctile(LatSampA{i},[25 75]);
                LatCIB{i}=prctile(LatSampB{i},[25 75]);
                
                disp(LatCIA{i})
                disp(LatCIB{i})
                a=histogram(LatSampA{i}*1000,1000*(0.01:.005:.1),'Normalization','probability')
                
                
                text(10,0+i/10,sprintf('%0.0f   ',1000*LatCIA{i}))
            end
            
            xlabel('Time of significant modulation (ms)')
            ylabel('Probability')
        end
        
        function PlotPWFit(obj,cvAccuracy,timeWindow,varargin)
            
            
            plt.fig('units','inches','width',6,'height',3,'font','Helvetica','fontsize',16)
            
            pnl = panel();  pnl.margin=15; pnl.pack('h',{2/3 1/3}); pnl.fontname='arial';pnl.fontsize=12;
            
            pnl(1).select()
            
            clrs=lines(length(obj.CILat1));
                        tmp=find(obj.timeWindow>.08);

             P0=[0     0     0    tmp(1)    30     0];
            lb = -inf*[1 1 1 1 1 1];
            ub  = inf*[1 1 1 1 1 1];
            plusfun = @(x) max(x,0);
            
            model = @(P,x) P(1) + P(6)*x + (P(6)-P(2))*plusfun(P(4)-x) + (P(3)-P(6))*plusfun(x-(P(4)+P(5)));
            options = optimoptions('lsqcurvefit','Display','off');
            
            for i=1:length(obj.CILat1)
                
                
                CI=obj.CILat1{i};
                CI=CI-obj.MuLat1{i};
%                 my_findchangepts( obj.MuLat1{i},'Statistic',obj.stati,'MaxNumChanges',obj.MaxNumChanges,'color',clrs(i,:));
                tmp=findchangepts( obj.MuLat1{1},'Statistic',obj.stati,'MaxNumChanges',obj.MaxNumChanges);
                Utilities.boundedline(1:length(obj.timeWindow),obj.MuLat1{i},mean(abs(CI'),2),'alpha','cmap',clrs(i,:),'linewidth',2);
                %     Utilities.boundedline(timeWindow,MuLat1{1},abs(CI'),'alpha','cmap',clrs(1,:),'linewidth',2);
                 Pfit = lsqcurvefit(model,P0,1:length( obj.MuLat1{i}), obj.MuLat1{i},lb,ub,options);
                    
                    modelpred = model(Pfit,sort(1:length( obj.MuLat1{i})));
                    plot(sort(1:length( obj.MuLat1{i})),modelpred,'color',clrs(i,:),'linewidth',2)
                    
                    
                      tmp=find(obj.timeWindow<0);
                    tmp=find(modelpred>prctile( obj.MuLat1{i}(obj.timeWindow<0),95));
                    plt.vline([(tmp(1))],{'--','color',clrs(i,:),'linewidth',2})
% plt.hline(prctile( obj.MuLat1{i}(obj.timeWindow<0),95),{'--','color',clrs(i,:),'linewidth',2})
                      
%                     plt.vline([Pfit(4) Pfit(4)+Pfit(5)],{'color',clrs(i,:)},'linewidth',2)
            end
            %%
            set(gca,'XTick',1:25:length(obj.timeWindow))
            
            obj.timeWindow( cellfun(@(x)str2num(x),get(gca,'XTickLabel')))
            
%             set(gca,'XTickLabel',obj.timeWindow( 1:25:length(obj.timeWindow)))
            
            set(gca,'XTickLabel',1000*obj.timeWindow( cellfun(@(x)str2num(x),get(gca,'XTickLabel'))))
            axis tight
            xlabel('Time relative to probe contact (ms)')
            ylabel('First Latent Dim')
            
            pnl(2).select(); hold on
            obj.BootStrapPWFit()
            %%
            %
            % pnl = panel();  pnl.margin=20; pnl.pack('h',length(obj.CILat1)); pnl.fontname='arial';pnl.fontsize=12;
            %
            % for i=1:length(obj.CILat1)
            % pnl(i).select(); hold on
            %
            % MaxNumChanges=2;
            % stati='linear'
            % CI=obj.CILat1{i};
            %     CI=CI-obj.MuLat1{i};
            %  my_findchangepts( obj.MuLat1{i},'Statistic',stati,'MaxNumChanges',MaxNumChanges);
            % tmp=findchangepts( obj.MuLat1{1},'Statistic',stati,'MaxNumChanges',MaxNumChanges);
            % Utilities.boundedline(1:length(obj.timeWindow),obj.MuLat1{i},mean(abs(CI'),2),'alpha','cmap',clrs(1,:),'linewidth',2);
            % %     Utilities.boundedline(timeWindow,MuLat1{1},abs(CI'),'alpha','cmap',clrs(1,:),'linewidth',2);
            %
            % end
            %
            
            %             if ~isempty(basename)
            %                 SaveID=sprintf('%s_FaceScratch', basename);
            %                 Basedir=fullfile(env.get('result'),'FaceScratch','ImagExp','Figs');
            %                 colormap(jet)
            %                 plt.SaveFigure(1,Basedir,SaveID,'PDF')
            %             end
        end
        
        
        
        
        
        
         function PlotDiscrimLat(obj,FRTrue,timeWindow,varargin)
            
            plt.fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16)
            clrs=lines(length(obj.CILat1));
                        tmp=find(obj.timeWindow>.08);

             P0=[0     0     0    tmp(1)    30     0];
            lb = -inf*[1 1 1 1 1 1];
            ub  = inf*[1 1 1 1 1 1];
            plusfun = @(x) max(x,0);
            
            model = @(P,x) P(1) + P(6)*x + (P(6)-P(2))*plusfun(P(4)-x) + (P(3)-P(6))*plusfun(x-(P(4)+P(5)));
            options = optimoptions('lsqcurvefit','Display','off');
            
            for i=1:length(obj.CILat1)
                
                
                CI=obj.CILat1{i};
                CI=CI-obj.MuLat1{i};
%                 my_findchangepts( obj.MuLat1{i},'Statistic',obj.stati,'MaxNumChanges',obj.MaxNumChanges,'color',clrs(i,:));
                tmp=findchangepts( obj.MuLat1{1},'Statistic',obj.stati,'MaxNumChanges',obj.MaxNumChanges);
                Utilities.boundedline(1:length(obj.timeWindow),obj.MuLat1{i},mean(abs(CI'),2),'alpha','cmap',clrs(i,:),'linewidth',2);
                %     Utilities.boundedline(timeWindow,MuLat1{1},abs(CI'),'alpha','cmap',clrs(1,:),'linewidth',2);
                 Pfit = lsqcurvefit(model,P0,1:length( obj.MuLat1{i}), obj.MuLat1{i},lb,ub,options);
                    
                    modelpred = model(Pfit,sort(1:length( obj.MuLat1{i})));
                    plot(sort(1:length( obj.MuLat1{i})),modelpred,'color',clrs(i,:),'linewidth',2)
                    plt.vline([Pfit(4) Pfit(4)+Pfit(5)],{'color',clrs(i,:)},'linewidth',2)
            end
            %%
            set(gca,'XTick',1:25:length(obj.timeWindow))
            
            obj.timeWindow( cellfun(@(x)str2num(x),get(gca,'XTickLabel')))
            set(gca,'XTickLabel',obj.timeWindow( cellfun(@(x)str2num(x),get(gca,'XTickLabel'))))
            
            %             if ~isempty(basename)
            %                 SaveID=sprintf('%s_FaceScratch', basename);
            %                 Basedir=fullfile(env.get('result'),'FaceScratch','ImagExp','Figs');
            %                 colormap(jet)
            %                 plt.SaveFigure(1,Basedir,SaveID,'PDF')
            %             end
        end
        
        
        
        
    end
    
end