%%Tuning

%%Notes To plot Figure 1 Tuning

% step 1. put a debugging red dot at line 161 in DataPlot.SUP_BasicTuning; this is
% the line that begins pnl(1,2).select();

%step 2. open testdPrimeTouch.m and run lines 3,4,7,8,9,11;

%step 3. next run the lines that follow this:


plt.fig('units','inches','width',12,'height',5,'font','Arial','fontsize',9);
% Basic Percent of units responsive to a given action
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1).select();

NC=13;

ActionTuned=ActionPvals<PThresh;
PercentTuned=sum(ActionTuned,1)/size(ActionTuned,1);
CI=bootci(2000,@sum,double(ActionTuned))/size(ActionTuned,1);
h = plt.barpatch([1:13],{PercentTuned*100,CI*100},'patchbar', 0, 'yl','Responsiveness (%)','fontsize',6,'t','Action Responsive',...
'fontsize',6,'barname',Labels);

clr=lines(8);

set(h.patch([3 5 7 9 11]),'facecolor',clr(1,:))
set(h.patch([4 6 8 10 12]),'facecolor',clr(2,:))
set(h.patch(13),'facecolor',[0 0 0])

% h = plt.barpatch([],{PercentTuned*100,CI*100},'barcmap','lines','patchbar',0,'yl','Responsiveness (%)','fontsize',6,'t','Action Responsive',...
%     'barname',Labels);
axis tight; xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold');
ax1 = gca;
ax1.FontWeight = 'bold';
set(ax1,'box','off');
ax1.XRuler.Axle.LineStyle = 'none';

%

% pnl(1,2).select(); gca; hold on;
% xlim([.25 NC+.25]);
% xtickangle(45); set(gca, 'fontweight','bold');
% % plt.barpatch(1:10,abs(dp(:,1:10)),'yl','d''','fontsize',6','t','Depth of modulation','barcmap','lines')
% 
% dp=cat(1,PlotData{1}.dPrime{:});
% pVals=cat(2,PlotData{1}.Pvals{:})';

%

pnl(1,2).select(); gca; hold on;
xlim([.25 NC+.25]);
xtickangle(45); set(gca, 'fontweight','bold');

h = plt.barpatch([1:10],abs(dp(:,1:10)),'yl','DI','fontsize',6,'t','Discriminability Index',...
'fontsize',6);

clr=lines(8);
set(h.patch([3 5 7 9 ]),'facecolor',clr(1,:))
set(h.patch([4 6 8 10 ]),'facecolor',clr(2,:))

ylim([1.2 2.2])

axis tight
ax1 = gca;
ax1.FontWeight = 'bold';
set(ax1,'box','off');
ax1.XRuler.Axle.LineStyle = 'none';
labelmat = {Labels{1},Labels{2},Labels{3},Labels{4},Labels{5},Labels{6},Labels{7},Labels{8},Labels{9},Labels{10}};
xticklabels(labelmat)
hold on
plot([3 4],nanmean(abs(dp(:,[3 4])),1),'k.-')
plot([5 6],nanmean(abs(dp(:,[5 6])),1),'k.-')
plot([7 8],nanmean(abs(dp(:,[7 8])),1),'k.-')
plot([9 10],nanmean(abs(dp(:,[9 10])),1),'k.-')

xlim([.25 NC+.25]);
ylim([1.2 2.2])


%% New figures
% Histgrams of left and right; number of fields to which each unit responds

plt.fig('units','inches','width',12,'height',5,'font','Arial','fontsize',9);
% Basic Percent of units responsive to a given action
pnl = panel();  pnl.margin=15; pnl.pack(2,6);
pnl.fontsize=11; pnl.fontname='Arial';

clr=lines(8);

% left side (ipsilateral)

pnl(1,1).select();

tempMatleft = ActionTuned(:,[1 2 3 5 7 9]);
tempMatleftsum = sum(tempMatleft,2);
h = histcounts(tempMatleftsum(tempMatleftsum~=0));

hfig = bar([1:6],h./sum(h)*100);

axis image square; axis tight;
ylim([0 40]);
xlabel('Number of Fields', 'fontweight', 'bold','fontsize',6);
ylabel('Neurons (%)', 'fontweight', 'bold','fontsize',6);
title('Left, Ipsilateral','fontweight', 'bold','fontsize',6);

set(gca, 'fontweight', 'bold','fontsize',6);
set(gca, 'xticklabel', [1:1:6], 'xtick', [1:1:6])
set(hfig, 'FaceColor',clr(1,:), 'FaceAlpha',.75);
set(hfig, 'BarWidth', 0.7)

% right side (contralateral)


pnl(1,2).select();

tempMatright = ActionTuned(:,[1 2 4 6 8 10]);
tempMatrightsum = sum(tempMatright,2);
h = histcounts(tempMatrightsum(tempMatrightsum~=0));

hfig = bar([1:6],h./sum(h)*100);

axis image square; axis tight;
ylim([0 40]);
xlabel('Number of Fields', 'fontweight', 'bold','fontsize',6);
ylabel('Neurons (%)', 'fontweight', 'bold','fontsize',6);
title('Right, Contralateral','fontweight', 'bold','fontsize',6);

set(gca, 'fontweight', 'bold','fontsize',6);
set(gca, 'xticklabel', [1:1:6], 'xtick', [1:1:6])
set(hfig, 'FaceColor', clr(2,:), 'FaceAlpha',.75);
set(hfig, 'BarWidth', 0.7)


%%

%%Correlation

%Notes To plot Figure 1 Correlation

% step 1. put a debugging red dot at line 70 in FaceScratch.CorrAnalysis.m
% the line that begins plt.fig

%step 2. run the lines that follow :

NC =13;

plt.fig('units','inches','width',12,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1).select();

LabelsAbbrev = {'Forehead','Vertex','L Back','R Back','L Cheek','R Cheek','L Shoulder','R Shoulder','L Neck','R Neck','L Arm','R Arm','Null'};

imagesc(rho);
axis image;

title('Correlation between conditions','fontsize',6);
set(gca,'XTick',1:NConds,'XTickLabel',LabelsAbbrev);
set(gca,'YTick',1:NConds,'YTickLabel',LabelsAbbrev);
set(gca,'YDir','normal')
colorbar;
colormap(jet);

axis square;
axis tight;
xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',6');

%
pnl(1,2).select(); ax1 = gca; hold on;

DiffVarNames = {'Laterality','Stimulation Site'};
pairs{1} = [3 4; 5 6; 7 8; 9 10];
pairs{2} = [3 5; 3 7; 3 9; 5 7; 5 9; 7 9; 4 6; 4 8; 4 10; 6 8; 6 10; 8 10];
for i = 1:length(DiffVarNames)
    corrVals{i} = [];
    for j = 1:length(pairs{i})
        corrVals{i} = [corrVals{i} rho(pairs{i}(j,1),pairs{i}(j,2))];
    end
end


for i = 1:length(corrVals)
    [mu(i),lb(i),ub(i)] = MeanAndCI(corrVals{i});
end

data = {mu'; [lb; ub]};
groupidx = 1:2;
%     xl = 'Differing Condition';
yl = 'Correlation';
t = 'Correlation by Condition';
clr = lines(7);
h  = plt.barpatchGroup(data, 'groupidx', groupidx, ...
    'yl', yl,             ...
    'fontsize',6,         ...
    'groupspace', 2,      ...
    'patchbar', 0,        ...
    'barname', DiffVarNames,      ...
    'barcmap',clr([7 4],:));
%     set(h.leg,'FontSize',9);



xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',6');
  set(h.leg,'FontSize',6, 'FontWeight', 'bold');



%% Decoding stuff

%%decoding

%Notes To plot Figure 1 decoding

% step 1. put a debugging red dot at line 79 in
% Analyze.FaceScratch.DecodeAnalysis.m
% the line that begins plt.fig

%step 2. run the lines that follow :

NC =13;

plt.fig('units','inches','width',12,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1).select();

ConfMatMu = zeros(NConds,NConds);
for d = 1:length(UDates)
    ConfMatMu = ConfMatMu + ConfMat{d};
end
ConfMatMu = ConfMatMu/length(UDates)*100;
Labels = {'Forehead','Vertex','L Back','R Back','L Cheek','R Cheek','L Shoulder','R Shoulder','L Neck','R Neck','L Arm','R Arm','Null'};
DecodeLabels = Labels;

imagesc(ConfMatMu);
colorbar;
axis([0 NConds 0 NConds]+.5);
%     set(gca);
set(gca,'XTick',1:NConds,'XTickLabel',DecodeLabels,StyleArgs{:},'XTickLabelRotation',45);
%     set(gca);
set(gca,'YTick',1:NConds,'YTickLabel',DecodeLabels,StyleArgs{:});
xlabel('Predicted Labels');
ylabel('True Labels');
title([Tag ' Confusion Matrix'],StyleArgs{:});

axis square;
axis tight;
xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',6');
colormap(cool);

%% Sample Units


PlotConds = {{'Condition','BackLeft'},...
    {'Condition','BackRight'},...
    {'Condition','LeftCheek'},...
    {'Condition','RightCheek'},...
    {'Condition','LeftShoulder'},...
    {'Condition','RightShoulder'},...
    {'Condition','LeftNeck'},...
    {'Condition','RightNeck'},...
    {'Condition','LeftArm'},...
    {'Condition','RightArm'},...
    {'Condition','XX'},...
    {'Condition','XX'}};
% Phases = {'CueTarget','Delay','Go'};
% TimeWindows = [-.5 1.5; 0 1; 0 3.5];
Phases = {'Go'};
TimeWindows = [-1 3.5];
SubplotIDs = [1 2 1 2 1 2 1 2 1 2 1 2];
% size(Units)
Units = [1 19 1 1975;
    1 19 1 1983;
    1 25 1 1975
    1 30 1 1987];

Analyze.FaceScratch.SCPlotSingleUnitERAs(AllTrialData,PlotConds,Phases,TimeWindows,OutDir,'Units',Units,...
    'SubplotIds',SubplotIDs,'FigSize',[6 9]);


%% Figure 2

%Notes To plot Figure 2

%first run the general code from the
%Analyze.FaceScratch.AnalyzeFaceScratch2 including the initial stuff with
%the Dates, phase, globaltrialcriteria, etc. then run the alltrialdata
%loading section. then, skip down to the plotanalysisthroughtime within the
%end of the first (XTouchRF) function to the area listed as within side
%analysis. next, either through the config file, or by just opening the
%file:

% step 1. open the file Analyze.FaceScratch.Plot_WithinAcross and put a
% debug marker at line 153, the line that begins plt.fig. then run the
% following code:



tmp = R2Within(:,R).*R2Across{R}(:,4);
tmp(isnan(tmp)) = 0;

[aa,ii] = sort(tmp,'descend');
BestTable = table(TaskDate(ii)',Units(ii,1),Units(ii,2),Units(ii,3),Units(ii,4),R2Within(ii,R),R2Across{R}(ii,4),...
    'VariableNames',{'Date','ASF1','ASF2','ASF3','ASF4','R2Within','R2Across'});
BestTable(1:10,:)


plt.fig('units','inches','width',12,'height',5.25,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';


pnl(1,1).select();gca; hold on;
plot(R2Within(isValid,R),R2Within(isValid,L),'k.');
axis([0 1 0 1]);
axis image square
xlabel(sprintf('R^{2} %s',Sides{R}), 'fontsize',6, 'fontweight', 'bold');
ylabel(sprintf('R^{2} %s',Sides{L}), 'fontsize',6, 'fontweight', 'bold');
hold on;
set(gca, 'fontweight', 'bold', 'fontsize',6);
plot([0 1],[0 1],'r');
title('R^{2}_{Within} x R^{2}_{Within}','fontweight', 'bold','fontsize',6);

% Usign each side as the "within"
clrs = lines(length(BPLabels));
for ss = 1:2
    % Plot within vs across
    cSide = Sides{ss};
    oSide = Sides{-(ss-2)+1};
    % Plot R2withinR v R2across
    pnl(1,ss+1).select();gca; hold on;
    plot(R2Within(isValid,ss),R2Across{ss}(isValid,2),'k.');
    axis([0 1 0 1]);
    axis image square
    % xlim([0 max(R2CV(:))]);
    % ylim([0 max(R2CV(:))]);
    xlabel(sprintf('R^{2} %s',cSide), 'fontsize',6, 'fontweight', 'bold');
    ylabel(sprintf('R^{2} across (to %s)',oSide), 'fontsize',6, 'fontweight', 'bold');
    set(gca, 'fontweight', 'bold','fontsize',6);
    hold on;
    plot([0 1],[0 1],'r');
    title(sprintf('R^{2}_{Within} x R^{2}_{Across} (%s to %s)',cSide,oSide),'fontsize',6,'fontweight','bold');
end

plt.fig('units','inches','width',12,'height',5.25,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';


% Also plot specificity for each of the scatter plots
A = R2Within(isValid,R);
B = R2Within(isValid,L);
specWithin = (abs(A) - abs(B))./(abs(A)+abs(B));
pnl(2,1).select(); gca; hold on;
figtmp = histogram(specWithin,-1:.1:1);
axis image square; axis tight;
ylim([0 45]);
% histogram(specWithin,-1:.1:1,'Normalization','pdf');
xlabel('L <-> R', 'fontweight', 'bold','fontsize',6);
[b,a]=signtest(specWithin);
title(sprintf('Median = %.2f \n p = %2.1g',median(specWithin),b),'fontweight', 'bold','fontsize',6);
set(gca, 'fontweight', 'bold','fontsize',6);
% [a,b]=ttest(specWithin);
% title(sprintf('mean = %.2f, p = %.2f',mean(specWithin),b));


A = R2Within(isValid,L);
B = R2Across{L}(isValid,2);
specAcross{L} = (abs(A) - abs(B))./(abs(A)+abs(B));
pnl(2,2).select();
histogram(specAcross{L},-1:.1:1);
axis image square; axis tight;
ylim([0 45]);
% histogram(specAcross{L},-1:.1:1,'Normalization','pdf');
xlabel('R <-> L', 'fontweight', 'bold','fontsize',6);
[b,a]=signtest(specAcross{L});
title(sprintf('Median = %.2f \n p = %2.1g',median(specAcross{L}),b),'fontweight', 'bold','fontsize',6);
set(gca, 'fontweight', 'bold','fontsize',6);
% [a,b]=ttest(specAcross{L});
% title(sprintf('mean = %.2f, p = %.2f',mean(specAcross{L}),b));

A = R2Within(isValid,R);
B = R2Across{R}(isValid,2);
specAcross{R} = (abs(A) - abs(B))./(abs(A)+abs(B));
pnl(2,3).select();
histogram(specAcross{R},-1:.1:1);
axis image square; axis tight;
ylim([0 45]);
% histogram(specAcross{R},-1:.1:1,'Normalization','pdf');
xlabel('L <-> R', 'fontweight', 'bold','fontsize',6);
[b,a]=signtest(specAcross{R});
title(sprintf('Median = %.2f \n p = %2.2g',median(specAcross{R}),b),'fontweight', 'bold','fontsize',6);
set(gca, 'fontweight', 'bold','fontsize',6);

%% Figure 4 Tuning

plt.fig('units','inches','width',12,'height',5,'font','Arial','fontsize',9);
% Basic Percent of units responsive to a given action
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1).select();

NC=13;

Labels{10} = 'Null';

ActionTuned=ActionPvals<PThresh;
PercentTuned=sum(ActionTuned,1)/size(ActionTuned,1);
CI=bootci(2000,@sum,double(ActionTuned))/size(ActionTuned,1);
plt.barpatch([],{PercentTuned*100,CI*100},'patchbar',0,'yl','Responsiveness (%)','fontsize',6,'t','Action Responsive',...
    'barname',Labels);
axis tight; xlim([.25 NC+.25]); xtickangle(0); set(gca, 'fontweight','bold');
ax1 = gca;
ax1.FontWeight = 'bold';
set(ax1,'box','off');
ax1.XRuler.Axle.LineStyle = 'none';


%

dp=cat(1,PlotData{1}.dPrime{:});
pVals=cat(2,PlotData{1}.Pvals{:})';


H=Utilities.MultipleComparisonsCorrection(pVals,'method','fdr');
sum(H)
len = length(sum(H));

dp(~H)=nan;

pnl(1,2).select(); gca; hold on;
xlim([.25 NC+.25]);
xtickangle(45); set(gca, 'fontweight','bold');
plt.barpatch(1:9,abs(dp(:,1:9)),'yl','DI','fontsize',6','t','Discriminability Index')
labelmat = {'1','2','3','4','5','6','7','8','9'};
xticklabels(labelmat)

axis tight; xlim([.25 NC+.25]); xtickangle(0); set(gca, 'fontweight','bold');
ax1 = gca;
ax1.FontWeight = 'bold';
set(ax1,'box','off');
ax1.XRuler.Axle.LineStyle = 'none';
ylim([1.0 2.3])


%% Figure 4 Correlation

%%Correlation

%Notes To plot Figure 1 Correlation

% step 1. put a debugging red dot at line 70 in FaceScratch.CorrAnalysis.m
% the line that begins plt.fig

%step 2. run the lines that follow :

NC =13;

plt.fig('units','inches','width',12,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1).select();

Labels2Incl = Labels;
LabelsAbbrev = Labels;
LabelsAbbrev{10} = 'Null';

imagesc(rho);
axis image;

title('Correlation between conditions','fontsize',6);
set(gca,'XTick',1:NConds,'XTickLabel',LabelsAbbrev);
set(gca,'YTick',1:NConds,'YTickLabel',LabelsAbbrev);
set(gca,'YDir','normal')
colorbar;
colormap(jet);

axis square;
axis tight;
xlim([.25 NC+.25]); xtickangle(0); set(gca, 'fontweight','bold','fontsize',6');

%% Figure 4 Sample Units

Coef = horzcat(PlotData{1}.Coef{:})';
CoefCI = horzcat(PlotData{1}.CoefCI{:})';
tmp = diff(CoefCI);
CoefCI = tmp(1:2:end,:)/2;

QualThresh = 1;
% meeting criteria
validUnits=find(PlotData{1}.R2 > .7 & PlotData{1}.UnitQuality<=2);
% top 25
validUnits = find(PlotData{1}.UnitQuality==QualThresh);
[r2sorted,r2sorti] = sort(PlotData{1}.R2(validUnits),'descend');
validUnits = validUnits(r2sorti(1:8));

[nx]=plt.getSubPlotDimensions(length(validUnits));
nx = [2 4];
plt.fig('units','inches','width',9,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15;


id=1;
for unit=validUnits;
    plt.subaxis(nx(1),nx(2),id); gca; hold on;
    
    plt.shadedErrorBar(1:9,Coef(unit,1:9),CoefCI(unit,1:9));
    id=id+1;
    axis tight; axis square; set(gca, 'fontweight', 'bold');
    %     title(sprintf('R2=%.2f',PlotData{1}.R2(unit)));
    ylabel('Firing Rate (Hz)', 'fontweight', 'bold', 'fontsize', 8);
    xlabel('Stimulation Site', 'fontweight', 'bold', 'fontsize', 8);
    
end
suptitle(sprintf('Best quality (%d)',QualThresh));

%% Figure 5 Distribution of Correlations

figH = plt.fig('units','inches','width',10,'height',4,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';



edges = 0:.01:.4;
edges2 = -.02:.005:.09;
diffLabel = {'Touch x Loom Near', 'Touch x Loom Far', 'Loom Near x Loom Far'};

    % Plot true distributions of correlations
    pnl(1,1).select(); gca; hold on; ax1=gca;
    histogram(RelCorrs(:,1),edges, 'FaceColor',[0.3 0.3 0.7],'FaceAlpha',0.4);
    hold on;
    histogram(RelCorrs(:,2),edges,'FaceColor',[0.1 0.4 0.1880],'FaceAlpha',0.4);
%     histogram(RelCorrs(:,3),edges,'FaceColor',[0.45 0.45 0.25],'FaceAlpha',0.4);    
    title('Distribution of Correlations','fontweight', 'bold', 'fontsize', 6);
    Leg =  legend(diffLabel{1},diffLabel{2});
    if i==1
        set(legend,'Position',[0.5446 0.3774 0.0729 0.0339]);
    elseif i==2
        set(legend,'Position',[0.5446 0.8019 0.0729 0.0339]);
    end
    legend boxoff;
    
    
    xlabel('Correlation','fontweight', 'bold', 'fontsize', 6);
    axis tight;
    set(gca,'fontweight', 'bold', 'fontsize', 6);
    ax1.FontWeight = 'bold';
    xlim([0 0.35])
    ylim([0 50])




    

%% Figure 4 Receptive Field Size 

% Run the file in the Desktop called RFSizeModified. Note, if the middle
% section doesn't run (if it errors out) at first, comment out the title
% line (with the sprintf). run the third section, and then come back to the
% middle section (and remember to uncomment the sprintf line in the title).
% this is needed because the fwhm (full width half max) isn't calculated
% until the last section

%% Figure 6 Tuning Curves

plt.fig('units','inches','width',12,'height',5,'font','Arial','fontsize',9);
% Basic Percent of units responsive to a given action
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1).select();

NC=13;

ActionTuned=ActionPvals<PThresh;
PercentTuned=sum(ActionTuned,1)/size(ActionTuned,1);
CI=bootci(2000,@sum,double(ActionTuned))/size(ActionTuned,1);
h = plt.barpatch([1:3 5:7 9:11 13],{PercentTuned*100,CI*100},'patchbar', 0, 'yl','Responsiveness (%)','fontsize',6,'t','Action Responsive',...
'fontsize',6,'barname',Labels);

clr=lines(8);

set(h.patch([1 4 7]),'facecolor',clr(1,:))
set(h.patch([1 4 7]+1),'facecolor',clr(2,:))
set(h.patch([1 4 7]+2),'facecolor',clr(7,:))
set(h.patch(10),'facecolor',[0 0 0])

% h = plt.barpatch([],{PercentTuned*100,CI*100},'barcmap','lines','patchbar',0,'yl','Responsiveness (%)','fontsize',6,'t','Action Responsive',...
%     'barname',Labels);
axis tight; xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold');
ax1 = gca;
ax1.FontWeight = 'bold';
set(ax1,'box','off');
ax1.XRuler.Axle.LineStyle = 'none';
%%

%

dp=cat(1,PlotData{1}.dPrime{:});
pVals=cat(2,PlotData{1}.Pvals{:})';


H=Utilities.MultipleComparisonsCorrection(pVals,'method','fdr');
sum(H)
len = length(sum(H));

dp(~H)=nan;

pnl(1,2).select(); gca; hold on;
xlim([.25 NC+.25]);
xtickangle(45); set(gca, 'fontweight','bold');

h = plt.barpatch([1:3 5:7 11],abs(dp(:,[1:6 9])),'yl','DI','fontsize',6,'t','Discriminability Index',...
'fontsize',6);

clr=lines(8);


set(h.patch([1 4]),'facecolor',clr(1,:))
set(h.patch([2 5]),'facecolor',clr(2,:))
set(h.patch([3 6 7]),'facecolor',clr(7,:))


% plt.barpatch([1:6 9],abs(dp(:,[1:6 9])),'yl','d''','fontsize',6','t','Depth of modulation')
labelmat = {Labels{1},Labels{2},Labels{3},Labels{4},Labels{5},Labels{6},Labels{9},Labels{10}};
xticklabels(labelmat)
xtickangle(45);
axis tight; xlim([.25 NC+.25]); set(gca, 'fontweight','bold');
ax1 = gca;
ax1.FontWeight = 'bold';
set(ax1,'box','off');
ax1.XRuler.Axle.LineStyle = 'none';
ylim([1.0 2.0])
hold on
% plot([1 2],nanmean(abs(dp(:,[1 2])),1),'k.-')
% plot([5 6],nanmean(abs(dp(:,[4 5])),1),'k.-')

%% Figure 6 Correlation

%%Correlation

%Notes To plot Figure 1 Correlation

% step 1. put a debugging red dot at line 70 in FaceScratch.CorrAnalysis.m
% the line that begins plt.fig

%step 2. run the lines that follow :

NC =13;

plt.fig('units','inches','width',12,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1).select();

Labels2Incl = Labels;
LabelsAbbrev = Labels;

NConds = 10;

imagesc(rho);
axis image;

title('Correlation between conditions','fontsize',6);
set(gca,'XTick',1:NConds,'XTickLabel',LabelsAbbrev);
set(gca,'YTick',1:NConds,'YTickLabel',LabelsAbbrev);
set(gca,'YDir','normal')
colorbar;
colormap(jet);

axis tight;
xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',6);

%% Figure 6 Distribution of Correlations

figH = plt.fig('units','inches','width',10,'height',4,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';


    
    edges = 0:.01:.4;
    edges2 = -.02:.005:.09;
    diffLabel = {'Touch Cheek x Imagine Cheek','Touch Cheek x Imagine Shoulder','Touch Shoulder x Imagine Shoulder','Touch Shoulder x Imagine Cheek','Imagine Shoulder x Imagine Cheek','Imagine Shoulder x Imagine Cheek'};
    for i = 1:2
        % Plot true distributions of correlations
        pnl(i,1).select(); gca; hold on; ax1=gca;
        histogram(RelCorrs(:,(i-1)*2+1),edges, 'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',0.4);
        hold on;
        histogram(RelCorrs(:,i*2),edges,'FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.4);
%         histogram(RelCorrs(:,i+4),edges,'FaceColor',[0.4 0.4 0.6], 'FaceAlpha',0.4);
        title('Distribution of Correlations','fontweight', 'bold', 'fontsize', 6);
       Leg =  legend(diffLabel{(i-1)*2+1},diffLabel{i*2},diffLabel{i+4});
        if i==1
            set(legend,'Position',[0.5446 0.3774 0.0729 0.0339]);
        elseif i==2
            set(legend,'Position',[0.5446 0.8019 0.0729 0.0339]);
        end
        legend boxoff; 
    
    
    xlabel('Correlation','fontweight', 'bold', 'fontsize', 6);
    axis tight;
    set(gca,'fontweight', 'bold', 'fontsize', 6);
    ax1.FontWeight = 'bold';
    ylim([0 70])
    end
               
        
%         % Plot null dist for corrdiffs
%         pnl(i,2).select(); gca; hold on; ax1=gca;
%         histogram(NullDiffs(i,:),edges2);
%         Leg = legend;
%         plt.vline(truediff(i));
%         title('NullDist');
%         xlabel(sprintf('Correlation differences (%s-%s)',diffLabel{i*2},diffLabel{(i-1)*2+1}));
%         axis tight; 
%         set(gca,'fontweight', 'bold', 'fontsize', 6);
%         ax1.FontWeight = 'bold';
 
% pnl(i,3).select()        
%         [p,df,n,chi2] = kruskalwallis2(RelCorrs(:,(i-1)*2+1),RelCorrs(:,i*2));
%         x = diff(xlim)/4 + max(xlim);
%         y = diff(ylim)/2 + min(ylim);
%         legend([x y],diffLabel{(i-1)*2+1},diffLabel{i*2})
%         text(x,y,sprintf('Kruskal-Wallis\np=%.2e\ndf=%d\nn=%d\nchi2=%.2f',p,df,n,chi2));
    

%% Figure 5 Tuning

plt.fig('units','inches','width',12,'height',5,'font','Arial','fontsize',9);
% Basic Percent of units responsive to a given action
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1).select();

NC=13;

ActionTuned=ActionPvals<PThresh;
PercentTuned=sum(ActionTuned,1)/size(ActionTuned,1);
CI=bootci(2000,@sum,double(ActionTuned))/size(ActionTuned,1);
plt.barpatch([],{PercentTuned*100,CI*100},'patchbar',0,'yl','Responsiveness (%)','fontsize',6,'t','Action Responsive',...
    'barname',Labels);
axis tight; xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold');
ax1 = gca;
ax1.FontWeight = 'bold';
set(ax1,'box','off');
ax1.XRuler.Axle.LineStyle = 'none';


%

dp=cat(1,PlotData{1}.dPrime{:});
pVals=cat(2,PlotData{1}.Pvals{:})';


H=Utilities.MultipleComparisonsCorrection(pVals,'method','fdr');
sum(H)
len = length(sum(H));

dp(~H)=nan;


pnl(1,2).select(); gca; hold on;
xlim([.25 NC+.25]);
xtickangle(45); set(gca, 'fontweight','bold');
plt.barpatch([1:3],abs(dp(:,[1:3])),'yl','DI','fontsize',6','t','Discriminability Index')
labelmat = Labels;
xticklabels(labelmat)
xtickangle(45);
 xlim([.25 NC+.25]); set(gca, 'fontweight','bold');
ax1 = gca;
ax1.FontWeight = 'bold';
set(ax1,'box','off');
ax1.XRuler.Axle.LineStyle = 'none';
ylim([1.2 2.1])

%% %% Figure 5 Correlation

%%Correlation

%Notes To plot Figure 1 Correlation

% step 1. put a debugging red dot at line 70 in FaceScratch.CorrAnalysis.m
% the line that begins plt.fig

%step 2. run the lines that follow :

NC =13;

plt.fig('units','inches','width',12,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1).select();

Labels2Incl = Labels;
LabelsAbbrev = Labels;


imagesc(rho);
axis image;

title('Correlation between conditions','fontsize',6);
set(gca,'XTick',1:NConds,'XTickLabel',LabelsAbbrev);
set(gca,'YTick',1:NConds,'YTickLabel',LabelsAbbrev);
set(gca,'YDir','normal')
colorbar;
colormap(jet);

axis tight;
xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',6);

%% Figure 5 Specificity

NC =13;

plt.fig('units','inches','width',9,'height',5,'font','Arial','fontsize',9);
pnl2 = panel();  pnl2.margin=15; pnl2.pack(2,3);
pnl2.fontsize=11; pnl2.fontname='Arial';

tmpInd=1;

for i = 1:NConds
    for j = (i+1):NConds
        leftCond = LabelsAbbrev{j};
        rightCond = LabelsAbbrev{i};
        pnl2(1,tmpInd).select();
        histogram(D{i,j},spc,'Normalization','pdf');
        
        [b1,~,stats] = signtest(D{i,j});
        b1
        stats
        m1 = median(D{i,j});
         plt.vline(m1,'r');
         if b1<0.05         
        title(sprintf('Median = %.2f \n p = %0.3g',m1,b1),'fontweight', 'bold');
         else
             title(sprintf('Median = %.2f \n p = %0.2f',m1,b1),'fontweight', 'bold');
         end
        xlabel(sprintf('%s < --- > %s',leftCond,rightCond),'fontweight', 'bold','fontsize',11);
        ylim([0 1.1]);
                axis tight; axis square;
        set(gca,'fontweight', 'bold', 'fontsize', 6);
        ax1.FontWeight = 'bold';
        tmpInd=tmpInd+1;
    end
end

%% Figure 3  Latency Main

% very important note. for this figure, make sure that around lines 29 or
% 30, use the version that says ResultsDir =
% fullfile(env.get('results'),'FaceScratch','SUAnal',[Tag
% '-Go'],'PopData'); and not % ResultsDir = fullfile(env.get('results'),'FaceScratch','SUAnal',[Tag '-Go'],'PopDataNewMod');

% To plot this figure go to the ImageryLatencyLR(varargin) section of
% Analyze.FaceScratch2.m. Run the first portion with Dates, etc. Then load
% alltrialdata. then, run the latency portion, , then debug marker at 499
% where it starts plt.fig. then run the following


plt.fig('units','inches','width',9,'height',3,'font','Arial','fontsize',9);
pnl = panel(); pnl.margin = 20; 
pnl.pack('h',{2/3 1/3});
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1).select();
ax1 = gca; hold on;

% Analyze.plotEventRelatedAverage(perfCVLRBoot,'nan','TimeVec',(tWs(:,1)+.025));

for i=1:2
    for j=1:8
%    A{i}(j,:)=smooth(perfCVLRDay{i}(j,:),10); 
   A{i}(j,:)=perfCVLRDay{i}(j,:); 
    end
end

Analyze.plotEventRelatedAverage(A,{'L','R'},'TimeVec',(tWs(:,1)+.025));

tmpLeg = legend('Left','Right','Location','southeast');
legend boxoff; set(tmpLeg,'FontWeight','bold');

% decodeLatMuLR = mean(decodeLatLR);

clrs = lines(2);
% tmp1 = plt.vline(decodeLatLRBoot(1),{'--','Color',clrs(1,:)});
% tmp2 = plt.vline(decodeLatLRBoot(2),{'--','Color',clrs(2,:)});

plt.vline(decodeLatLRBoot(1),{'--','Color',clrs(1,:)},sprintf(' %.2f ms',decodeLatLRBoot(1)*1000),[0 1.08],{'HorizontalAlignment','Left'});
plt.vline(decodeLatLRBoot(2),{'--','Color',clrs(2,:)},sprintf('%.2f ms ',decodeLatLRBoot(2)*1000),[0 1.08],{'HorizontalAlignment','Right'});


plt.vline(0,'k');
xlabel('Time relative to Touch Onset (s)', 'fontweight', 'bold', 'fontsize', 6);
ylabel('Decode performance', 'fontweight', 'bold', 'fontsize', 6);
title('Decode Performance Aligned by Touch Onset', 'fontweight', 'bold');
tmpTitleGCA = get(gca,'title');
tmpTitleGCA.Position(1) = tmpTitleGCA.Position(1)+0.25;
tmpTitleGCA.Position(2) = tmpTitleGCA.Position(2)+0.03;
xlim([-0.5 2.0])
set(gca,'fontweight', 'bold', 'fontsize', 6);
ax1.FontWeight = 'bold';
hold on; 
tmpP = patch([-0.1 0.3 0.3 -0.1],[0.2 0.2 0.98 0.98],'r','HandleVisibility','off');
set(tmpP,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.15);

%% Figure 3 Latency Inset

% very important note. for this figure, make sure that around lines 29 or
% 30, use the version that says  % ResultsDir = fullfile(env.get('results'),'FaceScratch','SUAnal',[Tag '-Go'],'PopDataNewMod');

% To plot this figure go to the ImageryLatencyLR(varargin) section of
% Analyze.FaceScratch2.m. Run the first portion with Dates, etc. Then load
% alltrialdata. then, run the latency portion, then debug marker at 499
% where it starts plt.fig. then run the following

plt.fig('units','inches','width',9,'height',3,'font','Arial','fontsize',9);
pnl = panel(); pnl.margin = 20; 
pnl.pack('h',{2/3 1/3});
pnl.fontsize=11; pnl.fontname='Arial';


pnl(2).select();
ax1 = gca; hold on;

for i=1:2
    for j=1:8
%    A{i}(j,:)=smooth(perfCVLRDay{i}(j,:),10); 
   A{i}(j,:)=perfCVLRDay{i}(j,:); 
    end
end

Analyze.plotEventRelatedAverage(A,{'L','R'},'TimeVec',(tWs(:,1)+.025)); hold on;

plot([0.05; 0.15], [0.35; 0.35], '-k', 'LineWidth', 2,'color',[0.4 .4 .4]);
text(.1,0.30, '100 ms', 'HorizontalAlignment','center')

% Analyze.plotEventRelatedAverage(perfCVLRBoot,'nan','TimeVec',(tWs(:,1)+.025));

% decodeLatMuLR = mean(decodeLatLR);

clrs = lines(2);
% plt.vline(decodeLatLRBoot(1),{'--','Color',clrs(1,:)},sprintf(' %.2f ms',decodeLatLRBoot(1)*1000),[0 1.1],{'HorizontalAlignment','Left'});
% plt.vline(decodeLatLRBoot(2),{'--','Color',clrs(2,:)},sprintf('%.2f ms ',decodeLatLRBoot(2)*1000),[0 1.1],{'HorizontalAlignment','Right'});
% plt.vline(decodeLatMeanFirst,'r--',sprintf('%.2f ms',decodeLatMeanFirst*1000));

set(gca,'XTick',[])

% plt.vline(0,'k');
xlabel('Time', 'fontweight', 'bold', 'fontsize', 6);
ylabel('Decode performance', 'fontweight', 'bold', 'fontsize', 6);
title('Shaded Region Expanded', 'fontweight', 'bold');
xlim([-0.25 0.20])
set(gca,'fontweight', 'bold', 'fontsize', 6);
ax1.FontWeight = 'bold';


%% Figure 3 Latencies

% very important note. for this figure, make sure that around lines 29 or
% 30, use the version that says ResultsDir =
% fullfile(env.get('results'),'FaceScratch','SUAnal',[Tag
% '-Go'],'PopData'); and not % ResultsDir = fullfile(env.get('results'),'FaceScratch','SUAnal',[Tag '-Go'],'PopDataNewMod');

% To plot this figure go to the ImageryLatencyLR(varargin) section of
% Analyze.FaceScratch2.m. Run the first portion with Dates, etc. Then load
% alltrialdata. then, run the latency portion, , then debug marker at 499
% where it starts plt.fig. then run the following

ttt = (-1.5:.005:1.5);

NShuffs = 100;
zzidx = find(ttt==0);
% FR = Analyze.getNeuralData(AlignedTrialData,ASF,ttt);
FR = Analyze.getNeuralData(AlignedTrialData,ASF,ttt,'smooth',10);
Conds = Analyze.returnFieldValues(AlignedTrialData,'Condition');
if all(LabelIdx2Incl)
    UnitLatencyFileName = fullfile(ResultsDir,'AlignedData','UnitLatencyResults3');
else
    UnitLatencyFileName = fullfile(ResultsDir,'AlignedData','UnitLatencyNoHand');
end

if anyRecalc || ~exist([UnitLatencyFileName '.mat'],'file')
    for i = 1:size(ASF,1)
        fprintf('Unit %d/%d...\n',i,size(ASF,1));
        for c = 1:length(UConds)
            cIdx = strcmp(Conds,UConds{c});
            sigg = nanmean(FR(cIdx,:,i));
            [iptAll{i,c},res(i,c)] = findchangepts(sigg,'statistic','linear','MaxNumChanges',1);
            [iptMean{i,c},resMean(i,c)] = findchangepts(sigg,'statistic','mean','MaxNumChanges',1);
            
            for j=1:NShuffs
                sigg2=Shuffle(sigg);
                [~,resShuff(i,j,c)]=findchangepts(sigg2,'statistic','linear','MaxNumChanges',1);
            end
            prct(i,c)=nnz(resShuff(i,:,c)<res(i,c))/NShuffs;
            SStot(i,c) = sum((sigg-mean(sigg)).^2);
            mdl = fitlm(ttt(1:end-1),sigg);
            linSS(i,c) = mdl.SSE;
            mdlSave{i,c} = mdl;
            
            FRtmp2=FR(cIdx,:,i);
            FRtmp2=FRtmp2(~isnan(FRtmp2(:,1)),:);
            
            [R2_2(i,c),SNR(i,c),cvR2(i,c)]=AU.VarianceExplainedByMean(FRtmp2');

        end
    end
    save(UnitLatencyFileName,'ttt','iptAll','res','resShuff','prct',...
        'SStot','linSS','mdlSave','R2_2','SNR','cvR2','iptMean','resMean');
%     save(UnitLatencyFileName,'ttt','iptAll','res','resShuff','prct','UnitLatR2');
%     save(UnitLatencyFileName,'ttt','LatencyUnitIdx','LatencyUnitIdxFirst');
else
    fprintf('Found unit latency file...loading\n')
    load(UnitLatencyFileName);
end
[IsSigPrct,alpha] = Utilities.MultipleComparisonsCorrection(prct,'method','fdr','perColumn');

R2VMean = 1-res./SStot;
R2VLin = 1-res./linSS;
R2Meas = {R2VMean, R2VLin, R2_2, cvR2};
R2MeasLabel = {'R2VMean','R2VLin','R2_2','cvR2'};
R2Cutoff = [.5 .2 .1 .05];

R2Meas = {R2_2};
R2MeasLabel = {'R2_2'};
R2Cutoff = [.2];


clear InclIdx;
clear latTime;
i=1;
    plt.fig('units','inches','width',8,'height',8,'font','Arial','fontsize',9);
    pnl = panel(); pnl.margin = 20; pnl.pack(3,2);
    
    for c = 1:length(UConds)
        pnl(floor((c-1)/2)+1,mod(c+1,2)+1).select(); ax1 = gca; hold on;
%         subplot(3,2,c);
        emptyCell = cellfun(@isempty,iptAll(:,c));
%         InclIdx{i}(:,c) = ~emptyCell & (R2Meas{i}(:,c) > R2Cutoff(i));
        InclIdx{i}(:,c) = ~emptyCell & (R2Meas{i}(:,c) > R2Cutoff(i)) & IsSigPrct(:,c);
        latTime{i}{c} = ttt([iptAll{InclIdx{i}(:,c),c}]);
        histogram(latTime{i}{c}+.025,-1.5:.05:1.5);
        hold on
        set(gca,'fontweight', 'bold', 'fontsize', 6);
ax1.FontWeight = 'bold';

        if ~isempty(latTime{i}{c})
            ksdensity(latTime{i}{c}+.025,'bandwidth',.025);
        end
   
        switch UConds{c}
            case 'CheekL'
                tmpCond{c} = 'Left Cheek';
            case 'CheekR'
                tmpCond{c} = 'Right Cheek';
            case 'HandL'
                tmpCond{c} = 'Left Hand';
            case 'HandR'
                tmpCond{c} = 'Right Hand';
            case 'ShoulderL'
                tmpCond{c} = 'Left Shoulder';
            case 'ShoulderR'
                tmpCond{c} = 'Right Shoulder';                
        end
        
        title(tmpCond{c},'fontweight', 'bold', 'fontsize',6);
        xlabel('Response Latency (s)','fontweight', 'bold', 'fontsize',6);
        ylabel('Number of Units (%)','fontweight', 'bold', 'fontsize',6);
        
              axis tight; xlim([-.25 1]); 
        tmp=latTime{i}{c}+.025;
        tmp(tmp>1 | tmp <-.25)=[];
        md=median(tmp);
        ylim([0 21]);
        plt.vline(md,{'k'});
        if strcmp(UConds{c},'CheekL') || strcmp(UConds{c},'CheekR') ||strcmp(UConds{c},'ShoulderR')
        text(md+0.04,19.5,sprintf('Median = %0.1f ms',md*1000),'fontsize',6', 'fontweight','bold')
        else
            text(md+0.07,19.5,sprintf('Median = %0.1f ms',md*1000),'fontsize',6', 'fontweight','bold')
        end
    end
    suptitle(['Latencies: ' R2MeasLabel{i}]);

   