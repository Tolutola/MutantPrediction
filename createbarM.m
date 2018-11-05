function createbarM(sChoice,RMSE,geneChoice,geneLabels,rChoice)

for i =1:length(geneLabels)
    geneLabels{i}=strcat('\Delta',geneLabels{i});
end

figure1 = figure;
ymatrix1=RMSE([1 2 rChoice 4+sChoice],geneChoice)';

if size(ymatrix1,1)==1
    ymatrix1=ymatrix1';
end

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to bar
bar1 = bar(ymatrix1,'Parent',axes1);
set(bar1(1),'DisplayName','FBA','FaceColor',[1 0 0]);
set(bar1(2),'DisplayName','MOMA','FaceColor',[0 0 1]);
set(bar1(3),'DisplayName','RELATCH');
set(bar1(4),'DisplayName','REMEP','FaceColor',[0 1 0]);

% Create ylabel
ylabel('rmse');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',20,'TickLength',[0 0.025],'XTick',...
    1:length(geneChoice),'XTickLabel',geneLabels(geneChoice));
% Create legend
legend(axes1,'show');

