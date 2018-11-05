function createplotsM(sChoice,SimD,TestData,PCC,RMSE,geneChoice,geneLabels,rChoice)

rmse=round(RMSE([1 2 rChoice 4+sChoice],geneChoice),1);
pcc=round(PCC([1 2 rChoice 4+sChoice],geneChoice),2);

ngenes=length(geneChoice);
figure1 = figure;
m=1; %FBA
for i=1:ngenes
    k=geneChoice(i);
    currentPlot = subplot(4,ngenes,i,'Parent',figure1);
    hold(currentPlot,'on');
    scatter(TestData(k,:),SimD(:,k,m)','Parent',currentPlot,'MarkerFaceColor',[1 0 0],...
        'MarkerEdgeColor',[1 0 0],...
        'LineWidth',0.8);
    title(geneLabels{k})
    %     minval=min(min(TestData(k,:),SimD(:,k,m)'));
    %     maxval=max(max(TestData(k,:),SimD(:,k,m)'));
    minval=-20;
    maxval=20;
    
    %     xlim([minval,maxval])
    %     ylim([minval,maxval])
    xlim([-20,20])
    ylim([-20,20])
    xlabel('Expt');
    ylabel('Model');
    box(currentPlot,'on');
    
    %textbox for PCC and rmse
    x1 = minval+ .2*(maxval-minval);
    y1 = minval+ .8*(maxval-minval);
    %PCC
    txt1 = strcat('r = ', num2str(pcc(m,i)));
    text(x1,y1,txt1,'FontWeight','bold','FontSize',20,'FontName','Times New Roman')
    %rmse
    y1 = minval+ .7*(maxval-minval);
    txt1 = strcat('rmse = ', num2str(rmse(m,i)));
    text(x1,y1,txt1,'FontWeight','bold','FontSize',20,'FontName','Times New Roman')
    
    % diagonal line
    plot(linspace(minval,maxval,10),linspace(minval,maxval,10),'Parent',currentPlot,'MarkerFaceColor',[1 0 0],...
        'MarkerEdgeColor',[1 0 0],...
        'LineWidth',2);
    
    % Set the remaining axes properties
    set(currentPlot,'FontName','Times New Roman','FontSize',20);
end

m=2; %MOMA
for i=1:ngenes
    k=geneChoice(i);
    currentPlot = subplot(4,ngenes,i+ngenes,'Parent',figure1);
    hold(currentPlot,'on');
    
    
    scatter(TestData(k,:),SimD(:,k,m)','Parent',currentPlot,'MarkerFaceColor',[1 0 1],...
        'MarkerEdgeColor',[1 0 1],...
        'LineWidth',0.8);
    %     minval=min(min(TestData(k,:),SimD(:,k,m)'));
    %     maxval=max(max(TestData(k,:),SimD(:,k,m)'));
    
    minval=-20;
    maxval=20;
    %     [minval, maxval]=deal(xlim);
    
    %     xlim([minval,maxval])
    %     ylim([minval,maxval])
    xlim([-20,20])
    ylim([-20,20])
    
    xlabel('Expt');
    ylabel('Model');
    box(currentPlot,'on');
    
    %textbox for PCC and rmse
    x1 = minval+ .2*(maxval-minval);
    y1 = minval+ .8*(maxval-minval);
    %PCC
    txt1 = strcat('r = ', num2str(pcc(m,i)));
    text(x1,y1,txt1,'FontWeight','bold','FontSize',20,'FontName','Times New Roman')
    %rmse
    y1 = minval+ .7*(maxval-minval);
    txt1 = strcat('rmse = ', num2str(rmse(m,i)));
    text(x1,y1,txt1,'FontWeight','bold','FontSize',20,'FontName','Times New Roman')
    
    % diagonal line
    plot(linspace(minval,maxval,10),linspace(minval,maxval,10),'Parent',currentPlot,'MarkerFaceColor',[1 0 1],...
        'MarkerEdgeColor',[1 0 1],...
        'LineWidth',2);
    
    % Set the remaining axes properties
    set(currentPlot,'FontName','Times New Roman','FontSize',20);
end

m=3; %RELATCH
for i=1:ngenes
    k=geneChoice(i);
    currentPlot = subplot(4,ngenes,i+2*ngenes,'Parent',figure1);
    hold(currentPlot,'on');
    scatter(TestData(k,:),SimD(:,k,m)','Parent',currentPlot,'MarkerFaceColor',[0 0 1],...
        'MarkerEdgeColor',[0 0 1],...
        'LineWidth',0.8);
    %     minval=min(min(TestData(k,:),SimD(:,k,m)'));
    %     maxval=max(max(TestData(k,:),SimD(:,k,m)'));
    minval=-20;
    maxval=20;
    %     [minval, maxval]=deal(xlim);
    
    %     xlim([minval,maxval])
    %     ylim([minval,maxval])
    xlim([-20,20])
    ylim([-20,20])
    
    xlabel('Expt');
    ylabel('Model');
    box(currentPlot,'on');
    
    %textbox for PCC and rmse
    x1 = minval+ .2*(maxval-minval);
    y1 = minval+ .8*(maxval-minval);
    %PCC
    txt1 = strcat('r = ', num2str(pcc(m,i)));
    text(x1,y1,txt1,'FontWeight','bold','FontSize',20,'FontName','Times New Roman')
    %rmse
    y1 = minval+ .7*(maxval-minval);
    txt1 = strcat('rmse = ', num2str(rmse(m,i)));
    text(x1,y1,txt1,'FontWeight','bold','FontSize',20,'FontName','Times New Roman')
    
    % diagonal line
    plot(linspace(minval,maxval,10),linspace(minval,maxval,10),'Parent',currentPlot,'MarkerFaceColor',[0 0 1],...
        'MarkerEdgeColor',[0 0 1],...
        'LineWidth',2);
    
    % Set the remaining axes properties
    set(currentPlot,'FontName','Times New Roman','FontSize',20);
end


m=4+sChoice; %sMEP
for i=1:ngenes
    k=geneChoice(i);
    currentPlot = subplot(4,ngenes,i+3*ngenes,'Parent',figure1);
    hold(currentPlot,'on');
    scatter(TestData(k,:),SimD(:,k,m)','Parent',currentPlot,'MarkerFaceColor',[0 1 0],...
        'MarkerEdgeColor',[0 1 0],...
        'LineWidth',0.8);
    %     minval=min(min(TestData(k,:),SimD(:,k,m)'));
    %     maxval=max(max(TestData(k,:),SimD(:,k,m)'));
    minval=-20;
    maxval=20;
    %     [minval, maxval]=deal(xlim);
    
    %     xlim([minval,maxval])
    %     ylim([minval,maxval])
    xlim([-20,20])
    ylim([-20,20])
    
    xlabel('Expt');
    ylabel('Model');
    box(currentPlot,'on');
    
    %textbox for PCC and rmse
    x1 = minval+ .2*(maxval-minval);
    y1 = minval+ .8*(maxval-minval);
    %PCC
    txt1 = strcat('r = ', num2str(pcc(m-sChoice,i)));
    text(x1,y1,txt1,'FontWeight','bold','FontSize',20,'FontName','Times New Roman')
    %rmse
    y1 = minval+ .7*(maxval-minval);
    txt1 = strcat('rmse = ', num2str(rmse(m-sChoice,i)));
    text(x1,y1,txt1,'FontWeight','bold','FontSize',20,'FontName','Times New Roman')
    
    % diagonal line
    plot(linspace(minval,maxval,10),linspace(minval,maxval,10),'Parent',currentPlot,'MarkerFaceColor',[0 1 0],...
        'MarkerEdgeColor',[0 1 0],...
        'LineWidth',2);
    
    % Set the remaining axes properties
    set(currentPlot,'FontName','Times New Roman','FontSize',20);
end

