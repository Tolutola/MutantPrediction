%% comparing intracellular fluxes
TestDataM=model.datasetM.val(2:end,:);

nMeas=size(TestDataM,2);
tempData=zeros(nMeas,nMod,numDels);
nMeasVec=findRxnIDs(model,model.datasetM.rxns);

cnt=1;
for i=1: length(Results)
    if ~isempty(Results{i})
        tempData(:,:,cnt)=Results{i}(nMeasVec,:);
        cnt=cnt+1;
    end

end


R=zeros(nMod,numDels);
RMSE=R;
PCC=R;
SimDM=[];
relDeviations=cell(nMod,1);
TestDataM2=TestDataM';
TestDataM2(abs(TestDataM2)<=1e-6)=5e-4;
for i=1:nMod
    simDataM=squeeze(tempData(:,i,:));
    [R(i,:),RMSE(i,:),PCC(i,:)]=myCompareFunM(TestDataM,simDataM');
    SimDM(:,:,i)=simDataM;
    relDeviations{i}=(simDataM-TestDataM2)./TestDataM2;
end

mRMSE=mean(RMSE,2);

createplotsM(sChoice,SimDM,TestDataM,PCC,RMSE,geneChoice,geneLabels,rChoice)
createbarM(sChoice,RMSE,geneChoice,geneLabels,rChoice)
% createboxplot1(sChoice,relDeviations,geneChoice,geneLabels,rChoice)


