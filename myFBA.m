function [modelPert,Flux,err]=myFBA(model,deletedGene)


modelPert=model;
if strcmp(model.datasetE.geneDeleted{1},'no_glucose')
    modelPert.lb(model.csourceID)=0;
    modelPert.ub(model.csourceID)=0;
    modelPert.modeIrrev.lb(model.rev2irrev{model.csourceID})=0;
    modelPert.modeIrrev.ub(model.rev2irrev{model.csourceID})=0;

    new_csourceID=findRxnIDs(model,deletedGene); %new carbon source
    modelPert.lb(new_csourceID)=-1000;
    modelPert.ub(new_csourceID)=0;

    modelPert.modeIrrev.lb(model.rev2irrev{new_csourceID})=0;
    modelPert.modeIrrev.ub(model.rev2irrev{new_csourceID})=1000;
    delR=[];
else
    [modelPert,~,delR] = deleteModelGenes(modelPert,deletedGene);
end

fSolution=optimizeCbModel(modelPert,[],'one');
if fSolution.f==0
    fSolution=optimizeCbModel(modelPert);
end
Flux=fSolution.x;
if ~isempty(delR)
    [~, deleteID]=ismember(delR,model.rxns);
else
    deleteID=[];
end

if fSolution.f==0
    err=1e5;
else
    err=0;
end

modelPert.delR=delR;
modelPert.deleteID=deleteID;
modelPert.FBAsol=Flux;


