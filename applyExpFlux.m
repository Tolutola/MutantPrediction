function model=applyExpFlux(model,External_Flux,MFA_Flux,allFlag)

solutionFBA = optimizeCbModel(model,[],'one'); % Solve FBA problem
if isempty(solutionFBA.x)
    solutionFBA = optimizeCbModel(model); % Solve FBA problem
end

if ~isempty(External_Flux)
    [~, Ext_ID]=ismember(External_Flux.rxns,model.rxns);
    mask= ~ismember(External_Flux.rxns,{'EX_o2(e)','EX_co2(e)'});
    if find(Ext_ID==0)
        error('Exchange reactions not found in the model');
    else
        exID=External_Flux.err==0 & External_Flux.val~=0;
        External_Flux.err(exID)=max(abs(0.001*External_Flux.val(exID)),1e-5);
%         Nouse_Exch_ID=find(findExcRxns(model) & (solutionFBA.x==0));
%         model.lb(Nouse_Exch_ID)=0;
%         model.ub(Nouse_Exch_ID)=0;
        model.lb(Ext_ID(mask)) = External_Flux.val(mask) - External_Flux.err(mask);
        model.ub(Ext_ID(mask)) = External_Flux.val(mask) + External_Flux.err(mask);
    end
    disID2=(~model.rev(Ext_ID))& model.lb(Ext_ID)<0;
    model.lb(Ext_ID(disID2))=0;
    
end

if model.chemostat~=0
    model=changeRxnBounds(model,'Biomass',model.chemostat,'b');
end

solutionFBA = optimizeCbModel(model,[],'one'); % Solve FBA problem
if allFlag==1
    if ~isempty(MFA_Flux)

        for k=1:length(MFA_Flux.rxns)
            MFA_ID(k)=findRxnIDs(model,MFA_Flux.rxns{k});

                
        end
        if find(MFA_ID==0)
            error('MFA reactions not found in the model');
        else
            mfID=MFA_Flux.err==0 & MFA_Flux.val~=0;
          
            MFA_Flux.err(mfID)=max(abs(.1*MFA_Flux.val(mfID)),1e-3);% x% error if no error is reported
           
            model.lb(MFA_ID(mfID)) = (MFA_Flux.val(mfID) - MFA_Flux.err(mfID));
            model.ub(MFA_ID(mfID)) = (MFA_Flux.val(mfID) + MFA_Flux.err(mfID));
        end
        disID=(~model.rev(MFA_ID))& model.lb(MFA_ID)<0;
        model.lb(MFA_ID(disID))=0;
    end
    
solutionFBA = optimizeCbModel(model,[],'one'); % Solve FBA problem 
nr=length(model.modeIrrev.rxns);
model.flux_irrev = convertRevFluxDistribution(solutionFBA.x,model.rev2irrev,nr);
end
