% general preprocessing script
%% solver set up
changeCobraSolver('ibm_cplex','LP');
changeCobraSolver('ibm_cplex','MILP');
changeCobraSolver('gurobi6','QP');
changeCobraSolver('gurobi6','MIQP');
%%
% get wild type data
External_Flux.rxns=model.datasetE.rxns;
External_Flux.val=model.datasetE.val(1,:);
External_Flux.err=model.datasetE.err(1,:);

MFA_Flux.rxns=model.datasetM.rxns;
MFA_Flux.val=model.datasetM.val(1,:);
MFA_Flux.err=model.datasetM.err(1,:);


%apply experimental fluxes to model
model_exp=applyExpFlux(model,External_Flux,MFA_Flux,1);

FBAsol=optimizeCbModel(model_exp,[],'one');
model.FBAsol_wildtype=FBAsol.x;
model.FBAsol_wildtype_ir=convertRevFluxDistribution(model.FBAsol_wildtype,model.rev2irrev,length(model.modeIrrev.rxns));

% RELATCH reference solution
solutionRef = RELATCH_Reference(model,model.gene_expression,External_Flux,MFA_Flux,'cplex_direct');

model.RELATCHsol_wildtype=solutionRef.w;
model.RELATCHsol_wildtype_ir=convertRevFluxDistribution(solutionRef.w,model.rev2irrev,length(model.modeIrrev.rxns));