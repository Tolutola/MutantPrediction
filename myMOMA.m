function [modelPert,Flux,err]=myMOMA(modelPert)

Aineq = [];
bineq = [];

Aeq=modelPert.modeIrrev.S;
beq=modelPert.modeIrrev.b;

%block deleted reactions
if ~isempty(modelPert.deleteID)
    for i=1:length(modelPert.deleteID)
        modelPert.modeIrrev.lb(modelPert.rev2irrev{modelPert.deleteID(i)})=0;
        modelPert.modeIrrev.ub(modelPert.rev2irrev{modelPert.deleteID(i)})=0;
    end
end

lb=modelPert.modeIrrev.lb;
ub=modelPert.modeIrrev.ub;

H=modelPert.Cmoma'*modelPert.Cmoma;
f=-modelPert.Cmoma'*modelPert.dmoma;

options = cplexoptimset;
options.Display = 'off';
options.TolFun=1e-9;
options.TolRLPFun=1e-9;
options.Algorithm='primal';

Aeq(:,modelPert.Excl_ID)=[];
lb(modelPert.Excl_ID)=[];
ub(modelPert.Excl_ID)=[];

[x_, fval, err, output] = cplexqp (H, f, Aineq, bineq, Aeq, beq, lb,...
    ub, [ ], options);
x=zeros(modelPert.nr_r,1);
x(modelPert.ntExcl_ID)=x_;
% convert irreversible flux distribution to reversible
Flux = convertIrrevFluxDistribution(x,modelPert.matchRev);
modelPert.MOMAsol=Flux;


