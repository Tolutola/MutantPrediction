function solutionPert = RELATCH_Perturbed(model,solutionRef,alpha,gamma,solver)
%RELATCH_Perturbed predicts the flux distribution in a perturbed state
%
% solutionPert = RELATCH_perturbed(model,solutionRef,alpha,gamma,solver)
%
%INPUTS
% model					Perturbed metabolic model (CobraToolBox model with GPR)
% solutionRef
%	solutionRef.w		Reference flux distribution
%	solutionRef.W		Reference enzyme contributions
% alpha					penalty for latent pathway activation
% gamma					limit on active enzyme contribution increases
% solver					QP solver ('CobraQP' (default) or 'cplex_direct')
%
%OUTPUTS
% solutionPert
%	solutionPert.v		Perturbed flux distribution
%	solutionPert.stat	Perturbed solution status
%
%Notes:
% RELATCH_Perturbed requires CobraToolBox 2.0 and a QP solver
%
% Joonhoon Kim and Jennifer L. Reed 9/13/2012

if (nargin < 2)
    error('Reference flux distribution and enzyme contributions not found');
end
if (nargin < 3)
    alpha=1;
end
if (nargin < 4)
    gamma=1.1;
end
if (nargin < 5)
    solver = 'CobraQP';
end

% Generate GPR mapping
fprintf('Generating GPR mapping..\n');

[nMets,nRxns] = size(model.S);
% Find reactions with GPR
GPR_ID=find(~cellfun(@isempty,model.rules));
% Exclude Porin transport reactions and spontaneous reactions
Porin_ID=find(~cellfun(@isempty,strfind(model.subSystems,'Porin')));
Spont_ID=find(strcmp(model.grRules,'s0001'));
Excl_Rxn_ID=union(Porin_ID,Spont_ID);
GPR_ID=setdiff(GPR_ID,Excl_Rxn_ID);
% Find reversible reactions with GPR
[GPRrev_ID, GPRrev_ID2]=intersect(GPR_ID,find(model.rev));
nRxnsGPR=size(GPR_ID,1);
nRxnsGPRrev=size(GPRrev_ID,1);

% Generate GPR matrix
GPR=regexp(model.rules(GPR_ID),'\|','split');
GPR=cellfun(@(c) regexp(c(:),'x\((\d+)\)','tokens'), GPR, 'UniformOutput', 0);
Enz=cellfun(@(c) str2double([c{:}]), vertcat(GPR{:}), 'UniformOutput',0);
nEnzs=size(Enz,1);
nGenes=size(model.genes,1);

numEnz=zeros(nRxnsGPR,1);
for i=1:nRxnsGPR; numEnz(i)=size(GPR{i},1); end;
Rxn2Enz=sparse(nRxnsGPR,nEnzs);
for i=1:nRxnsGPR; Rxn2Enz(i,sum(numEnz(1:i-1))+1:sum(numEnz(1:i)))=ones(1,numEnz(i)); end;
numSub=zeros(nEnzs,1);
for i=1:nEnzs; numSub(i)=numel(Enz{i}); end;
Enz2Gene=sparse(nEnzs,nGenes);
for i=1:nEnzs; Enz2Gene(i,Enz{i})=1; end;

% Setting up the problem
% Variables in the following problem are
% x = [v;V]
% where v = flux distribution in the perturbed state
%       V = Enzyme contribution in the perturbed state
fprintf('Setting up the problem..\n');

model.lb(model.lb==-1000)=-inf;
model.ub(model.ub==1000)=inf;

F1 = sparse(nRxns,nRxns);
c1 = zeros(nRxns,1);
c2 = zeros(nEnzs,1);
Enz_Limit = 1000*ones(nEnzs,1);

Act_Rxn_ID = find(solutionRef.w ~= 0);
Exch_Rxn_ID = findExcRxns(model);
Act_Rxn_ID = setdiff(Act_Rxn_ID,union(Exch_Rxn_ID,Excl_Rxn_ID));
Act_Enz_ID = find(solutionRef.W ~= 0);
Inact_Enz_ID = find(solutionRef.W == 0);

Gene_status = true(size(model.genes));
Gene_status(~cellfun('isempty',(regexp(model.genes,'_deleted')))) = false;
Enz_status = true(nEnzs,1);
Enz_status = ~(Enz2Gene*~Gene_status);

F1(Act_Rxn_ID,Act_Rxn_ID) = diag(solutionRef.w(Act_Rxn_ID).^(-2));
c1(Act_Rxn_ID) = -solutionRef.w(Act_Rxn_ID).^(-1);
c2(Inact_Enz_ID) = 0.5*alpha;
if gamma ~= inf
    Enz_Limit(Act_Enz_ID) = gamma*solutionRef.W(Act_Enz_ID);
end
Enz_Limit(~Enz_status) = 0;

QP.F = [ F1 sparse(nRxns,nEnzs) ;
         sparse(nEnzs,nRxns+nEnzs) ];
QP.c = [ c1 ; c2 ];
QP.A = [ model.S sparse(nMets,nEnzs) ;
         sparse(1:nRxnsGPR,GPR_ID,-1,nRxnsGPR,nRxns) Rxn2Enz ;
         sparse(1:nRxnsGPRrev,GPRrev_ID,1,nRxnsGPRrev,nRxns) Rxn2Enz(GPRrev_ID2,:) ];
QP.b = zeros(nMets+nRxnsGPR+nRxnsGPRrev,1);
QP.lb = [ model.lb ; zeros(nEnzs,1) ];
QP.ub = [ model.ub ; Enz_Limit ];
QP.csense(1:nMets) = 'E';
QP.csense(nMets+1:nMets+nRxnsGPR+nRxnsGPRrev) = 'G';
QP.osense = 1;

% Solve the problem
fprintf('Solving perturbed state prediction problem..\n');

if strcmp(solver,'CobraQP')
	% Solve the reference state estimation problem using CobraQP solver
	QPsol = solveCobraQP(QP, 'printLevel',1);
	if (QPsol.stat > 0)
		solutionPert.v = QPsol.full(1:nRxns);
		solutionPert.stat = QPsol.stat;
		fprintf('Done\n\n');
	else
		solutionPert.stat = QPsol.stat;
		warning('Perturbed state prediction problem is infeasible or unbounded');
	end
elseif strcmp(solver,'cplex_direct')
	% Solve the reference state estimation problem using ILOG CPLEX
	[x,obj,exitflag,output] = cplexqp(QP.F,QP.c,-QP.A(nMets+1:end,:),-QP.b(nMets+1:end),QP.A(1:nMets,:),QP.b(1:nMets),QP.lb,QP.ub,[solutionRef.w;solutionRef.W]);
	if (exitflag > 0)
		solutionPert.v = x(1:nRxns);
		solutionPert.stat = exitflag;
		solutionPert.output = output;
		fprintf('Done\n\n');
	else
		solutionPert.stat = exitflag;
		solutionPert.output = output;
		warning('Perturbed state prediction problem is infeasible or unbounded');
	end
else error('Unknown QP solver');
end

end