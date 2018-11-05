% Thermodynamic framework for mutant phenotype prediction
% raw script to be functionalized

% Tola Oyetunde (10/08/2018)
clc
clear
load case1

% case studies
%case #1 4 e coli
%case #3 36 yeast


% inputs (loaded as part of the cases)
%************************
% gibbs free energy of reactions at biochemical conditions (dGr0)
% gibbs free energy of formation at biochemcial conditions (dGf0)
% measured metabolite distribution (old_mets)
%genome scale model (model)

% (all the below are saved as datasetE and datasetM)
% old wild type flux distribution
% measured mutant flux distribution


gen_preprocess % general data preprocessing and set up script

model=REMEP_preprocess(model); % generate variables needed for REMEP
model=MOMA_preprocess(model); % generate variables needed for MOMA
%%
numDels=length(model.datasetE.geneDeleted)-1;
Results=cell(numDels,1);


nMod=5; % number of test algorithms
dont_plot=[]; % indices of genes not found in the model

%REMEP options

modeFlag=3; %  3: classic REMEP
n_iter=1e3; % number of iterations for iterative flux-met computation
Mets=zeros(sum(model.goodMx),numDels);
wild_dGr=model.dGr0(model.goodRx)+ model.RT*model.small_S'*model.adj_old_mets;
dGs=zeros(sum(model.goodRx),numDels+1);
dGs(:,1)=wild_dGr;
for k=1:numDels
    deletedGene=model.datasetE.geneDeleted{k+1};
    tempResults =zeros(length(model.rxns),nMod);
    
    if ~isempty(find(ismember(model.genes,deletedGene))) || strcmp(model.datasetE.geneDeleted{1},'no_glucose')
        %FBA
        [modelPert, tempResults(:,1),err]=myFBA(model,deletedGene);


        [modelPert,Flux,err]=myMOMA(modelPert);
        tempResults(:,2)=Flux;
         

        %RELATCH
        % before adaptive laboratory evolution
        alpha = 10; gamma = 1.1;

        unevolved = RELATCH_Perturbed(modelPert,solutionRef,alpha,gamma,'cplex_direct');
        if unevolved.stat>0
            tempResults(:,3)=unevolved.v;
        else
            tempResults(:,3)=solutionRef.w;

        end
        modelPert.RELATCHsol_unevolved=tempResults(:,3);

        % after adaptive laboratory evolution
        alpha = 1; gamma = inf;
        evolved = RELATCH_Perturbed(modelPert,solutionRef,alpha,gamma,'cplex_direct');
        if evolved.stat>0
            tempResults(:,4)=evolved.v;
        else
            tempResults(:,4)=solutionRef.w;
        end
        
        modelPert.RELATCHsol_evolved=tempResults(:,4);


        %REMEP
        [tempResults(:,5),Mets(:,k),dGs(:,k+1)]=REMEP(modelPert,n_iter,modeFlag);
%         tempResults(:,5)=solutionDel.x; % dummy

        
    else
        dont_plot=[dont_plot;k];
    end
    Results{k,1}=tempResults;
    disp(k)
    
end

%% plots
sChoice=1;% REMEP
rChoice=3; % 3 unevolved 4 evolved RE LATCH

geneChoice=1:4;
geneChoice(dont_plot)=[];
geneLabels=model.datasetM.geneDeleted2(2:end);

compareInt % compare algorithms

