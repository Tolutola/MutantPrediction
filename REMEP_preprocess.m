function model=REMEP_preprocess(model)
%% compute REMEP matrices
% number of metabolites and reactions
[nm,nr]=size(model.S);

% number of irreverisble reactions
nr_r=length(model.modeIrrev.rxns);


% get rid of unwanted metabolites
%#1  external metabolites
% m,n, x,e g v f
metIDa=strfind(model.mets,'[e]');
metIDb=strfind(model.mets,'_e');
metIDc=strfind(model.mets,'[p]');
metIDd=strfind(model.mets,'_p');
metIDe=strfind(model.mets,'_m');
metIDf=strfind(model.mets,'_n');
metIDg=strfind(model.mets,'_x');
metIDh=strfind(model.mets,'_g');
metIDi=strfind(model.mets,'_v');
metIDj=strfind(model.mets,'_f');
metIDx=cell2mat(cellfun(@(a,b,c,d,e,f,g,h,i,j) isempty(a)& isempty(b)& isempty(c) & ...
    isempty(d)& isempty(e)&isempty(f)& isempty(g)& isempty(h)&isempty(i) & isempty(j), ...
    metIDa,metIDb,metIDc,metIDd,metIDe,metIDf,metIDg,metIDh,metIDi,metIDj,'UniformOutput',false));

goodNM=find(metIDx);
% #2 metabolites with no info on gibbs free energy of formation 
goodNM2= intersect(goodNM,model.thermo_met_match(:,1));


metThresh=0; % using both active and inactive metabolites for now
goodNM3=[];
% form C and d matrices, line up thermodynamic coefficients


C=[];
for i=1:nm
    if ismember(i,goodNM2)
        % generate PP, PC, CP and CC relationships
        tempVec=model.modeIrrev.S(i,:);
        tempVecP=tempVec; tempVecC=tempVec;
        tempVecP(tempVecP<0)=0;
        tempVecC(tempVecC>0)=0;
        if abs(tempVecP*model.RELATCHsol_wildtype_ir)>=metThresh % current metabolite is active

            C=[C;tempVecP];
            goodNM3=[goodNM3; i];
        end
    end
end

goodRx=false(nr,1);

% find reactions that have gibbs free energy of formation and for whom all
% the metabolites are in the 'good' list.
for i=1:size(model.S,2)
    if ismember(i,model.thermo_rxn_match(:,1))
         if isempty(setdiff(find(model.S(:,i)), goodNM3))
             goodRx(i)=true;
         end
    end
end

goodMx=false(nm,1);
goodMx(goodNM3)=true;


%% adjust given wild type metabolite concentrations
gnm=sum(goodMx);
gnr=sum(goodRx);

options = cplexoptimset;
options.Display = 'off';
options.TolFun=1e-9;
options.TolRLPFun=1e-9;
options.Algorithm='primal';

temp_mets=log(model.modeIrrev.old_mets(goodMx));

% compute gibbs free energy of formation that is consistent
% with metabolite concentrations

old_fluxes=model.FBAsol_wildtype(goodRx);

S=model.S;
hBool = strcmp(model.metFormulas,'H');

S(hBool,:) = 0; % Set proton coefficients to 0
S=S(goodMx,:);
S=S(:,goodRx);

used_mets=model.mets(goodMx);

lb_m=log(1e-8*ones(gnm,1));
ub_m=log(0.1*ones(gnm,1));

min_thermo_flux=5e-1;

% remove water and protons from consideration
h2o_ID=find(ismember(used_mets,'h2o_c'));
h_ID=find(ismember(used_mets,'h_c'));

lb_m([h2o_ID,h_ID])=0;
ub_m([h2o_ID,h_ID])=0;


neg_flag=old_fluxes<0 & abs(old_fluxes)>=min_thermo_flux;
nul_flag=abs(old_fluxes)<min_thermo_flux;

Aineq_m=model.RT*S';
dGr0=model.dGr0(goodRx);
bineq_m=-dGr0+ model.beta_g;

bineq_m(neg_flag)=-bineq_m(neg_flag);
Aineq_m(neg_flag,:)=-Aineq_m(neg_flag,:);

Aineq_m(nul_flag,:)=[];
bineq_m(nul_flag)=[];

d_m=temp_mets;
C_m=eye(gnm);

[adj_old_mets, resnorm, residual, exitflag, output, lambda] = cplexlsqlin (C_m, d_m, Aineq_m, bineq_m,...
    [], [], lb_m, ub_m, [], []);


%% load up model with new data
model.nr=nr;
model.nr_r=nr_r;
model.nm=nm;
model.gnm=gnm;
model.gnr=gnr;
model.goodMx=goodMx;
model.goodRx=goodRx;
model.C_remep=C;
model.d_remep=C*model.FBAsol_wildtype_ir;
model.min_thermo_flux=min_thermo_flux;
model.adj_old_mets=adj_old_mets(1:gnm);
model.small_S=S;
model.lb_m=lb_m;
model.ub_m=ub_m;
model.metThresh=metThresh;