function [Flux,Mets,dG,err]=REMEP(model,n_iter,modeFlag)
Mets=zeros(model.gnm,1); dG=zeros(model.gnr,1); err=0;
Aeq=model.modeIrrev.S;
beq=model.modeIrrev.b;

%block deleted reactions
if ~isempty(model.deleteID)
    model.lb(model.deleteID)=0;
    model.ub(model.deleteID)=0;
    for i=1:length(model.deleteID)
        modelIrrev.lb(model.rev2irrev{model.deleteID(i)})=0;
        modelIrrev.ub(model.rev2irrev{model.deleteID(i)})=0;
    end
end

lb=model.modeIrrev.lb;
ub=model.modeIrrev.ub;

x0=[];

options = cplexoptimset;
options.Display = 'off';
options.TolFun=1e-9;
options.TolRLPFun=1e-9;
options.Algorithm='primal';

switch modeFlag
   
        
    case 3 %thermo REMEP
        
        % solve the classic remep case first
        H=model.C_remep'*model.C_remep;
        f=-model.C_remep'*model.d_remep;
        
        [c_x, fval, exitflag, output] = cplexqp (H, f, [], [], Aeq, beq,...
            lb, ub,x0);
        
        used_mets=model.mets(model.goodMx);
        curr_mets=model.adj_old_mets;
        lb_m=log(1e-8*ones(model.gnm,1));
        ub_m=log(0.1*ones(model.gnm,1));
        
        h2o_ID=find(ismember(used_mets,'h2o_c'));
        h_ID=find(ismember(used_mets,'h_c'));
        
        lb_m([h2o_ID,h_ID])=0;
        ub_m([h2o_ID,h_ID])=0;
        
        Aineq_m=model.RT*model.small_S';
        dGr0=model.dGr0(model.goodRx);
        bineq_m=-dGr0+ model.beta_g;
        
        curr_mets_coefs=abs(model.dGf0(model.goodMx)+ (curr_mets*model.RT))*1e-3;
        d1=curr_mets_coefs.*model.d_remep;
        %         model.C_remep
        options = cplexoptimset;
        options.Display = 'off';
        options.TolFun=1e-9;
        options.TolRLPFun=1e-9;
        options.Algorithm='primal';
  
        
        Nouse_Exch_ID=find(findExcRxns(model.modeIrrev) & (model.FBAsol_wildtype_ir==0));
        lb(Nouse_Exch_ID)=0;
        ub(Nouse_Exch_ID)=0;
        for k=1:n_iter
            %flux
            C1=model.C_remep.*repmat(curr_mets_coefs,1,model.nr_r);
          
            H=C1'*C1;
            f=-C1'*d1;
            
            
            [x, fval, exitflag, output] = cplexqp (H, f, [], [], Aeq, beq, lb,...
                ub, [],[]);
            
            if exitflag~=1
                [x, fval, exitflag, output] = cplexqp (H, f, [], [], Aeq, beq, lb,...
                    ub, [],options);
            end
            
            
            Flux = convertIrrevFluxDistribution(x,model.matchRev);
            temp_fluxes=Flux(model.goodRx);
            
            neg_flag=temp_fluxes<0 & abs(temp_fluxes)>=model.min_thermo_flux;
            nul_flag=abs(temp_fluxes)<model.min_thermo_flux;
            
            % metabolites
            
            bineq_m_temp=bineq_m;
            Aineq_m_temp=Aineq_m;
            
            bineq_m_temp(neg_flag)=-bineq_m_temp(neg_flag);
            Aineq_m_temp(neg_flag,:)=-Aineq_m_temp(neg_flag,:);
            
            Aineq_m_temp(nul_flag,:)=[];
            bineq_m_temp(nul_flag)=[];
            
            
            K=model.RT*eye(model.gnm);
            
            
            CV=model.C_remep*x;
            C_m=1e-3*diag(CV'*K);
            d_m=-(CV.*model.dGf0(model.goodMx)*1e-3)+ d1;
          
            
            H_m=C_m'*C_m;
            f_m=-C_m'*d_m;
            H_m=[H_m;eye(model.gnm)];
            
            [x_m, fval_m, exitflag_m, output_m] = cplexqp (H_m, f_m, Aineq_m_temp, bineq_m_temp, [], [], lb_m,ub_m);
            
            if exitflag_m~=1
                [x_m, fval_m, exitflag_m, output_m] = cplexqp (H_m, f_m, Aineq_m_temp, ...
                    bineq_m_temp, [], [], lb_m,ub_m,[],options);
            end
            
            if ~isempty(x_m)
                curr_mets=x_m;
            end
            curr_mets_coefs=abs(model.dGf0(model.goodMx)+ (curr_mets*model.RT))*1e-3;
            
    
                min_flux=x;
                min_mets=curr_mets;

        end
        
        
        Flux = convertIrrevFluxDistribution(min_flux,model.matchRev);
        err=output.cplexstatusstring;
        Mets=min_mets;
        dG=model.dGr0(model.goodRx)+ model.RT*model.small_S'*min_mets;
        
    case 4
        options=optimoptions('particleswarm','Display','iter','MaxTime',...
            100*60*60,'MaxStallIterations',100,'UseParallel',false);
        
        
        %variables (fluxes, metabolites)
        LB=[model.lb;model.lb_m];
        UB=[model.ub;model.ub_m];
        
        
        %figure out external reactions with no flux and block them
        Nouse_Exch_ID=find(findExcRxns(model) & (model.FBAsol_wildtype==0));
        LB(Nouse_Exch_ID)=0;
        UB(Nouse_Exch_ID)=0;
        
        model.Aineq_m=model.RT*model.small_S';
        dGr0=model.dGr0(model.goodRx);
        model.bineq_m=-dGr0+ model.beta_g;
        curr_mets=model.adj_old_mets;
        curr_mets_coefs=abs(model.dGf0(model.goodMx)+ (curr_mets*model.RT))*1e-3;
        model.old_met_sum=curr_mets_coefs.*model.d_remep;
        
     
        
        [x_sol,fval,err]= particleswarm(@(x) remep_obj(x,model),length(LB),LB,UB,options);
        Flux=x_sol(1:model.nr);
        Mets=x_sol(model.nr+1:end);
        
        
        
    case 5
        options=optimoptions('fmincon','Display','iter','UseParallel',false,...
            'MaxFunctionEvaluations',1e5);
        options=optimoptions('patternsearch','Display','iter','UseParallel',false,...
            'MaxFunctionEvaluations',1e5);
        
        
        %variables (fluxes, metabolites)
        LB=[model.lb;model.lb_m];
        UB=[model.ub;model.ub_m];
        
        
        %figure out external reactions with no flux and block them
        Nouse_Exch_ID=find(findExcRxns(model) & (model.FBAsol_wildtype==0));
        LB(Nouse_Exch_ID)=0;
        UB(Nouse_Exch_ID)=0;
        
        model.Aineq_m=model.RT*model.small_S';
        dGr0=model.dGr0(model.goodRx);
        model.bineq_m=-dGr0+ model.beta_g;
        curr_mets=model.adj_old_mets;
        curr_mets_coefs=abs(model.dGf0(model.goodMx)+ (curr_mets*model.RT))*1e-3;
        model.old_met_sum=curr_mets_coefs.*model.d_remep;
        
       
        x0=[model.FBAsol;curr_mets];
        A_eq=[model.S zeros(model.nm,model.gnm)];
        b_eq=model.b;
        [x_sol,fval,err]= patternsearch(@(x) remep_fmincon(x,model),x0,[],[],...
            A_eq,b_eq,LB,UB,@(x) remep_nonlcon(x,model),options);
        
        Flux=x_sol(1:model.nr);
        Mets=x_sol(model.nr+1:end);
        
        
end

function total_err=remep_obj(x, model)
fluxes=x(1:model.nr)';
mets=x(model.nr+1:end)';

% equality constraints for fluxes
flux_err=(model.S*fluxes).^2;

%inequality constraints for mets
bineq_m_temp=model.bineq_m;
Aineq_m_temp=model.Aineq_m;

temp_fluxes=fluxes(model.goodRx);
neg_flag=temp_fluxes<0 & abs(temp_fluxes)>=model.min_thermo_flux;
nul_flag=abs(temp_fluxes)<model.min_thermo_flux;

bineq_m_temp(neg_flag)=-bineq_m_temp(neg_flag);
Aineq_m_temp(neg_flag,:)=-Aineq_m_temp(neg_flag,:);

Aineq_m_temp(nul_flag,:)=[];
bineq_m_temp(nul_flag)=[];

met_err=Aineq_m_temp*mets - bineq_m_temp;
met_err=(met_err(met_err>0)).^2;

% objective function
fluxes_ir=convertRevFluxDistribution(fluxes,model.rev2irrev,length(model.modeIrrev.rxns));
curr_mets_coefs=abs(model.dGf0(model.goodMx)+ (mets*model.RT))*1e-3;
obj_err=(model.old_met_sum - (curr_mets_coefs.*(model.C_remep*fluxes_ir))).^2;

total_err=100*sum(met_err) + 1e3*sum(flux_err) + sum(obj_err);

function obj_err=remep_fmincon(x, model)
fluxes=x(1:model.nr);
mets=x(model.nr+1:end);
% objective function
fluxes_ir=convertRevFluxDistribution(fluxes,model.rev2irrev,length(model.modeIrrev.rxns));
curr_mets_coefs=abs(model.dGf0(model.goodMx)+ (mets*model.RT))*1e-3;
obj_err=sum((model.old_met_sum - (curr_mets_coefs.*(model.C_remep*fluxes_ir))).^2);

function [c,c_eq]=remep_nonlcon(x,model)
c_eq=[];

fluxes=x(1:model.nr);
mets=x(model.nr+1:end);

%inequality constraints for mets
bineq_m_temp=model.bineq_m;
Aineq_m_temp=model.Aineq_m;

temp_fluxes=fluxes(model.goodRx);
neg_flag=temp_fluxes<0 & abs(temp_fluxes)>=model.min_thermo_flux;
nul_flag=abs(temp_fluxes)<model.min_thermo_flux;

bineq_m_temp(neg_flag)=-bineq_m_temp(neg_flag);
Aineq_m_temp(neg_flag,:)=-Aineq_m_temp(neg_flag,:);

Aineq_m_temp(nul_flag,:)=[];
bineq_m_temp(nul_flag)=[];

if isempty(Aineq_m_temp)
    c=0;
else
    c=Aineq_m_temp*mets - bineq_m_temp;
end

