function model=MOMA_preprocess(model)

old_fluxes=model.RELATCHsol_wildtype;
old_fluxes_ir=model.RELATCHsol_wildtype_ir;
Porin_ID=find(~cellfun(@isempty,strfind(model.subSystems,'Porin')));
Spont_ID=find(strcmp(model.grRules,'s0001'));
Excl_Rxn_ID=union(Porin_ID,Spont_ID);
Act_Rxn_ID = find(old_fluxes ~= 0);
Exch_Rxn_ID = findExcRxns(model);
Act_Rxn_ID = setdiff(Act_Rxn_ID,union(Exch_Rxn_ID,Excl_Rxn_ID));

Act_ID=[];
for i=1:length(Act_Rxn_ID)
    if length(model.rev2irrev{Act_Rxn_ID(i)})~=1
        Act_ID(end+1:end+2,1)=model.rev2irrev{Act_Rxn_ID(i)};
    else
        Act_ID(end+1,1)=model.rev2irrev{Act_Rxn_ID(i)};
    end
end

Nouse_Exch_ID=find(findExcRxns(model) & (old_fluxes==0));

Excl_ID=[];
for i=1:length(Nouse_Exch_ID)
    if length(model.rev2irrev{Nouse_Exch_ID(i)})~=1
        Excl_ID(end+1:end+2,1)=model.rev2irrev{Nouse_Exch_ID(i)};
    else
        Excl_ID(end+1,1)=model.rev2irrev{Nouse_Exch_ID(i)};
    end
end


fluxThresh=1e-5;


% MOMA
Cmoma=zeros(model.nr_r);%regular MOMA
act_C=abs(old_fluxes_ir(Act_ID))>=fluxThresh;
Cmoma(Act_ID(act_C),Act_ID(act_C))=diag(ones(length(Act_ID(act_C)),1));
dmoma=zeros(model.nr_r,1);
dmoma(Act_ID(act_C))=old_fluxes_ir(Act_ID(act_C));


Cmoma(Excl_ID,:)=[];
Cmoma(:,Excl_ID)=[];
dmoma(Excl_ID)=[];


model.Cmoma=Cmoma;
model.dmoma=dmoma;
model.Excl_ID=Excl_ID;
model.ntExcl_ID=setdiff(1:model.nr_r,Excl_ID)';
model.fluxThresh=fluxThresh;

