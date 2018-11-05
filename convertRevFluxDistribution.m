function vOldIrrev = convertRevFluxDistribution(vOld,rev2irrev,nr)

vOldIrrev=zeros(nr,1);
for i=1:length(vOld)
    if length(rev2irrev{i})==1
        vOldIrrev(rev2irrev{i})=vOld(i);
         
    else
        if vOld(i)>0
            vOldIrrev(rev2irrev{i})=[vOld(i) 0];
        else
            vOldIrrev(rev2irrev{i})=[0 -vOld(i)];
        end
    end
end
