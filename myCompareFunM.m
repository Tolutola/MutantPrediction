function [r2,rmse,pcc]=myCompareFunM(TestData,expData)
ngenes=size(TestData,1);
r2=zeros(ngenes,1);
rmse=zeros(ngenes,1);
for i=1: ngenes
    y=TestData(i,:);
    yCalc=expData(i,:);
    [r2(i), rmse(i)]=rsquare(y',yCalc',true);
end

pcc=diag(corr(TestData',expData'))';
r2=r2';
rmse=rmse';