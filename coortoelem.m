% 本程序将结果由节点存储改为单元存储，astrain按照单元存储
function [astrain,astress]=coortoelem(selement,avestrain,avestress,element)
% 将应变应力按照单元存储
astrain=zeros(selement,4,3); astress=zeros(selement,4,3);
for i=1:selement
    for j=1:4
        for k=1:3
           astrain(i,j,k)=avestrain(element(i,j),k);
           astress(i,j,k)=avestress(element(i,j),k);
        end
    end
end