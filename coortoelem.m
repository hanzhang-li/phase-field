% �����򽫽���ɽڵ�洢��Ϊ��Ԫ�洢��astrain���յ�Ԫ�洢
function [astrain,astress]=coortoelem(selement,avestrain,avestress,element)
% ��Ӧ��Ӧ�����յ�Ԫ�洢
astrain=zeros(selement,4,3); astress=zeros(selement,4,3);
for i=1:selement
    for j=1:4
        for k=1:3
           astrain(i,j,k)=avestrain(element(i,j),k);
           astress(i,j,k)=avestress(element(i,j),k);
        end
    end
end