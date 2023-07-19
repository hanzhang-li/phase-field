% 本段程序用于求解平均应力应变值
function [avestrain,avestress]=avestressfunop(element,strain,stress,linknode,snode,nnode)
avestrain=zeros(snode,3); avestress=zeros(snode,3);
[m,n]=size(strain);
for i=1:m/nnode
    for j=1:nnode
        for k=1:n
           avestrain(element(i,j),k)=avestrain(element(i,j),k)+strain(4*(i-1)+j,k)/linknode(element(i,j));
           avestress(element(i,j),k)=avestress(element(i,j),k)+stress(4*(i-1)+j,k)/linknode(element(i,j));
        end
    end
end  

