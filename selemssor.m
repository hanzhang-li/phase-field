% 本程序用于建立每个单元的节点自由度矩阵
function [elmssor]=selemssor(selement,nnode,ndofn,issor,element)
elmssor=zeros(selement,nnode*ndofn);
for i=1:selement
    elmssor(i,ndofn*(1:nnode)-1)=issor(element(i,1:nnode),1);    
    elmssor(i,ndofn*(1:nnode))=issor(element(i,1:nnode),2); 
end