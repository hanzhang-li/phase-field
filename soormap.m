% 本程序用于建立总结点与自由度对应矩阵
function [degssor]=soormap(sdegfree,snode,ndofn,issor)
degssor=zeros(sdegfree,1);
for i=1:snode
    for j=1:ndofn
        if issor(i,j)~=0
            degssor(issor(i,j))=2*(i-1)+j;          
        end
    end
end