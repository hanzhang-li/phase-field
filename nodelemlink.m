% 寻找每个节点的关系矩阵
function [nband, mband]=nodelemlink(snode,nnode,element,selement)
nband=zeros(snode,1);
mband=zeros(snode,1);
for i1=1:selement
    for j1=1:nnode
        ip1=element(i1,j1);
        for j2=1:nnode
            ip2=element(i1,j2);
            if ip1>=ip2
                if sum(ip2==mband(ip1,:))==0
                     nband(ip1)=nband(ip1)+1;
                     mband(ip1,nband(ip1))=ip2;
                end
            end
        end
    end
end

% 对索引值进行排序
mband(mband==0)=10e9;
mband=sort(mband,2);
mband(mband==10e9)=0;