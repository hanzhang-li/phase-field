% 寻找每个节点的关系矩阵 节点关联节点数以及每个节点关联节点
function [nbandp, mbandp]=nodelemlinkdp(snode,nnode,element)
nbandp=zeros(snode,1);
mbandp=zeros(snode,1);
for i1=1:length(element)
    for j1=1:nnode
        ip1=element(i1,j1);
        for j2=1:nnode
            ip2=element(i1,j2);
            if sum(ip2==mbandp(ip1,:))==0
                nbandp(ip1)=nbandp(ip1)+1;
                mbandp(ip1,nbandp(ip1))=ip2;
            end
        end
    end
end

% 对索引值进行排序
mbandp(mbandp==0)=10e9;
mbandp=sort(mbandp,2);
mbandp(mbandp==10e9)=0;