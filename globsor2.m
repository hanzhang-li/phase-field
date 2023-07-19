% 本程序用于建立整体自由度矩阵和节点的对应矩阵 globsor
% 时表示约束，为1时为自由
function [issor,sdegfree]=globsor2(degfree,snode)
issor=zeros(snode,2);
issor(1,1)=degfree(1,1);
issor(1,2)=degfree(1,2)+issor(1,1);
for i=2:snode
    issor(i,1)=degfree(i,1)+issor(i-1,2);
    issor(i,2)=degfree(i,2)+issor(i,1);
end
for i=1:snode
    if degfree(i,1)==0
        issor(i,1)=0;
    end
    if degfree(i,2)==0
        issor(i,2)=0;
    end
end
% 节点整体自由度sdegfree
sdegfree=max(max(issor));