% 本程序用于寻找施加位移的单元
% 本程序用于添加固定约束，x1,x2为横坐标，y1，y2为纵坐标，a b为约束为0时表示约束，为1时为自由
function [choelem]=poselemall(forcexy,nodecoor,selement,element)
% err为网格误差
err1=0.0000001;
% position用于存储选中的节点坐标
choelem=[];
for st=1:size(forcexy,1)
    for i=1:selement
        for j=1:4
            if nodecoor(element(i,j),1)>=forcexy(st,1)-err1&&nodecoor(element(i,j),1)<=forcexy(st,2)+err1&& ...
                nodecoor(element(i,j),2)>=forcexy(st,3)-err1&&nodecoor(element(i,j),2)<=forcexy(st,4)+err1  % 选择相应的节点
                 choelem=[choelem; i element(i,j)];
            end
        end
    end
end

