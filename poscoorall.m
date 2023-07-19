% 本程序用于寻找施加约束的节点坐标
% 本程序用于添加固定约束，x1,x2为横坐标，y1，y2为纵坐标，a b为约束为0时表示约束，为1时为自由
function [position]=poscoorall(loadu,nodecoor,snode)
% err为网格误差
err1=0.000000001;
% position用于存储选中的节点坐标
max=snode;
position=zeros(max,3);
num=0;
for j=1:size(loadu,1)
    for i=1:snode
        % 选择相应的节点
        if loadu(j,5)==0&&nodecoor(i,1)>=loadu(j,1)-err1&&nodecoor(i,1)<=loadu(j,2)+err1 ...
                &&nodecoor(i,2)>=loadu(j,3)-err1&&nodecoor(i,2)<=loadu(j,4)+err1  
            num=num+1;
            position(num,:)=[i loadu(j,6) loadu(j,7)];
        end
    end
end
position(num+1:max,:)=[];