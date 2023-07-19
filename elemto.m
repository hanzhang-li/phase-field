% 本程序用于调整计算单元局部节点坐标以计算
function [newelement]=elemto(element,nodecoor)
% 存储新的单元
newelement=element;
% 坐标轴主向量
svert=[1 0];
for i=1:size(element,1)
    % 计算单元形心
    delement=element(i,:);
    elemcen=0.25*sum(nodecoor(delement,:),1);
    
    % 装载单元节点
    elemface=nodecoor(delement,:);
    
    % 计算向量
    vert=elemface-elemcen;
    
    % 计算向量夹角
    theta=zeros(4,1);  % 余弦
    angel=zeros(4,1);  % 角度
    for j=1:4
        theta(j)=acos(dot(vert(j,:),svert)/norm(vert(j,:)));
        if vert(j,2)>=0
            angel(j)=theta(j);
        else
            angel(j)=2*pi-theta(j);
        end
    end
    
    % 排序
    [~,index]=sort(angel);
    
    % 装入新单元中
    newelement(i,:)=delement(index);
end