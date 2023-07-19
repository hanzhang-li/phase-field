% 本程序用于寻找施加位移的单元
% 本程序用于添加固定约束，x1,x2为横坐标，y1，y2为纵坐标
function [dispelem,fpelem,fpnode]=poslinkelem2(position2,element,dispelem,nnode)
    for i=1:nnode
        [~,row]=ismember(position2,element(:,i));
        % 找到包含红点的单元，写入红单元矩阵
        dispelem=[dispelem; row];
    end
    % 删去重复单元
    dispelem=unique(dispelem);
    dispelem(dispelem==0)=[];
    fpelem=element(dispelem,:);
    fpnode=unique(fpelem);
end
