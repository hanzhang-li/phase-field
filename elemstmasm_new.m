% 本程序用于求解单元劲度矩阵        
function [delemk]=elemstmasm_new(nnode,tndjmart,dnodecoor,dcmat)
    % 单元劲度矩阵初始化
    % delemk矩阵赋初值，单元刚度矩阵
    delemk=zeros(8,8);
    % 计算单元劲度矩阵
    for j=1:nnode
        cmat=dcmat(3*(j-1)+(1:3),:);
        %建立B矩阵即应变与位移关系矩阵
        [detjac,b,~]=pstrainlinkdisp(tndjmart(2*(j-1)+(1:2),:),dnodecoor);
        % 计算delemk    添加加速度项
        delemk=delemk+detjac*b'*cmat*b;
    end
end