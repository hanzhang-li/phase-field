% 本程序用于求解单元劲度矩阵        
function [delemk]=elemstmasm2(nnode,tnijmart,tndjmart,dnodecoor,ephi,xk,dmat)
    % 单元劲度矩阵初始化
    % delemk矩阵赋初值，单元刚度矩阵
    delemk=zeros(8,8);
    % 计算单元劲度矩阵
    for j=1:nnode
        %建立B矩阵即应变与位移关系矩阵
        [detjac,b,~]=pstrainlinkdisp(tndjmart(2*(j-1)+(1:2),:),dnodecoor);
        % 计算节点phi值
        phi=min(1,tnijmart(j,:)*ephi);
        % 计算delemk    添加加速度项
        delemk=delemk+detjac*((1-phi)^2+xk)*b'*dmat*b;
    end
end