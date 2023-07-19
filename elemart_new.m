% 本程序用于装载单元矩阵
function [delemk,rhsk,delemphi,rhsp]=elemart_new ...
    (tdetjac,tbmat,tndx,ephi,ehm,estres,dcmat,tnijmart,gc,xlc,nnode)
    % delemk矩阵赋初值，单元刚度矩阵
    delemk=zeros(8,8);
    % delemphi矩阵赋初值，存储单元相场
    delemphi=zeros(4,4);
    % 右端项矩阵初始化
    rhsk=zeros(8,1);
    rhsp=zeros(4,1);
    % 计算单元劲度矩阵
    for j=1:nnode
        % 导入计算所需的单元矩阵
        nij=tnijmart(j,:);
        detjac=tdetjac(j);
        bmat=tbmat(3*(j-1)+1:3*(j-1)+3,:);
        ndx=tndx(2*(j-1)+1:2*(j-1)+2,:);
        cmat=dcmat(3*(j-1)+(1:3),:);

        % 计算节点phi值
        phi=min(1,nij*ephi);

        % 防止裂缝闭合，选取历史变量中最大值
        psi=ehm(j);

        % 计算delemk    添加加速度项
        delemk=delemk+detjac*bmat'*cmat*bmat;

        % 计算rhsk  包括相场造成的力，以及体积力和面力
        rhsk=rhsk-detjac*bmat'*estres(j,1:3)';

        % 计算delemphi
        delemphi=delemphi+detjac*(gc*xlc*(ndx'*ndx)+(gc/xlc+2.0*psi)*(nij'*nij));

        % 计算rhsp
        rhsp=rhsp-detjac*(gc*xlc*(ndx'*(ndx*ephi))+((gc/xlc+2.0*psi)*phi-2.0*psi)*nij');
    end
end