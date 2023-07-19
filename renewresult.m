% 本程序用于更新力场以及相场计算值
function [totalph,totalstress,totaldisp,totalhm,loadforcest,pristress,pristrain]=renewresult ...
         (snode,degssor,xdetaf,xdetap,sudisp,totaldisp,totalph,totalhm,element,dmat,guassnode, ...
         totbmart,nnode,linknode,totalstress,foreleminf,lemada,meu,modulus,poisson)
    % 提取节点位移和相场 udisp存储位移，resphi存储相场值
    udisp=zeros(snode*2,1);
    udisp(degssor)=xdetaf;
    udisp=udisp+sudisp';

    % 存储每步的节点位移和相场增量
    ddisp(1:snode,1)=udisp(2*(1:snode)-1);
    ddisp(1:snode,2)=udisp(2*(1:snode)); 
    totaldisp=totaldisp+ddisp;
    totalph=totalph+xdetap;

    % 按照单元存储节点与单元坐标
    [dispx,dispy]=elemdispcoor(ddisp,element);

    % 已知位移推求应变 ---=为方便后续采用高斯点计算，故此处采用高斯点应变推求节点应变
    [strain,stress]=stresstrainop2(dmat,dispx,dispy,guassnode,totbmart,nnode);

    % 求取节点平均应力应变值,avestrain按照节点存储
    [~,avestress]=avestressfunop(element,strain,stress,linknode,snode,nnode);

    % 在phi和1之间取小值
    totalph(1:snode)=min(totalph(1:snode),1);

    % 形成总应力和应变矩阵
    % totalstrain=totalstrain+avestrain;
    totalstress=totalstress+avestress;

    % 通过应变计算反力: 节点，y方向，荷载步，    长度，     force=应力乘以面积
    loadforcest(1)=sum(totalstress(foreleminf(:,2),1).*foreleminf(:,3)/2/1000);
    loadforcest(2)=sum(totalstress(foreleminf(:,2),2).*foreleminf(:,3)/2/1000);

    % 记录phi值
    % totalphi(:,st1)=totalph;

    % 该段程序用于求解主应力和主应变
    [pristress]=pricipstress(snode,totalstress);
    [pristrain]=pricipstrain3(snode,pristress,modulus,poisson);

    %% 拉压分解，更新hm矩阵，使得每一个节点的phi值均取为最大
    % 进行拉压分解
    prisum=pristrain(:,1)+pristrain(:,2);
    % 进行麦考利括号的计算
    pristrain(pristrain<0)=0;
    prisum(prisum<0)=0;
    allhm=lemada/2.*prisum.^2+meu.*(pristrain(:,1).^2+pristrain(:,2).^2);
    % 选择历史大值
    totalhm=max(totalhm,allhm);
end