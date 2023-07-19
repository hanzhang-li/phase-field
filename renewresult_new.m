% 本程序用于更新力场以及相场计算值
function [totalph,totalstress,totaldisp,totalhm,loadforcest,totalstrain,Cmat,tstress_p,tstrain_en] ...
    =renewresult_new(snode,degssor,xdetaf,xdetap,sudisp,totaldisp,totalph,totalhm,element, ...
    totalstrain,guassnode,totbmart,nnode,linknode,foreleminf,lemada,meu,xk)
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
    [strain]=nodestrain(dispx,dispy,guassnode,totbmart,nnode);

    % 求取节点平均应变值,avestrain按照节点存储  
    [avestrain]=ave_strain(element,strain,linknode,snode,nnode);
    totalstrain=totalstrain+avestrain;
    
    % 在phi和1之间取小值
    totalph(1:snode)=min(totalph(1:snode),1);
    
    % 依据总应变以及总相场值求解每一个节点的应力（显示应力、计算应力），C（应力应变关系），应变能
    [Cmat,tstress_p,totalstress,tstrain_en,totalhm2]=phi_stress(totalstrain,totalph,snode,lemada,meu,xk);

    % 拉应变能选择历史大值
    totalhm=max(totalhm,totalhm2);
    
    % 通过应变计算反力: 节点，y方向，荷载步，    长度，     force=应力乘以面积
    loadforcest(1)=sum(totalstress(foreleminf(:,2),1).*foreleminf(:,3)/2/1000);
    loadforcest(2)=sum(totalstress(foreleminf(:,2),2).*foreleminf(:,3)/2/1000);

end