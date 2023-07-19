% 本程序依据应变求解每一个节点的应力（显示应力、计算应力），C（应力应变关系），应变能
function [Cmart,tstress_p,totalstress,tstrain_en,totalhm]=phi_stress(totalstrain,totalph,snode,lemada,meu,k)
    % 定义初始值
    Cmart=zeros(3*snode,3);      % 节点应力应变关系
    tstress_p=zeros(snode,3);    % 总拉应力
    totalstress=zeros(snode,3);  % 总应力值
    tstrain_en=zeros(snode,1);   % 总应变能
    totalhm=zeros(snode,1);      % 拉应力应变能（hm）
    % 计算每个节点的数值
    for i=1:snode
        strain_tensor(1,1)=totalstrain(i,1);   strain_tensor(2,2)=totalstrain(i,2);
        strain_tensor(1,2)=totalstrain(i,3)/2;   strain_tensor(2,1)=totalstrain(i,3)/2;
        phase=totalph(i);
        % C--Fourth order tensor（四阶张量，应力应变关系）；strain_en--strain energy（应变能）；
        % stress_vector 应力向量；  stress_p 拉应力;  phi_p 拉应变能
        [C,stress_p,stress_vector,strain_en,phi_p]=Miehe_model(strain_tensor,phase,lemada,meu,k);    
        Cmart(3*(i-1)+(1:3),:)=C;
        totalstress(i,:)=stress_vector;
        tstress_p(i,:)=stress_p;
        tstrain_en(i)=strain_en;
        totalhm(i)=phi_p;
    end
end