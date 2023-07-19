% 本程序由于求解高斯点应变，同时采用高斯点应变推求节点应变
function [strain]=nodestrain(dispx,dispy,guassnode,totbmart,nnode)
% 数组初始化
[selem,selen]=size(dispx);
elementdisp=zeros(1,8); 
strain=zeros(nnode*selem,3); 
for i=1:selem
    % 单元节点的x，y方向的位移
    elementdisp(2*(1:4)-1)=dispx(i,1:4);
    elementdisp(2*(1:4))=dispy(i,1:4);
    estrain2=totbmart(12*(i-1)+(1:12),:)*elementdisp';
    % 建立高斯点存储
    gstrain=zeros(4,3);
    for j=1:selen
        % 单元节点高斯点的应变strain=b*disp
        estrain=estrain2(3*(j-1)+(1:3),:);
        % 三维矩阵存储,第一维为x方向，第二维为y方向，第三维为切应力，gstrain为高斯点应变，gstress为高斯点应力
        gstrain(j,:)=estrain(:)';
    end
    % 利用elemn外推节点的应变，strain为节点应变，stress为节点应力
    dstrain=guassnode*gstrain;
    strain(4*(i-1)+(1:4),:)=dstrain;
end