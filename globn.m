% 本程序用于提前计算每个单元的nij,ndj矩阵
function [tnijmart,tndjmart]=globn(gaur,gaus)
tnijmart=zeros(4,4);
tndjmart=zeros(4*2,4);
for j=1:4
    % 建立形函数以及形函数对于局部坐标的导数，采用三维等参单元
    [nij,ndj]=nxhs(gaur(j),gaus(j));
    tnijmart(j,:)=nij;
    tndjmart(2*(j-1)+(1:2),:)=ndj;
end
