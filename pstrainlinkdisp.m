% 本程序用于建立B矩阵
function [detjac,b,ndx]=pstrainlinkdisp(ndj,dnodecoor)
        % 建立雅克比矩阵
        jacob=ndj*dnodecoor;
        % 雅克比矩阵行列式
        detjac=det(jacob);
        % 形函数对于整体坐标导数
        ndx=jacob\ndj;
        % 建立B矩阵
        b=zeros(3,8);
        for k=1:4
            b(1,2*k-1)=ndx(1,k);b(2,2*k)=ndx(2,k);
            b(3,2*k-1)=ndx(2,k);b(3,2*k)=ndx(1,k);
        end