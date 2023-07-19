%  本程序用于计算y0=w(-1)*g0，其中w为上三角或者下三角矩阵，y0为要求的逆
function [y0]=inverse(y0,w,sk,ma,id,n,label)
if (label==1)
    y0(1)=y0(1)/sk(1)*w;
    for ni=2:n
        low=ma(ni-1)+1;
        up=ma(ni)-1;
        for j=low:up
            nj=id(j);
            y0(ni)=y0(ni)-sk(j)*y0(nj);
        end
        y0(ni)=y0(ni)/sk(ma(ni))*w;
    end
else
    for ni=n:-1:2
        y0(ni)=y0(ni)/sk(ma(ni))*w;
        low=ma(ni-1)+1;
        up=ma(ni)-1;
        for j=low:up
            nj=id(j);
            y0(nj)=y0(nj)-sk(j)*y0(ni);
        end
    end
    y0(1)=y0(1)/sk(1)*w;
end