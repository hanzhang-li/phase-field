% 本程序用于求解变带宽存储时：A*x0
function [b]=initvalue(sk,x,ma,id,sdegfree,label)
% 主对角元素相乘
b=zeros(sdegfree,1);

for i=1:sdegfree
  b(i)=sk(ma(i))*x(i);
end
if (label==1)     % 下三角矩阵
    for i=2:sdegfree
        low=ma(i-1)+1; 
        up=ma(i)-1;
        temp1=x(i);
        for j=low:up
            bj=id(j);
            b(i)=b(i)+sk(j)*x(bj);
            b(bj)=b(bj)+sk(j)*temp1;
        end
    end
else
    for i=1:sdegfree-1
        low=ma(i)+1;
        up=ma(i+1)-1;
        temp1=x(i);
        for j=low:up
            bj=id(j);
            b(i)=b(i)+sk(j)*x(bj);
            b(bj)=b(bj)+sk(j)*temp1;
        end
    end
end
    
