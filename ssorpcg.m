% 本程序采用 预处理共轭梯度法 求解大型稀疏矩阵 ssor_pcg方法
function [x,k]=ssorpcg(sk,fload,ma,id,sdegfree,xinitial)
maxk=8000;    % 最大迭代次数
zero2=1.0e-30;  % 误差限
w=1;        % 超松弛系数
w1=(2-w)/w;
k=0;   % 记录循环次数

label=1;
v=zeros(1,sdegfree); z0=zeros(1,sdegfree);
for i=1:sdegfree
    v(i)=sk(ma(i))*w1;
end
x=xinitial;
% 赋初值，采用前一步的温度值，作为此步温度迭代的初值
[temp1]=initvalue(sk,xinitial,ma,id,sdegfree,1);
g0=temp1-fload;
[y0]=inverse(g0,w,sk,ma,id,sdegfree,label);
for i=1:sdegfree
    z0(i)=-v(i)*y0(i);
end
[d0]=inverse(z0,w,sk,ma,id,sdegfree,0);
toler=0;
for i=1:sdegfree
    toler=toler+y0(i)*y0(i)*v(i);
end
for i=1:maxk
    c=0;
    for j=1:sdegfree
        c=c+d0(j)*(2*z0(j)-d0(j)*v(j));
    end
    tk=toler/c;
    for j=1:sdegfree
        x(j)=x(j)+tk*d0(j);
    end
    if (abs(toler)<=zero2)
        return;
    end
    for j=1:sdegfree
        temp1(j)=z0(j)-d0(j)*v(j);
    end
    [temp2]=inverse(temp1,w,sk,ma,id,sdegfree,1);
    for j=1:sdegfree
        y0(j)=y0(j)+tk*(temp2(j)+d0(j));
    end
    toler1=0;
    for j=1:sdegfree
        toler1=toler1+y0(j)*y0(j)*v(j);
    end
    beta=toler1/toler;
    for j=1:sdegfree
        z0(j)=-y0(j)*v(j)+beta*z0(j);
    end
    [d0]=inverse(z0,w,sk,ma,id,sdegfree,0);
    toler=toler1;
    k=k+1;
end
k=k-1;
        
