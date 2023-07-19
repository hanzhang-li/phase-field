% 本程序用于按照单元存储单元坐标
function [elemphi]=elemdispcoorsm(totalphi,element)
[m,n]=size(element);
elemphi=zeros(m,n);
% elementy(1:m,1:n)=nodecoor(element(1:m,1:n),2);

for j=1:n
    elemphi(1:m,j)=totalphi(element(1:m,j));
end


