% ���������ڰ��յ�Ԫ�洢��Ԫ����
function [elementx,elementy]=elemdispcoor(nodecoor,element)
[m,n]=size(element);
elementx=zeros(m,n);
elementy=zeros(m,n);
% elementy(1:m,1:n)=nodecoor(element(1:m,1:n),2);

for j=1:n
    elementx(1:m,j)=nodecoor(element(1:m,j),1);
    elementy(1:m,j)=nodecoor(element(1:m,j),2);
end


