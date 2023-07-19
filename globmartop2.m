function [totbmart,tdetmart,tndxmart]=globmartop2(tndjmart,elementx,elementy)
[selem,selen]=size(elementx);
totbmart=zeros(12*selem,8);
tdetmart=zeros(4*selem,1);
tndxmart=zeros(8*selem,selen);

for i=1:selem
    %建立B矩阵即应变与位移关系矩阵
    dnodecoor(:,1)=elementx(i,:);
    dnodecoor(:,2)=elementy(i,:);
    for j=1:selen
        %建立B矩阵即应变与位移关系矩阵
        [detjac,b,ndx]=pstrainlinkdisp(tndjmart(2*(j-1)+(1:2),:),dnodecoor);
        % 将b矩阵装载至totbmart矩阵中
        totbmart(12*(i-1)+3*(j-1)+(1:3),:)=b;
        tdetmart(4*(i-1)+j)=detjac;
        tndxmart(8*(i-1)+2*(j-1)+(1:2),:)=ndx;
    end
end