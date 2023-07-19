% 本程序用于计算所选单元的长度
function [elemleng,foreleminf]=elemlength(choelem,nodecoor,a)
elemleng=[];
foreleminf=zeros(length(choelem),3);
foreleminf(:,1:2)=choelem;
for i=1:length(choelem)/2
    if a==1
         leng=abs(nodecoor(choelem(2*i-1,2),a)-nodecoor(choelem(2*i,2),a));
         elemleng=[elemleng; choelem(2*i,1),leng];
         foreleminf(2*i-1,3)=leng;
         foreleminf(2*i,3)=leng;
    end
end


