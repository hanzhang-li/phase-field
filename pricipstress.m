% 该段程序用于求解主应力
function [pristress]=pricipstress(snode,a1)
pristress=zeros(snode,2);
for i=1:snode
    pristress(i,1)=(a1(i,1)+a1(i,2))/2+sqrt((a1(i,1)-a1(i,2))^2/4+a1(i,3)^2);
    pristress(i,2)=(a1(i,1)+a1(i,2))/2-sqrt((a1(i,1)-a1(i,2))^2/4+a1(i,3)^2);    
end