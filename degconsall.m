% 本程序用于添加固定约束，x1,x2为横坐标，y1，y2为纵坐标，a b为约束为0
% 时表示约束，为1时为自由
function [degfree]=degconsall(posdeg,degfree)
for i=1:length(posdeg)
    if posdeg(i,2)==0
        degfree(posdeg(i,1),1)=0;
    end
    if posdeg(i,3)==0
        degfree(posdeg(i,1),2)=0;
    end
end