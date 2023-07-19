% 本程序用于展现荷载位移曲线 
clc,clear
% 输入路径名称
Path='combine_II_crack_sb\';
File=dir(fullfile(Path,'荷载-位移.dat'));
struct=[Path  File.name];
Tval0=readtable(struct);
Tval=table2array(Tval0);
loaddisp=Tval(:,1:2);
loadforce=Tval(:,3:4);
plot(loaddisp(:,2),loadforce(:,2));