% 建立高斯积分点和权重函数
function [s,t,wt]=guass(gpt)
s(1)=-gpt;  t(1)=-gpt;
s(2)= gpt;  t(2)=-gpt;
s(3)= gpt;  t(3)= gpt;
s(4)=-gpt;  t(4)= gpt;
% 权重函数
wt=[1 1 1 1];  