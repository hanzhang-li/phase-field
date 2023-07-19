% 本程序建立高斯点与节点对应矩阵
function [guassnode]=guasslinknode(gpt)
a1=1+gpt/2;
a2=1-gpt/2;
a3=-0.5;
guassnode=[a1 a3 a2 a3; a3 a1 a3 a2; a2 a3 a1 a3; a3 a2 a3 a1];