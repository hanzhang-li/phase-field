% 建立弹性矩阵D，对于平面应力或平面应变问题
function [d]=flatd(modulus,poisson,flat)
if flat==0
    d=modulus/(1-poisson^2)*[1 poisson 0; poisson 1 0; 0 0 (1-poisson)/2];
else
    modulus=modulus/(1-poisson^2); 
    poisson=poisson/(1-poisson);
    d=modulus/(1-poisson^2)*[1 poisson 0; poisson 1 0; 0 0 (1-poisson)/2];
end