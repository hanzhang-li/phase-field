% 该段程序用于求解主应变
function [pristrain]=pricipstrain3(snode,pristress,modulus,poisson)
pristrain=zeros(snode,2);
pristrain(:,1)=1/modulus*(pristress(:,1)-poisson*pristress(:,2));
pristrain(:,2)=1/modulus*(pristress(:,2)-poisson*pristress(:,1));   