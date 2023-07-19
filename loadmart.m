% 本程序用于计算位移荷载对应的单元矩阵
function [msdisp]=loadmart(dispelem,element,nodecoor,totalph,nnode, ...
         tnijmart,tndjmart,xk,dmat,msdisp,dispmpiv,disprelat,sudisp)  

    for i2=1:length(dispelem)
        % 单元所含节点矩阵
        elemnode3=element(dispelem(i2),:);
        % 初始单元矩阵
        dnodecoor=zeros(4,2);
        dnodecoor(1:4,:)=nodecoor(elemnode3(1:4),:);
       
        % 提取生成单元劲度矩阵所需的信息  3个总应力 3个总应变 1个phi值 1个hm值 
        ephi=totalph(elemnode3(1:4));
       
        % 计算单元劲度矩阵
        [delemk]=elemstmasm2(nnode,tnijmart,tndjmart,dnodecoor,ephi,xk,dmat);
        
        % 将位移荷载的单元矩阵装载至整体的位移荷载矩阵
        if i2==1
            for j1=1:dispmpiv(i2)
                msdisp(disprelat(j1,3))=msdisp(disprelat(j1,3))+ ...
                    delemk(disprelat(j1,4),disprelat(j1,5))*sudisp(disprelat(j1,2));
            end
        else
            for j1=dispmpiv(i2-1)+1:dispmpiv(i2)
                msdisp(disprelat(j1,3))=msdisp(disprelat(j1,3))+ ...
                    delemk(disprelat(j1,4),disprelat(j1,5))*sudisp(disprelat(j1,2));
            end
        end 
    end
end