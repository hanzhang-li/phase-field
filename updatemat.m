% 本程序用于更新计算所需的变量和值，采用线性插值法进行求解
function [updmatr]=updatemat(inimatr,redelem,greenelem,intelement,grecord,nnode,nodeindex,nodemark)
    %% 对红单元进行处理
    intsnode=max(max(intelement));
    num=size(inimatr,2);
    nmatr=zeros(nodemark(2),num);
    elemred=intelement(redelem,:);
    elemred(:,5:8)=elemred(:,1:4);
    for i=1:size(elemred,1)
        for j=1:nnode
            % allocate coordinate for newly formed nodes  below为红单元新节点赋予坐标
            nmatr(12*(i-1)+2*j-1,:)= (2*inimatr(elemred(i,j),:)+inimatr(elemred(i,j+1),:))/3;
            nmatr(12*(i-1)+2*j,:)  = (2*inimatr(elemred(i,j+1),:)+inimatr(elemred(i,j),:))/3;
            nmatr(12*(i-1)+8+j,:)  = (2*inimatr(elemred(i,j+1),:)+2*inimatr(elemred(i,j+3),:)+...
                4*inimatr(elemred(i,j),:)+inimatr(elemred(i,j+2),:))/9;
        end
    end


    %% 绿单元进行映射
    nnmatr=length(nmatr);
    elemgreen=intelement(greenelem,:);
    elemgreen(:,5:8)=elemgreen(:,1:4);
    for i=1:length(elemgreen)
        rgpoi=grecord(i,:);
         if sum(rgpoi)==3
            greypoi=find(rgpoi==0);
            for j=1:2
            % allocate coordinate for newly formed nodes  below为绿单元新节点赋予坐标
            nmatr(nnmatr+2*j-1,:)= (2*inimatr(elemgreen(i,greypoi+j),:)+inimatr(elemgreen(i,greypoi+j+1),:))/3;
            nmatr(nnmatr+2*j,:)  = (2*inimatr(elemgreen(i,greypoi+j+1),:)+inimatr(elemgreen(i,greypoi+j),:))/3;

            nmatr(nnmatr+5,:)  = (2*inimatr(elemgreen(i,greypoi+1),:)+2*inimatr(elemgreen(i,greypoi+3),:)+...
                                              4*inimatr(elemgreen(i,greypoi),:)+inimatr(elemgreen(i,greypoi+2),:))/9;   
            nmatr(nnmatr+6,:)  = (2*inimatr(elemgreen(i,greypoi+1),:)+2*inimatr(elemgreen(i,greypoi+3),:)+...
                                              4*inimatr(elemgreen(i,greypoi+2),:)+inimatr(elemgreen(i,greypoi),:))/9;  
            end
         end
         nnmatr=length(nmatr);
         if sum(rgpoi)==2
            greypoi=find(rgpoi==1);
            if greypoi(2)-greypoi(1)==3
                greypoi(1)=4;
                greypoi(2)=1;
            end
            if greypoi(2)-greypoi(1)~=2
              %allocate coordinate for newly formed nodes  below为绿单元新节点赋予坐标
              nmatr(nnmatr+1,:)= (2*inimatr(elemgreen(i,greypoi(1)),:)+inimatr(elemgreen(i,greypoi(2)),:))/3;
              nmatr(nnmatr+2,:)= (2*inimatr(elemgreen(i,greypoi(2)),:)+inimatr(elemgreen(i,greypoi(1)),:))/3;
              nmatr(nnmatr+3,:)  = (2*inimatr(elemgreen(i,greypoi(1)),:)+2*inimatr(elemgreen(i,greypoi(1)+2),:)+...
                                              4*inimatr(elemgreen(i,greypoi(2)+2),:)+inimatr(elemgreen(i,greypoi(2)),:))/9;   
              nmatr(nnmatr+4,:)  = (2*inimatr(elemgreen(i,greypoi(2)),:)+2*inimatr(elemgreen(i,greypoi(2)+2),:)+...
                                              4*inimatr(elemgreen(i,greypoi(1)+2),:)+inimatr(elemgreen(i,greypoi(1)),:))/9;   
            end
         end
         nnmatr=length(nmatr);
    end

    %% 单元映射
    nmatr2=zeros(length(nodeindex),size(inimatr,2));
    for i=1:length(nodeindex)
        nmatr2(i,:)=nmatr(nodeindex(i,1)-intsnode,:);
    end
    updmatr=[inimatr;nmatr2];

end