% 本子程序用于生成新的网格节点,并嵌入到主程序中
function [nodecoorc,elementc,elemark,nodeindex2,nodenum,redelem,greenelem,grecord] ...
    =adapmeshst2(nodecoor,element,snode,selement,nnode,totalph,crit,elemark,nodemark_0)
    intnodecoor=nodecoor;
    nodenum(1)=snode;
    %% 以下为判定红单元程序
    % 对每一种单元进行标记识别
    % 设置最大节点数
    nodemax=length(nodemark_0);
    % redpoi为红点矩阵
    redpoi=zeros(nodemax,1); 
    rednum=0;
    % 网格剖分条件
    for i=1:nodemax
        % 设置红点识别条件，即当前步达到φ的阈值，但上一步未达到
        if totalph(nodemark_0(i))>=crit
            rednum=rednum+1;
            % 找到红点并存储
            redpoi(rednum)=nodemark_0(i);
        end
    end
    
    redpoi=redpoi(1:rednum);

    % 定义红单元矩阵
    redelem=[];
    for i=1:nnode    
        row=ismember(element(:,i),redpoi);
        index1=find(row);
        %搜索边缘单元并存储
        redelem=[redelem; index1];
    end

    % 删去重复单元
    redelem=unique(redelem);
    redelem(redelem==0)=[];

    % elemred is the total nodes of red elements 找到所有红单元点
    elemred=element(redelem,:);
    % 定义不重复的红单元所对应的点
    unipoi=unique(elemred);


    %% 以下查找绿单元
    % 搜索红点与红单元点不重合的点
    [row1,~]=find(ismember(unipoi,redpoi)==0);
    % 绿节点矩阵，筛选绿点，绿点连接红单元与绿单元
    greenpoi=unipoi(row1); 

    % 定义边缘单元
    edgelem=[];
    for i=1:nnode    
        row2=ismember(element(:,i),greenpoi);
        index2=find(row2);
        %搜索边缘单元并存储
        edgelem=[edgelem; index2];
    end
    
    % 保留不重复边缘单元
    edgelem=unique(edgelem);
    % 找到所有边单元
    edgelem(edgelem==0)=[];
    % 搜索边缘单元中非红单元的单元，即为绿单元
    [row3,~]=find(ismember(edgelem,redelem)==0);
    % 筛选得到绿单元
    greenelem=edgelem(row3);
    greenelem=unique(greenelem);
    % elemgreen 存储所有绿单元的节点
    elemgreen=element(greenelem,:);

    %% 排除一种特殊情况，当绿单元中出现四个红点，则判断为红单元
    % 判断需要更改的单元
    exnum=0;
    rednum=size(redelem,1);
    cmax=length(greenelem);
    changelem=zeros(cmax,1);
    % 记录绿单元特征
    grecord=zeros(length(greenelem),4);
    for i=1:size(greenelem,1)
        rgpoi=ismember(element(greenelem(i),:),greenpoi);
        grecord(i,1:4)=rgpoi;
        if sum(rgpoi)==4
            exnum=exnum+1;
            changelem(exnum)=i;
       end
    end
%     changelem(exnum+1:cmax)=[];
%     % 对红绿单元进行重新编辑
%     if exnum~=0
%         redelem=[redelem; changelem];
%         greenelem(changelem)=[];
%         grecord(changelem,:)=[];
%         redelem=sort(redelem);
%         
%         elemred=element(redelem,:);
%         elemgreen=element(greenelem,:);     
%     end

    %% 对红绿单元进行细化
    %   此处显示详细说明
    newelem=zeros(8*rednum,4);  %newly formed elements with nodes定义新单元
    elemred(:,5:8)=elemred(:,1:4);

    %% 以下对红单元进行细化
    for i=1:size(elemred,1)
        for j=1:nnode
        % allocate coordinate for newly formed nodes  below为红单元新节点赋予坐标
        nodecoor(snode+12*(i-1)+2*j-1,:)= (2*nodecoor(elemred(i,j),:)+nodecoor(elemred(i,j+1),:))/3;
        nodecoor(snode+12*(i-1)+2*j,:)  = (2*nodecoor(elemred(i,j+1),:)+nodecoor(elemred(i,j),:))/3;
        nodecoor(snode+12*(i-1)+8+j,:)  = (2*nodecoor(elemred(i,j+1),:)+2*nodecoor(elemred(i,j+3),:)+...
                                           4*nodecoor(elemred(i,j),:)+nodecoor(elemred(i,j+2),:))/9;   
        end
        % 以下为新单元赋予节点
        % refill the existing elements保留原有单元号
        element(redelem(i),:)=[elemred(i,1) snode+12*(i-1)+1 snode+12*(i-1)+9 snode+12*(i-1)+8];
        % allocate new elements 为新单元分配节点
        newelem(8*(i-1)+1,:) =[snode+12*(i-1)+1 snode+12*(i-1)+2 snode+12*(i-1)+10 snode+12*(i-1)+9];
        newelem(8*(i-1)+2,:) =[snode+12*(i-1)+2 elemred(i,2) snode+12*(i-1)+3 snode+12*(i-1)+10];
        newelem(8*(i-1)+3,:) =[snode+12*(i-1)+10 snode+12*(i-1)+3 snode+12*(i-1)+4 snode+12*(i-1)+11];
        newelem(8*(i-1)+4,:) =[snode+12*(i-1)+11 snode+12*(i-1)+4 elemred(i,3) snode+12*(i-1)+5];
        newelem(8*(i-1)+5,:) =[snode+12*(i-1)+12 snode+12*(i-1)+11 snode+12*(i-1)+5 snode+12*(i-1)+6];
        newelem(8*(i-1)+6,:) =[snode+12*(i-1)+7 snode+12*(i-1)+12 snode+12*(i-1)+6 elemred(i,4)];
        newelem(8*(i-1)+7,:) =[snode+12*(i-1)+8 snode+12*(i-1)+9 snode+12*(i-1)+12 snode+12*(i-1)+7];
        newelem(8*(i-1)+8,:) =[snode+12*(i-1)+9 snode+12*(i-1)+10 snode+12*(i-1)+11 snode+12*(i-1)+12]; 
    end
    nelem=size(newelem,1);
    nsnode=size(nodecoor,1);
    nodenum(2)=nsnode-snode;
    % 对不同类型的单元进行标记，不改变单元为0，红单元为1，2个点绿单元为2，3个点绿单元为3
    % 原始单元标记
    elemark(redelem)=1;
    % 红单元标记
    elemark(selement+1:selement+nelem)=1;

    %% 以下对绿单元细化
    % greypoi 定义灰单元，当length(find(rgpoi==1))==3，灰单元为非红单元点，当length(find(rgpoi==1))==2，为红单元点
    elemgreen(:,5:8)=elemgreen(:,1:4);
    for i=1:size(elemgreen,1)
        rgpoi=grecord(i,:);
         if sum(rgpoi)==3
            greypoi=find(rgpoi==0);
            for j=1:2
            % allocate coordinate for newly formed nodes  below为绿单元新节点赋予坐标
            nodecoor(nsnode+2*j-1,:)= (2*nodecoor(elemgreen(i,greypoi+j),:)+nodecoor(elemgreen(i,greypoi+j+1),:))/3;
            nodecoor(nsnode+2*j,:)  = (2*nodecoor(elemgreen(i,greypoi+j+1),:)+nodecoor(elemgreen(i,greypoi+j),:))/3;
            end
            nodecoor(nsnode+5,:)  = (2*nodecoor(elemgreen(i,greypoi+1),:)+2*nodecoor(elemgreen(i,greypoi+3),:)+...
                                              4*nodecoor(elemgreen(i,greypoi),:)+nodecoor(elemgreen(i,greypoi+2),:))/9;   
            nodecoor(nsnode+6,:)  = (2*nodecoor(elemgreen(i,greypoi+1),:)+2*nodecoor(elemgreen(i,greypoi+3),:)+...
                                              4*nodecoor(elemgreen(i,greypoi+2),:)+nodecoor(elemgreen(i,greypoi),:))/9;  
            
            % 以下为新单元赋予节点
            % refill the existing elements保留原有单元号
            element(greenelem(i),:)=[elemgreen(i,greypoi) elemgreen(i,greypoi+1) nsnode+1 nsnode+5];
            %allocate new elements 为新单元分配节点
            newelem(nelem+1,:) =[nsnode+1 nsnode+2 nsnode+6 nsnode+5];
            newelem(nelem+2,:) =[nsnode+2 elemgreen(i,greypoi+2) nsnode+3 nsnode+6];
            newelem(nelem+3,:) =[nsnode+6 nsnode+3 nsnode+4 nsnode+5];
            newelem(nelem+4,:) =[nsnode+5 nsnode+4 elemgreen(i,greypoi+3) elemgreen(i,greypoi)];
            elemark(greenelem(i))=3;
            elemark(selement+(nelem+1:nelem+4))=3;  
            
            nsnode=length(nodecoor);
            nelem=length(newelem); 
            
         elseif  sum(rgpoi)==2
            greypoi=find(rgpoi==1);
            if greypoi(2)-greypoi(1)==3
                greypoi(1)=4;
                greypoi(2)=1;
            end
            if greypoi(2)-greypoi(1)~=2
              % allocate coordinate for newly formed nodes  below为绿单元新节点赋予坐标
              nodecoor(nsnode+1,:)= (2*nodecoor(elemgreen(i,greypoi(1)),:)+nodecoor(elemgreen(i,greypoi(2)),:))/3;
              nodecoor(nsnode+2,:)= (2*nodecoor(elemgreen(i,greypoi(2)),:)+nodecoor(elemgreen(i,greypoi(1)),:))/3;
              nodecoor(nsnode+3,:)  = (2*nodecoor(elemgreen(i,greypoi(1)),:)+2*nodecoor(elemgreen(i,greypoi(1)+2),:)+...
                                              4*nodecoor(elemgreen(i,greypoi(2)+2),:)+nodecoor(elemgreen(i,greypoi(2)),:))/9;   
              nodecoor(nsnode+4,:)  = (2*nodecoor(elemgreen(i,greypoi(2)),:)+2*nodecoor(elemgreen(i,greypoi(2)+2),:)+...
                                              4*nodecoor(elemgreen(i,greypoi(1)+2),:)+nodecoor(elemgreen(i,greypoi(1)),:))/9;   
              % 以下为新单元赋予节点
              % refill the existing elements保留原有单元号
              element(greenelem(i),:)=[elemgreen(i,greypoi(1)) nsnode+1 nsnode+3 elemgreen(i,greypoi(2)+2)];
              % allocate new elements 为新单元分配节点
              newelem(nelem+1,:) =[nsnode+1 nsnode+2 nsnode+4 nsnode+3];
              newelem(nelem+2,:) =[nsnode+2 elemgreen(i,greypoi(2)) elemgreen(i,greypoi(2)+1) nsnode+4];
              newelem(nelem+3,:) =[nsnode+4 elemgreen(i,greypoi(1)+2) elemgreen(i,greypoi(2)+2) nsnode+3];
              elemark(greenelem(i))=2;
              elemark(selement+(nelem+1:nelem+3))=2;
            end
            nsnode=length(nodecoor);
            nelem=length(newelem);
         end
    end
    nodenum(3)=length(nodecoor)-nodenum(2);

    nnodecoor(:,1:2)=nodecoor(snode+1:length(nodecoor),:); %形成新的节点存储矩阵
    for i=1:size(nnodecoor,1)        %新矩阵第三列编号
        nnodecoor(i,3)=snode+i;
    end
    % 按第一列升序排列
    nnodecoor=sortrows(nnodecoor,[1 2]);

    %% 对生成的节点和单元进行重新排序
    % 建立索引和新节点
    newnode=snode+1;
    % 假定最多生成的单元数
    newmax=size(nnodecoor,1);
    % 定义索引
    nodeindex=zeros(newmax,2);
    nodeindex2=zeros(newmax,2);
    nodeindex(1,:)=[nnodecoor(1,3),newnode];
    nodeindex2(1,:)=[nnodecoor(1,3),newnode];
    numindex=1;
    % 定义新节点
    nodecoornew=zeros(newmax,2);
    nodecoornew(1,:)=nnodecoor(1,1:2);

    for i=1:length(nnodecoor)-1
        if nnodecoor(i,1:2)==nnodecoor(i+1,1:2)
            numindex=numindex+1;
            nodeindex(numindex,:)=[nnodecoor(i+1,3) newnode];
        else
            numindex=numindex+1;
            newnode=newnode+1;
            nodeindex(numindex,:)=[nnodecoor(i+1,3) newnode];
            nodecoornew(newnode-snode,:)=nnodecoor(i+1,1:2);
            nodeindex2(newnode-snode,:)=[nnodecoor(i+1,3) newnode];
        end
    end
    nodecoornew(newnode-snode+1:newmax,:)=[];
    nodeindex2(newnode-snode+1:newmax,:)=[];

    % 生成新的节点
    nodecoorc=[intnodecoor;nodecoornew];
    nodeindex=sortrows(nodeindex,1);

    % 生成新单元elementc
    elementc=[element;newelem];
    for i=1:length(elementc)
        for j=1:4
            if elementc(i,j)>snode
                 elementc(i,j)=nodeindex(elementc(i,j)-snode,2);
            end
        end
    end
end



