% 本程序用于调试自适应网格改进，绿单元转化为红单元
function [nodecoorc,elementc,elemark,elemindex,nodeindex2,nodenum,redelem12,greenelem,grecord,rednode] ...
    =adapmeshstz(itnodecoor,itelement,nodecoor,element,elemark,elemindex,totalph,crit)

    %% 输入单元和节点信息
    itsnode=size(itnodecoor,1);
    itselement=size(itelement,1);
    nnode=4; % 四节点单元
    snode=size(nodecoor,1);
    %% 在初始网格中寻找所有的红单元
    % 上一步所有的红单元
    redelem0=find(elemark==1);
    redlast=find(elemark==1);
    rednode1=unique(element(redelem0,:));
    startelem=max(redelem0);
    redelem0=redelem0(redelem0<=itselement);
    
    % 绿单元转化的红单元
    redstage_1=find(elemark==2);
    redstage_2=find(elemark==3);
    redstage=[redstage_1;redstage_2];
    redelem1=zeros(size(redstage,1),1);
    num1=0;
    for i1=1:size(redstage,1)
        if max(totalph(element(redstage(i1),:)))>=crit
            num1=num1+1;
            redelem1(num1)=elemindex(redstage(i1));
        end
    end
    redelem1=unique(redelem1);
    redelem1(redelem1==0)=[];

    % 零单元变化为红单元
    redstage2=find(elemark==0);
    redelem2=zeros(size(redstage2,1),1);
    num2=0;
    for i1=1:size(redstage2,1)
        if max(totalph(element(redstage2(i1),:)))>=crit
            num2=num2+1;
            redelem2(num2)=elemindex(redstage2(i1));
        end
    end
    redelem2=unique(redelem2);
    redelem2(redelem2==0)=[];

    % 总的红单元为三者之和
    redelemz=[redelem0; redelem1; redelem2];
    redelemz=unique(redelemz);

    % 所有红单元
    elemredz=itelement(redelemz,:);
    % 定义不重复的红单元所对应的点
    unipoi=unique(elemredz);


    %% 依据总的红单元寻找总的绿单元
    % 总单元中去除红单元
    elementr=itelement;
    elementr(redelemz,:)=0;

    % 定义边缘单元
    edgelem=[];
    for i=1:nnode
        row2=ismember(elementr(:,i),unipoi);
        index1=find(row2);
        %搜索边缘单元并存储
        edgelem=[edgelem; index1];
    end
    % 保留不重复边缘单元
    edgelem=unique(edgelem);
    % 找到所有边单元
    edgelem(edgelem==0)=[];

    % 搜索边缘单元中非红单元的单元，即为绿单元
    [row3,~]=find(ismember(edgelem,redelemz)==0);
    % 筛选得到绿单元
    greenelem=edgelem(row3);
    greenelem=unique(greenelem);
    % elemgreen 存储所有绿单元的节点
    elemgreen=itelement(greenelem,:);

    %% 排除一种特殊情况，当绿单元中出现四个红点，则判断为红单元
    % 判断需要更改的单元
    exnum=0;
    cmax=length(greenelem);
    changelem=zeros(cmax,1);
    % 记录绿单元特征
    grecord=zeros(length(greenelem),4);
    for i=1:size(greenelem,1)
        rgpoi=ismember(itelement(greenelem(i),:),unipoi);
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
    %  红单元细化的范围为：1.绿单元转化的红单元，2.零单元转化的红单元
    redelem12=[redelem1;redelem2];
    elemred=itelement(redelem12,:);
    rednum=size(redelem12,1);
    newelem=zeros(8*rednum,4);  %newly formed elements with nodes定义新单元
    % 建立新单元的索引
    elemindexc=zeros(8*rednum,1);
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
        element(redelem12(i),:)=[elemred(i,1) snode+12*(i-1)+1 snode+12*(i-1)+9 snode+12*(i-1)+8];
        % allocate new elements 为新单元分配节点
        newelem(8*(i-1)+1,:) =[snode+12*(i-1)+1 snode+12*(i-1)+2 snode+12*(i-1)+10 snode+12*(i-1)+9];
        newelem(8*(i-1)+2,:) =[snode+12*(i-1)+2 elemred(i,2) snode+12*(i-1)+3 snode+12*(i-1)+10];
        newelem(8*(i-1)+3,:) =[snode+12*(i-1)+10 snode+12*(i-1)+3 snode+12*(i-1)+4 snode+12*(i-1)+11];
        newelem(8*(i-1)+4,:) =[snode+12*(i-1)+11 snode+12*(i-1)+4 elemred(i,3) snode+12*(i-1)+5];
        newelem(8*(i-1)+5,:) =[snode+12*(i-1)+12 snode+12*(i-1)+11 snode+12*(i-1)+5 snode+12*(i-1)+6];
        newelem(8*(i-1)+6,:) =[snode+12*(i-1)+7 snode+12*(i-1)+12 snode+12*(i-1)+6 elemred(i,4)];
        newelem(8*(i-1)+7,:) =[snode+12*(i-1)+8 snode+12*(i-1)+9 snode+12*(i-1)+12 snode+12*(i-1)+7];
        newelem(8*(i-1)+8,:) =[snode+12*(i-1)+9 snode+12*(i-1)+10 snode+12*(i-1)+11 snode+12*(i-1)+12];
        elemindexc(8*(i-1)+(1:8))=redelem12(i);
    end
    nelem=size(newelem,1);
    nsnode=size(nodecoor,1);
    nodenum=nsnode-snode;
    % 对不同类型的单元进行标记，不改变单元为0，红单元为1，2个点绿单元为2，3个点绿单元为3
    % 原始单元标记
    elemark(redelem12)=1;
    % 红单元标记
    elemark(startelem+1:startelem+nelem)=1;

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
            elemindexc(nelem+(1:4))=greenelem(i);
            elemark(greenelem(i))=3;
            elemark(startelem+(nelem+1:nelem+4))=3;

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
                elemindexc(nelem+(1:3))=greenelem(i);
                elemark(greenelem(i))=2;
                elemark(startelem+(nelem+1:nelem+3))=2;
            end
            nsnode=length(nodecoor);
            nelem=length(newelem);
        end
    end

    %% 生成总单元，并对内部节点进行重新排序

    % 上一步生成的红色节点并不需要排序，因此排序内容包括这一步生成的红色节点以及绿色节点
    % 上一步的红色节点
    rednode=rednode1(rednode1>=itsnode);
    nnodecoor1=nodecoor(rednode,:);
    % 建立上一步新生成的红节点的映射
    amax=max(rednode);
    nindex=zeros(amax,1);
    nindex(rednode)=1;
    anum=itsnode;
    for i1=1:amax
        if nindex(i1)==1
            anum=anum+1;
            nindex(i1)=anum;
        end
    end

    % 对新生成的单元的节点编号进行压缩
    for i=1:size(newelem)
        for j=1:nnode
            if newelem(i,j)>snode
                newelem(i,j)=newelem(i,j)-snode+anum;
            end
        end
    end

    for i=1:itselement
        for j=1:nnode
            if element(i,j)>snode
                element(i,j)=element(i,j)-snode+anum;
            end
        end
    end

    % 生成总单元
    elementc=[element(1:startelem,:); newelem];

    % 单元节点重新生成
    for i2=1:size(redlast,1)
        for j1=1:nnode
            if elementc(redlast(i2),j1)>itsnode 
                elementc(redlast(i2),j1)=nindex(elementc(redlast(i2),j1));
            end
        end
    end


    %% 对所有新节点进行排序
    nnodecoor(:,1:2)=[nnodecoor1; nodecoor(snode+1:length(nodecoor),:)]; %形成新的节点存储矩阵
    for i=1:size(nnodecoor,1)        %新矩阵第三列编号
        nnodecoor(i,3)=itsnode+i;
        if  i<=anum-itsnode
            nnodecoor(i,4)=0;
        else
            nnodecoor(i,4)=1;
        end
    end
    % 按第一列升序排列
    nnodecoor=sortrows(nnodecoor);

    % 建立索引和新节点
    newnode=itsnode+1;
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
    a11=max(nodecoor(:,1));
    a12=min(nodecoor(:,1));
    lengthx=a11-a12;
    err=lengthx/1000;

    for i=1:length(nnodecoor)-1
        if abs(nnodecoor(i,1)-nnodecoor(i+1,1))<err && abs(nnodecoor(i,2)-nnodecoor(i+1,2))<err  ...
                && nnodecoor(i+1,4)==1
            numindex=numindex+1;
            nodeindex(numindex,:)=[nnodecoor(i+1,3) newnode];
        else
            numindex=numindex+1;
            newnode=newnode+1;
            nodeindex(numindex,:)=[nnodecoor(i+1,3) newnode];
            nodecoornew(newnode-itsnode,:)=nnodecoor(i+1,1:2);
            nodeindex2(newnode-itsnode,:)=[nnodecoor(i+1,3) newnode];
        end
    end

    nodecoornew(newnode-itsnode+1:newmax,:)=[];
    nodeindex2(newnode-itsnode+1:newmax,:)=[];

    % 生成新的节点
    nodecoorc=[itnodecoor;nodecoornew];
    nodeindex=sortrows(nodeindex,1);

    % 更新初始单元矩阵elementc
    for i=1:size(elementc,1)
        for j=1:nnode
            if elementc(i,j)>itsnode
                elementc(i,j)=nodeindex(elementc(i,j)-itsnode,2);
            end
        end
    end
    
    % 更新单元索引矩阵
    elemindex=[elemindex(1:startelem); elemindexc];
end
