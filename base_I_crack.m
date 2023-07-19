%% 本程序为采用C的应力应变关系矩阵，对程序进行重新编辑
% 基础算例
clc,clear;

%% 作用及功能介绍
% 本程序采用二维平面四节点单元
% 计算相场断裂模型
% copyright：Xu Wenqiang, (Li Hanzhang), Qiang sheng, Li Yu, Wang Kang

%% 输入单元和节点信息
% 本程序采用ansys剖分的网格数据，导入matlab之中
str0='base-3';   %%--%%--%% 输入计算网格名
storefile='base_I_crack_new1';  %%--%%--%% 输入存储结果的文件
str2='\';  str3='.data'; 
strfile=[str0 '.txt'];
getdata=textread(strfile);  %%--%%--%%--%% 输入
mkdir(storefile);  % 创建结果存储文件夹
% snode总节点数，selement总单元数
snode=getdata(1,1); selement=getdata(1,2); 
% nodecoor节点坐标
nodecoor=getdata(2:snode+1,2:3)/1e3;
% element单元所含的节点
element=getdata(snode+2:snode+1+selement,2:5);
nnode=4; % 四节点单元
ndofn=2; % 力场计算自由度
phdofn=1; % 相场计算自由度
dnum=nnode*ndofn; % 力场单元节点
% 记录初始剖分
itnodecoor=nodecoor;  itelement=element;
% 输入信息
fea=max(nodecoor); fea2=min(nodecoor);
lengx=fea(1);   lengxs=fea2(1);
lengy=fea(2);   lengys=fea2(2);
maxk=2000;


%% 输入材料属性
lemada=121.154e9; % (N/m^2);      
meu=80.769e9;  % (N/m^2); 
poisson= lemada/(2*(meu+lemada)); %%--%%--%%--%% 输入泊松比
modulus=2*meu*(1+poisson);  %(Pa)  %%--%%--%%--%% 输入弹性模量 
dens=2300;  % (Kg/m^3) %%--%%--%%--%% 输入密度 
flat=0;   %%--%%--%%--%% 输入平面应力（0）还是应变（其他）
xlc=7.5e-6; % (m) 裂缝宽度 
gc=2.7e3;  % (N/m) 裂缝断裂能
xk=0; % 模型计算参数，小值
stresscirt=9/16*sqrt(modulus*gc/(3*xlc));
straincirt=sqrt(gc/(3*xlc*modulus));
cc=sqrt(stresscirt/(modulus*straincirt));
coffe=150;
cirthm=0.5*stresscirt*straincirt*coffe;
phi0=1;

% 建立弹性矩阵D，对于平面应力问题
[dmat]=flatd(modulus,poisson,flat);

%% 裂纹尖端处重新剖分
% 选择是否进行裂纹的重新剖分
adapchose=0;  
if adapchose==1
    elemark=zeros(selement,1);
    %%--%%--%%--%% 输入裂纹端点 0为主剖分端点，1为副端点
    crkx0=lengx/2; crky0=lengy/2; crkx1=0; crky1=lengy/2;
    % 进行网格剖分，剖分为自由剖分,elemark 标记单元类型
    [nodecoor,element,elemark,snode,selement,numb,nodets,elemindex] ...
        =adapmeshjd(nodecoor,element,snode,selement,nnode,crkx0,crky0,crkx1,crky1,4*xlc,elemark);
    % 更新网格单元
end

%% 输入约束信息
% 定义信息
% 自由度指示矩阵(degfree)若是自由度等于0则为约束
degfree=ones(snode,2);
% 定义相场自由度
phdegfree=snode;
% 定义相场指示矩阵
phissor=(1:snode)';
% 记录荷载
dispnum=0;
% 存储荷载，设置最大荷载为100个,前5个为坐标和形式，后2个为方向
dispmax=100;
dispxy=zeros(dispmax,7);

% 下端ux,uy等于0
x1=lengxs; x2=lengx; y1=lengys; y2=lengys;  %%--%%--%%--%% 输入
u1=1; u2=0; %%--%%--%%--%% 输入 0为约束 其他为未约束
dispnum=dispnum+1; dispxy(dispnum,:)=[x1 x2 y1 y2 0 u1 u2];

% 输入信息
% 左侧ux等于0
x1=lengxs; x2=lengxs; y1=lengys; y2=lengy;  %%--%%--%%--%% 输入
u1=1; u2=1; %%--%%--%%--%% 输入 u1为x方向，u2为y方向，0为约束 其他为未约束
dispnum=dispnum+1; dispxy(dispnum,:)=[x1 x2 y1 y2 0 u1 u2]; 

% 输入信息
% 右侧ux等于0
x1=lengx; x2=lengx; y1=lengys; y2=lengy;  %%--%%--%%--%% 输入
u1=1; u2=1; %%--%%--%%--%% 输入 u1为x方向，u2为y方向，0为约束 其他为未约束
dispnum=dispnum+1; dispxy(dispnum,:)=[x1 x2 y1 y2 0 u1 u2];

% 输入信息
% 上端ux等于0
x1=lengxs; x2=lengx; y1=lengy; y2=lengy;  %%--%%--%%--%% 输入
u1=1; u2=1; %%--%%--%%--%% 输入 u1为x方向，u2为y方向，0为约束 其他为未约束
dispnum=dispnum+1; dispxy(dispnum,:)=[x1 x2 y1 y2 0 u1 u2];

% 施加固定约束
dispxy(dispnum+1:dispmax,:)=[];
[posdeg]=poscoorall(dispxy,nodecoor,snode);
[degfree]=degconsall(posdeg,degfree);


%% 输入荷载信息  为计算方便:重力，集中力，面力均不考虑，仅考虑位移荷载
% 存储施加位移荷载的单元
dispelem=[];
% 存储施加荷载的主要信息：单元序号（循环用），连接节点编号（荷载施加用）
% 施载节点，对应的两个局部坐标
disprelat=[];
% 存储主变化信息
dispmpiv=[];

% 记录荷载
loadnum=0;
% 存储荷载，设置最大荷载为100个,前5个为坐标和形式，后2个为方向
loadxy=zeros(dispmax,7);
% 选择节点信息
x1=lengxs; x2=lengx; y1=lengy; y2=lengy;   %%--%%--%%--%% 输入
u1=1; u2=0; %%--%%--%%--%% 输入 u1为x方向，u2为y方向，0为施加位移荷载
loadnum=loadnum+1; loadxy(loadnum,:)=[x1 x2 y1 y2 0 u1 u2];

% 施加荷载
loadxy(loadnum+1:dispmax,:)=[];
[posdeg2]=poscoorall(loadxy,nodecoor,snode);
position2=posdeg2(:,1);
[degfree]=degconsall(posdeg2,degfree);

% 建立整体自由度矩阵和节点的对应矩阵2*i (issor)，平面问题每个节点两个自由度
[issor,sdegfree]=globsor2(degfree,snode);
% 建立总结点与自由度对应矩阵
[degssor]=soormap(sdegfree,snode,ndofn,issor);


% 输入位移荷载 位移荷载变化，1位置为x方向，2位置为y方向
% 可以对每一个点施加不同方向的荷载，若为0则默认不参与计算
step=2200;    %%--%%--%%--%% 输入计算总步数 %%--%%--%%--%% 输入计算总步数
% 荷载步，数值
loadstep=zeros(step,2); 
for st1=1:step
    if st1<=400
        loadstep(st1,1)=0;  %%--%%--%%--%% 输入x方向的荷载
        loadstep(st1,2)=1e-8; %%--%%--%%--%% 输入y方向的荷载
    elseif st1<=600 
        loadstep(st1,1)=0;  %%--%%--%%--%% 输入x方向的荷载
        loadstep(st1,2)=0.5e-8; %%--%%--%%--%% 输入y方向的荷载
    else
        loadstep(st1,1)=0;  %%--%%--%%--%% 输入x方向的荷载
        loadstep(st1,2)=0.1e-8; %%--%%--%%--%% 输入y方向的荷载      
    end    
end
% 荷载步综合 1为x方向 2为y方向
loaddisp=cumsum(loadstep(:,1:2));

%% 建立单元矩阵与整体矩阵的关联矩阵
% 建立位移荷载对应的关系矩阵
% 选取节点关联的单元
[dispelem,fpelem,fpnode]=poslinkelem2(position2,element,dispelem,nnode);
% 节点关联性矩阵
[nbandp, mbandp]=nodelemlinkdp(snode,nnode,fpelem);
% 施荷单元自由度指示矩阵dfsoor
dfssor=ones(snode,2);
[dfssor]=degconsall(posdeg2,dfssor);
% 建立施加荷载的详细信息
[disprelat,dispmpiv]=dispinforfp(fpelem,nnode,ndofn,disprelat,dfssor,dispmpiv);

% 寻找每个节点的关系矩阵
[nband, mband]=nodelemlink(snode,nnode,element,selement);
% 建立二维矩阵到一维矩阵的索引，力场 nndex
[nndex,mpiv]=matndx(sdegfree,snode,ndofn,issor,nband,mband);
% 建立二维矩阵到一维矩阵的索引，相场 nndex2
[nndex2,mpiv2]=matndx(phdegfree,snode,phdofn,phissor,nband,mband);
% 建立每个单元的节点自由度矩阵
[elmssor]=selemssor(selement,nnode,ndofn,issor,element);
% 该段程序用于确定同一节点出现在几个单元之中
[linknode]=nodelink(snode,element);
% 按照单元存储单元坐标
[elementx,elementy]=elemdispcoor(nodecoor,element);

%% 求解反力位移曲线
dir1=1;  %%--%%--%%--%% 输入需要计算的单元的方向 1为x方向 2为y方向
% 求解反力-位移曲线
forcexy=loadxy(:,1:4);
% 按照节点坐标选择相应单元  choelem 1为荷载节点编号 2为整体编号
[choelem]=poselemall(forcexy,nodecoor,selement,element);
% 按照单元选择相应的长度
[elemleng,foreleminf]=elemlength(choelem,nodecoor,dir1);
% 建立计算反力位移的矩阵
loadforce=zeros(step,2);


%% 二维四节点高斯积分点
gpt=1/sqrt(3); gpt2=sqrt(3);
[gaur,gaus,wt]=guass(gpt);
% 高斯点与节点对应矩阵
[guassnode]=guasslinknode(gpt2);
% 提前计算每个单元的nij,ndj矩阵
[tnijmart,tndjmart]=globn(gaur,gaus);
% 提前计算每个单元的b,nij,det,ndx矩阵
[totbmart,tdetmart,tndxmart]=globmartop2(tndjmart,elementx,elementy);
% 建立初始步的Cmat矩阵
Cmat=zeros(3*snode,3);
for i1=1:snode
    Cmat(3*(i1-1)+(1:3),:)=dmat;
end


%% 建立输出存储矩阵
% 计算的必要输出矩阵
ddisp=zeros(snode,2);  % 存储每步的节点位移
% ddispp=zeros(snode,step);    % 存储每步的节点相场
totaldisp=zeros(snode,2);   % 存储节点总位移场
totalph=zeros(snode,1);   % 存储节点相场值
totalphpr=zeros(snode,1);  % 自适应网格存储所需
totalstrain=zeros(snode,3);   % 存储节点总应变
totalstress=zeros(snode,3);   % 存储节点总应力
totalhm=zeros(snode,1);    % 记忆矩阵  每一次需要对其更新
% totalphi=zeros(snode,step);  % 记忆每一次的phi值
% 选择性输出矩阵
oputchse=0;   %%--%%--%%--%% 如果选择为0则默认不进行输出
if oputchse~=0
     dstrain=zeros(snode,3,step); % 存储每步的节点应变
     dstress=zeros(snode,3,step); % 存储每步的节点应力
     apristress=zeros(snode,2,step);  % 存储节点主应力
     apristrain=zeros(snode,2,step);  % 存储节点主应变
end

%% 设置初始相场值
% 选择是否插入裂纹
insertchose=0;
if insertchose~=0
    % 输入裂缝端点 %%--%%--%%--%% 输入
    epx0=lengxs; epy0=0.5*(lengy-lengys)+lengys;  epx1=0.5*(lengx-lengxs)+lengxs;  epy1=0.5*(lengy-lengys)+lengys;
    [totalph,totalhm]=initialph(epx0,epy0,epx1,epy1,phi0,cirthm,snode,nodecoor,xlc,totalph,totalhm);
    % 输出初始迭代结果
    pristress=zeros(snode,3);
    strout=[storefile str2 '相场初始无迭代结果0' str3];
    outputfile(strout,snode,selement,nodecoor,element,totalph,totalhm, ...
        totaldisp,totalstress,pristress);

    %% 进行相场初始平衡计算
    % 平衡计算即力场无位移，相场存在初始迭代
    [onetk,trhsk,oneph,trhsp]=golabkmat(element,totalstress,totalph,totalhm,tnijmart,tdetmart, ...
        totbmart,tndxmart,xk,dmat,gc,xlc,elmssor,mpiv,mpiv2,nnode,dnum,sdegfree,snode,nndex,nndex2);
    % 初始力场为零
    sudisp=zeros(1,2*snode);
    msdisp=zeros(1,2*snode); 
    % 形成计算矩阵并计算
    kload=msdisp(degssor);
    trhsk=trhsk-kload;
    finitial=zeros(length(mpiv),1)+0.0000000001;
    pinitial=zeros(snode,1)+0.000000001;
    % 计算求解
    onetk=onetk'; trhsk=trhsk'; fdegfree=length(mpiv);
    % 求解力场方程
    [xdetaf,kbushuf]=ssorpcg(onetk,trhsk,mpiv,nndex,fdegfree,finitial);
    % 相场计算前准备
    oneph=oneph'; trhsp=trhsp';
    % 求解相场方程
    [xdetap,kbushup]=ssorpcgp(oneph,trhsp,mpiv2,nndex2,snode,pinitial);
    abushuit(1)=kbushuf;  abushuit(2)=kbushup;
    % 更新相场值和力场值
    [totalph,totalstress,totaldisp,totalhm,loadforcest,pristress,pristrain]=renewresult ...
             (snode,degssor,xdetaf,xdetap,sudisp,totaldisp,totalph,totalhm,element,dmat,guassnode, ...
             totbmart,nnode,linknode,totalstress,foreleminf,lemada,meu,modulus,poisson);
    % 输出初始迭代结果
    strout=[storefile str2 '相场初始迭代结果' str3];
    outputfile(strout,snode,selement,nodecoor,element,totalph,totalhm, ...
        totaldisp,totalstress,pristress);
end
    
% 记录程序运行时间
atime=zeros(4,step);
% 记录程序运行荷载步
abushu=zeros(step,2);

%% 对每一步进行循环计算并且输出计算数值
% adapstep记录自适应网格剖分次数
adapstep=0;
% adapstart判断是否需要更新网格
adapstart=0;
% 记录荷载步，以便于初始值迭代
calinitial=0;
for st1=1:step
    tic   
    %% 自适应网格更新计算所需的准备
    if adapstart==1   
        tic
        % 节点单元更新
        nodecoor=nodecoorc;      element=elementc;
        intsnode=snode;          intelement=selement;
        snode=size(nodecoor,1);  selement=size(element,1);
        
        % 网格自动剖分后结果输出
        strout=[storefile str2 'adapt_' num2str(st1-1) '步__' 'x_' num2str(loaddisp(st1-1,1)) 'mm__' ...
           'y_' num2str(loaddisp(st1-1,2)) 'mm'  str3];
        outputfile(strout,snode,selement,nodecoor,element,totalph,totalhm, ...
            totaldisp,totalstress,pristress);
        
        % 更新荷载信息
        degfree=ones(snode,2);
        [posdeg]=poscoorall(dispxy,nodecoor,snode);
        [degfree]=degconsall(posdeg,degfree);
       
        % 更新荷载位移矩阵
        [posdeg2]=poscoorall(loadxy,nodecoor,snode);
        position2=posdeg2(:,1);
        [degfree]=degconsall(posdeg2,degfree);
        
        % 更新issor以及degsoor
        [issor,sdegfree]=globsor2(degfree,snode);
        [degssor]=soormap(sdegfree,snode,ndofn,issor);
        
        % 更新关联矩阵
        dispelem=[]; disprelat=[]; dispmpiv=[];
        [dispelem,fpelem,fpnode]=poslinkelem2(position2,element,dispelem,nnode);
        [nbandp, mbandp]=nodelemlinkdp(snode,nnode,fpelem);
        dfssor=ones(snode,2);  [dfssor]=degconsall(posdeg2,dfssor);
        [disprelat,dispmpiv]=dispinforfp(fpelem,nnode,ndofn,disprelat,dfssor,dispmpiv);
        [nband, mband]=nodelemlink(snode,nnode,element,selement);
        [nndex,mpiv]=matndx(sdegfree,snode,ndofn,issor,nband,mband);
        phissor=(1:snode)';
        [nndex2,mpiv2]=matndx(snode,snode,phdofn,phissor,nband,mband);
        [elmssor]=selemssor(selement,nnode,ndofn,issor,element);
        [linknode]=nodelink(snode,element);
        [elementx,elementy]=elemdispcoor(nodecoor,element);
        [totbmart,tdetmart,tndxmart]=globmartop2(tndjmart,elementx,elementy);
        
        % 更新反力位移曲线
        [choelem]=poselemall(forcexy,nodecoor,selement,element);
        [elemleng,foreleminf]=elemlength(choelem,nodecoor,dir1);
        
        % 更新存储矩阵
        ddisp=zeros(snode,2);  % 存储每步的节点位移
        totalphpr=totalph;     % 存储更新后的相场值以便下次迭代
              
        % 矩阵更新之后，置零等待下一次更新
        adapstart=0;
        atime(4,st1)=toc;   
    end
    
    %% 建立力场和相场单元矩阵并装载至整体矩阵
    [onetk,trhsk,oneph,trhsp]=golabkmat_new(element,totalstress,totalph,totalhm,tnijmart,tdetmart, ...
        totbmart,tndxmart,Cmat,gc,xlc,elmssor,mpiv,mpiv2,nnode,dnum,sdegfree,snode,nndex,nndex2);
    % 初始化位移矩阵
    sudisp=zeros(1,2*snode);       % 力场位移矩阵
    msdisp=zeros(1,2*snode);       % 力场位移产生的荷载矩阵    
    % 生成每步的初始位移
    dispnum=length(position2);
    % loadcal 每步计算时的位移荷载
    loadcal=zeros(2*dispnum,1);
    loadcal(2*(1:dispnum)-1)=loadstep(st1,1);
    loadcal(2*(1:dispnum))  =loadstep(st1,2);
    sudisp(2*position2(1:dispnum)-1)=loadcal(2*(1:dispnum)-1);
    sudisp(2*position2(1:dispnum))  =loadcal(2*(1:dispnum));
    % 计算位移荷载对应的单元矩阵
    [msdisp]=loadmart_new(dispelem,element,nodecoor,nnode, ...
         tndjmart,Cmat,msdisp,dispmpiv,disprelat,sudisp);
    atime(1,st1)=toc;
    
    %% 形成计算矩阵并计算
    tic   
    % 形成右端项荷载矩阵:svrsh
    % 由msdisp矩阵生成计算矩阵
    kload=msdisp(degssor);  
    trhsk=trhsk-kload;
    
    % 初始求解矩阵优化
    if st1==1 || st1==calinitial+1
        finitial=zeros(length(mpiv),1)+0.0000000001;
        pinitial=zeros(snode,1)+0.0000000001;
    else
        finitial=xdetaf;
        pinitial=xdetap; 
    end
    
    % 力场计算前准备
    onetk=onetk'; trhsk=trhsk'; fdegfree=length(mpiv);  
    % 求解力场方程
    [xdetaf,kbushuf]=ssorpcg(onetk,trhsk,mpiv,nndex,fdegfree,finitial);
    
    % 相场计算前准备
    oneph=oneph'; trhsp=trhsp'; 
    % 求解相场方程
    [xdetap,kbushup]=ssorpcgp(oneph,trhsp,mpiv2,nndex2,snode,pinitial);
    
    % 记录迭代次数
    abushu(st1,1)=kbushuf;  abushu(st1,2)=kbushup;
    
    if abushu(st1,1)>=maxk
        return;
    end
    
    atime(2,st1)=toc;

    
    %% 分步结果后处理程序----获取单元和节点的位移数据
    tic
    % 更新相场值和力场值
    [totalph,totalstress,totaldisp,totalhm,loadforcest,totalstrain,Cmat,tstress_p,tstrain_en] ...
    =renewresult_new2(snode,degssor,xdetaf,xdetap,sudisp,totaldisp,totalph,totalhm,element, ...
    guassnode,totbmart,nnode,linknode,foreleminf,lemada,meu,xk);
    loadforce(st1,:)=loadforcest;
    pristress=tstress_p;
    %% 带裂纹网格自适应网格剖分,剖分条件判断，条件设置为标记为0，2，3单元的节点最大值max(phi)>=0.5 
    % 选择是否进行网格剖分
    if adapchose==1
        % 选择0单元的节点
        elemark_n1=elemark==1;
        nodemark_n1=unique(element(elemark_n1,:));
        nodeall=1:1:snode;
        nodeall(nodemark_n1)=[];
        nodemark_c=nodeall;
        % 设置剖分变量
        totalph_c=totalph(nodemark_c);  

        if  max(totalph_c)>=0.5 
           %% 将未剖分前的数据存储至问价夹中，主要输出量 节点、单元、相场值、应变能、位移值、应力、主应力
            strout=[storefile str2 num2str(st1) '步__' 'x_' num2str(loaddisp(st1,1)) 'mm__' ...
               'y_' num2str(loaddisp(st1,2)) 'mm'  str3];
            outputfile(strout,snode,selement,nodecoor,element,totalph,totalhm, ...
                       totaldisp,totalstress,pristress);

           %% 进行网格剖分 
            % 设置剖分条件，当phi值大于0.1时即进行剖分
            crit=0.1;
            % 进行网格剖分，剖分为自由剖分
            [nodecoorc,elementc,elemark,elemindex,nodeindex2,nodenum,redelem12,greenelem,grecord,rednode] ...
             =adapmeshstz(itnodecoor,itelement,nodecoor,element,elemark,elemindex,totalph,crit);

            % 进行求解值的映射，计算所需，包括totalph,totalhm,totalstress (必须量)
            [totalph]=updatematz(totalph,redelem12,greenelem,itelement,grecord,nnode,nodeindex2,nodenum,rednode);
            [totalhm]=updatematz(totalhm,redelem12,greenelem,itelement,grecord,nnode,nodeindex2,nodenum,rednode);
            [totalstress]=updatematz(totalstress,redelem12,greenelem,itelement,grecord,nnode,nodeindex2,nodenum,rednode);

            % 输出选择量
            [totaldisp]=updatematz(totaldisp,redelem12,greenelem,itelement,grecord,nnode,nodeindex2,nodenum,rednode);
            [pristress]=updatematz(pristress,redelem12,greenelem,itelement,grecord,nnode,nodeindex2,nodenum,rednode);

            % 记录网格更新次数
            adapstep=adapstep+1;
            % 更新计算所需初始准备
            adapstart=1;
            % 记录荷载步，以便于初始值计算
            calinitial=st1;
        end
    end
    
    % 最后一步结果输出
    if adapchose==1
        if st1==1 || st1==step 
            % 定义输出格式为x步_xmm位移
            strout=[storefile str2  num2str(st1) '步__' 'x_' num2str(loaddisp(st1,1)) 'mm__' ...
               'y_' num2str(loaddisp(st1,2)) 'mm'  str3];
            outputfile(strout,snode,selement,nodecoor,element,totalph,totalhm, ...
                totaldisp,totalstress,pristress);
        end
    else 
        if st1==1 || mod(st1,50)==0 || st1==step
           % 定义输出格式为x步_xmm位移
            strout=[storefile str2  num2str(st1) '步__' 'x_' num2str(loaddisp(st1,1)) 'mm__' ...
               'y_' num2str(loaddisp(st1,2)) 'mm'  str3];
            outputfile(strout,snode,selement,nodecoor,element,totalph,totalhm, ...
                totaldisp,totalstress,pristress);
        end
    end
        
    
    atime(3,st1)=toc;

end

%% 后处理及输出文件
% 输出荷载反力
strout=[storefile str2  '荷载-位移'  '.txt'];
dispforce=[loaddisp loadforce];
fileid=fopen(strout,'wt');
formatSpec='%12.9f  %12.9f  %12.6f  %12.6f \n';
fprintf(fileid,formatSpec,dispforce');
%  plot(loaddisp(1:1760,2),loadforce(1:1760,2))

aa=sum(sum(atime));
