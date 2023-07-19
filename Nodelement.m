
%% 该程序用于提取Abaqus导出的inp文件中的节点坐标及单元编号并输出txt文件
%% 由李汉章，徐文强基于Emilio转化代码编写

clear variables
clc

% Open the file
fid = fopen('0228-30.inp','r');
%打开文件
file=textscan(fid,'%s','Delimiter','\n');
%读取文件
flines=(file{1});
%文件分行
lif=length(flines); 
%存储文件行数
FFLINES=upper(flines);
%文件中所有字符大写
Lin_nodele_cells=strfind(FFLINES,['*' upper('Node')]); 
%搜索出现*Node行
Lin_nodele_zeros=cellfun(@isempty,Lin_nodele_cells); 
%判定Lin_nodele_cells各行是否为空数组
start_nodele=find(Lin_nodele_zeros==0); 
%查找非零元素

Lin_element_cells=strfind(FFLINES,['*' upper('Element')]); 
%搜索出现*ELEMENT行
Lin_element_zeros=cellfun(@isempty,Lin_element_cells); 
%判定Lin_element_cells各行是否为空数组
start_elements=find(Lin_element_zeros==0); 
%查找非零元素

Lin_order_cells=strfind(FFLINES,'*');
%搜索出现*行
Lin_order_zeros=cellfun(@isempty,Lin_order_cells); 
%判定Lin_order_cells各行是否为空数组
Lin_orders=find(Lin_order_zeros==0);
%查找非零元素

for i=1:length(start_nodele)
 end_nodele(i)=Lin_orders(find(Lin_orders==start_nodele(i))+1);
 %i=1:4
 %当Lin_orders数组与start_elements数组元素值相同时，Lin_orders相同元素下一元素赋值给 end_nodele
 p=2;
 str_a=FFLINES{start_nodele(i)}(~isspace(FFLINES{start_nodele(i)}));
 %消除空格
 flag=0;
 if length(str_a)==(length('Node')+1)
  flag=1;    
 else
  if str_a(length('Node')+2)==','
   flag=1;
  end
 end
 %以上几行作用为确定是否为节点编号上一行

 if flag==1
  for b=(start_nodele(i)+1):(end_nodele(i)-1) 
 %b为节点坐标编号起止行
   Node{i}(p,:)=strread(FFLINES{b}, '%f', 'delimiter', ',');
   p=p+1;
  end
 end
end

LE=Node{1}; 
[NUME, ~] = size(LE);
Nnode = NUME-1;
LE1=zeros(NUME,5);
LE1(:,1:3)=LE(:,1:3);
%计算节点数及节点坐标数组

%Order=floor(log10(NUME));
%offset=10*10^Order;
%LEnew=zeros(NUME,5);
%end    

for i=1:length(start_elements)
 end_elements(i)=Lin_orders(find(Lin_orders==start_elements(i))+1);
 %i=1:2
 %当Lin_orders数组与start_elements数组元素值相同时，Lin_orders相同元素下一元素赋值给 end_elements
 q=1;
 str_b=FFLINES{start_elements(i)}(~isspace(FFLINES{start_elements(i)}));
 %消除空格
 flag=0;
 if length(str_b)==(length('Element')+1)
  flag=1;    
 else
  if str_b(length('Element')+2)==','
   flag=1;
  end
 end
 %以上几行作用为确定是否为单元编号上一行

 if flag==1
  for b=(start_elements(i)+1):(end_elements(i)-1) 
 %b为单元编号起止行
   Nodele{i}(q,:)=strread(FFLINES{b}, '%f', 'delimiter', ',');
   q=q+1;
  end
 end
end

LEnew=Nodele{1};
[NELE,~]=size(LEnew);
Linetotal=NUME+NELE;

LEtotal=zeros(Linetotal,5);
LEtotal(1:NUME,:)=LE1;
LEtotal(NUME+1:Linetotal,:)=LEnew;
LEtotal(1,1)=NUME-1;LEtotal(1,2)=NELE;


% 建立计算输出文件
tonoel=[NUME-1 NELE];
LE2=LE(2:length(LE),:);
fileid=fopen('0228-30.txt','wt');
formatSpec=' %5d  %5d\n';
fprintf(fileid,formatSpec,tonoel');
formatSpec='%5d   %12.7f  %12.7f\n';
fprintf(fileid,formatSpec,LE2');
formatSpec='%5d  %6d   %6d   %6d   %6d\n';
fprintf(fileid,formatSpec,LEnew');

fclose all;