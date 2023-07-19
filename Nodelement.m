
%% �ó���������ȡAbaqus������inp�ļ��еĽڵ����꼰��Ԫ��Ų����txt�ļ�
%% ����£�����ǿ����Emilioת�������д

clear variables
clc

% Open the file
fid = fopen('0228-30.inp','r');
%���ļ�
file=textscan(fid,'%s','Delimiter','\n');
%��ȡ�ļ�
flines=(file{1});
%�ļ�����
lif=length(flines); 
%�洢�ļ�����
FFLINES=upper(flines);
%�ļ��������ַ���д
Lin_nodele_cells=strfind(FFLINES,['*' upper('Node')]); 
%��������*Node��
Lin_nodele_zeros=cellfun(@isempty,Lin_nodele_cells); 
%�ж�Lin_nodele_cells�����Ƿ�Ϊ������
start_nodele=find(Lin_nodele_zeros==0); 
%���ҷ���Ԫ��

Lin_element_cells=strfind(FFLINES,['*' upper('Element')]); 
%��������*ELEMENT��
Lin_element_zeros=cellfun(@isempty,Lin_element_cells); 
%�ж�Lin_element_cells�����Ƿ�Ϊ������
start_elements=find(Lin_element_zeros==0); 
%���ҷ���Ԫ��

Lin_order_cells=strfind(FFLINES,'*');
%��������*��
Lin_order_zeros=cellfun(@isempty,Lin_order_cells); 
%�ж�Lin_order_cells�����Ƿ�Ϊ������
Lin_orders=find(Lin_order_zeros==0);
%���ҷ���Ԫ��

for i=1:length(start_nodele)
 end_nodele(i)=Lin_orders(find(Lin_orders==start_nodele(i))+1);
 %i=1:4
 %��Lin_orders������start_elements����Ԫ��ֵ��ͬʱ��Lin_orders��ͬԪ����һԪ�ظ�ֵ�� end_nodele
 p=2;
 str_a=FFLINES{start_nodele(i)}(~isspace(FFLINES{start_nodele(i)}));
 %�����ո�
 flag=0;
 if length(str_a)==(length('Node')+1)
  flag=1;    
 else
  if str_a(length('Node')+2)==','
   flag=1;
  end
 end
 %���ϼ�������Ϊȷ���Ƿ�Ϊ�ڵ�����һ��

 if flag==1
  for b=(start_nodele(i)+1):(end_nodele(i)-1) 
 %bΪ�ڵ���������ֹ��
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
%����ڵ������ڵ���������

%Order=floor(log10(NUME));
%offset=10*10^Order;
%LEnew=zeros(NUME,5);
%end    

for i=1:length(start_elements)
 end_elements(i)=Lin_orders(find(Lin_orders==start_elements(i))+1);
 %i=1:2
 %��Lin_orders������start_elements����Ԫ��ֵ��ͬʱ��Lin_orders��ͬԪ����һԪ�ظ�ֵ�� end_elements
 q=1;
 str_b=FFLINES{start_elements(i)}(~isspace(FFLINES{start_elements(i)}));
 %�����ո�
 flag=0;
 if length(str_b)==(length('Element')+1)
  flag=1;    
 else
  if str_b(length('Element')+2)==','
   flag=1;
  end
 end
 %���ϼ�������Ϊȷ���Ƿ�Ϊ��Ԫ�����һ��

 if flag==1
  for b=(start_elements(i)+1):(end_elements(i)-1) 
 %bΪ��Ԫ�����ֹ��
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


% ������������ļ�
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