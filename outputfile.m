% 本程序用于输出计算的文件
function outputfile(strout,snode,selement,nodecoor,element,totalph,totalhm,totaldisp,totalstress,pristress)
    fileid=fopen(strout,'wt');
    fprintf(fileid,'TITLE="PF-FEM to Tecplt 2D"\n');
    fprintf(fileid,'VARIABLES="X","Y","Phase1","Hm","Dispx","Dispy","SX","SY","SXY","S1"\n');
    fprintf(fileid,'ZONE T="FE QUADRILATERAL"\n');
    formatSpec='N= %d  E= %d  ZONETYPE=FEQUADRILATERAL\n';
    fprintf(fileid,formatSpec,snode,selement);
    fprintf(fileid,'DATAPACKING=POINT\n');
    fprintf(fileid,'DT=(SINGLE SINGLE SINGLE SINGLE )\n');

    % 按照节点输入节点坐标以及值的输出
    outputdata=zeros(snode,10);
    outputdata(:,1:2)=nodecoor*1000;               % 节点
    outputdata(:,3)=totalph;                  % 相场值
    outputdata(:,4)=totalhm;                  % 应变能
    outputdata(:,5)=totaldisp(:,1);           % x方向位移
    outputdata(:,6)=totaldisp(:,2);           % y方向位移
    outputdata(:,7)=totalstress(:,1);         % x方向应力
    outputdata(:,8)=totalstress(:,2);         % y方向应力
    outputdata(:,9)=totalstress(:,3);         % xy方向应力
    outputdata(:,10)=pristress(:,1);               % 主应力
    
    formatSpec='%12.6f  %12.6f  %10.6f  %20.5f  %10.7f  %10.7f  %15.5f  %15.5f  %15.5f  %15.5f\n';
    fprintf(fileid,formatSpec,outputdata');
    % 输出单元信息
    
    formatSpec='%6d  %6d  %6d  %6d \n';
    fprintf(fileid,formatSpec,element');
    fclose all;
end