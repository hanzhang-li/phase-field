%该程序用于在paraview中绘制网格，二维三维皆可。
fid = fopen('danzhoulashen75.vtk', 'w'); 
% VTK files contain five major parts
% 1. VTK DataFile Version
fprintf(fid, '# vtk DataFile Version 2.0\n');
% 2. Title
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid, ['POINTS ' num2str(length(nodecoor)) ' float']);
if size(nodecoor,2)==3
    formatSpec='\n%18.12f   %18.12f  %18.12f';
    fprintf(fid,formatSpec,nodecoor'); 
    sele=[selement; 9*selement];
    spec = '\nCELLS %d %d';
    fprintf(fid,spec,sele');
    elem=zeros(selement,9);
    elem(:,2:9)=element(:,[7:8,5:6,3:4,1:2]);
    elem(:,2:9)=elem(:,2:9)-1;
    elem(:,1)=8;
    formatSpec='\n%5d  %6d   %6d   %6d   %6d  %6d  %6d  %6d  %6d';
    fprintf(fid,formatSpec,elem');
    fprintf(fid, ['\nCELL_TYPES ' num2str(selement) '\n']);
    kelem(1:selement)=12;
    formatSpec='%5d\n';
    fprintf(fid,formatSpec,kelem);
    fprintf(fid, ['POINT_DATA ' num2str(length(nodecoor))]);
    fprintf(fid, ['\nVECTORS ', 'Displacement',' float\n']);
    formatSpec='%18.12f   %18.12f  %18.12f\n';
    fprintf(fid,formatSpec,totaldispf'); 
    fprintf(fid, ['SCALARS ', 'Phi',' float 1\n']);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    formatSpec='%24.20f\n';
    fprintf(fid,formatSpec,totalphi'); 
else
    nodecoor2=zeros(length(nodecoor),3);
    nodecoor2(:,1:2)=nodecoor(:,1:2);
    nodecoor2(:,3)=0;
    formatSpec='\n%18.12f   %18.12f  %18.12f';
    fprintf(fid,formatSpec,nodecoor2'); 
    sele=[selement; 5*selement];
    spec = '\nCELLS %d %d';
    fprintf(fid,spec,sele');
    elem=zeros(selement,5);
    elem(:,2:5)=element(:,1:4)-1;
    elem(:,1)=4;
    formatSpec='\n%5d  %6d   %6d   %6d   %6d';
    fprintf(fid,formatSpec,elem');
    fprintf(fid, ['\nCELL_TYPES ' num2str(selement) '\n']);
    kelem(1:selement)=9;
    formatSpec='%5d\n';
    fprintf(fid,formatSpec,kelem);
    fprintf(fid, ['POINT_DATA ' num2str(length(nodecoor))]);
    fprintf(fid, ['\nVECTORS ', 'Displacement',' float\n']);
    totaldisp(:,3)=0;
    formatSpec='%18.12f   %18.12f  %18.12f\n';
    fprintf(fid,formatSpec,totaldisp'); 
    fprintf(fid, ['SCALARS ', 'Phi',' float 1\n']);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    formatSpec='%24.20f\n';
    fprintf(fid,formatSpec,totalph(:)'); 
end
fclose(fid);