% ����������ʩ��λ�ƺ��صĽڵ�͵�Ԫ����ϸ��Ϣ
function [disprelat,dispmpiv]=dispinforfp(fpelem,nnode,ndofn,disprelat,dfsoor,dispmpiv)
% disprelat:ʩ�ؽڵ㣬���ӽڵ㣬�ֲ������Ӧ��ʩfpelem�ؽڵ�����ӽڵ㣬��Ӧ�ľֲ��������ɶ�
a1=length(fpelem);
dispnum=0;
num=0;
for i1=1:a1
    num=num+1;
    for j1=1:nnode
        for p1=1:ndofn
            itotv=dfsoor(fpelem(i1,j1),p1);
            a1=2*(fpelem(i1,j1)-1)+p1;
            b1=2*(j1-1)+p1;
            if itotv==0
                for j2=1:nnode
                    for p2=1:ndofn
                        ktotv=dfsoor(fpelem(i1,j2),p2);
                        a2=2*(fpelem(i1,j2)-1)+p2;
                        b2=2*(j2-1)+p2;
                        if ktotv~=0
                            dispnum=dispnum+1;
                            disprelat=[disprelat;num a1 a2 b1 b2];
                        end
                    end
                end
            end
        end
    end
    dispmpiv=[dispmpiv; dispnum];
end