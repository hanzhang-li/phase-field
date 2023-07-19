% 本程序用于设置初始相场值
function [totalph,totalhm]=initialph(epx0,epy0,epx1,epy1,phi0,hm0,snode,nodecoor,xlc,totalph,totalhm)
    vert0(1)=epx1-epx0; vert0(2)=epy1-epy0;
    legab=norm(vert0);
    % 对每一个单元节点进行循环
    for nd=1:snode
        vert1(1)=nodecoor(nd,1)-epx0; vert1(2)=nodecoor(nd,2)-epy0;
        % 依据向量进行判断点是否在线段内部
        rd=dot(vert0,vert1)/legab;
        if rd>=0 && rd<=legab
            % 求解点到直线的距离
            legd=abs(vert0(1)*vert1(2)-vert0(2)*vert1(1))/legab;
            if legd<=xlc
                totalph(nd)=phi0*exp(-legd/xlc);
                totalhm(nd)=hm0*exp(-2*legd/xlc);
            end
        end
    end
end