function [nndex,mpiv]=matndx(sdegfree,snode,ndofn,issor,nband,mband)
% 建立二维矩阵到一维矩阵的索引
mpiv=zeros(sdegfree,1);
nkk=0;
for ip1=1:snode
    for j1=1:ndofn
        itotv=issor(ip1,j1);
        if itotv~=0
            for iband=1:nband(ip1)
                ip2=mband(ip1,iband);
                for j2=1:ndofn
                    ktotv=issor(ip2,j2);
                    if itotv>=ktotv&&ktotv~=0
                        nkk=nkk+1;
                        nndex(nkk)=ktotv;
                    end
                end
            end
            mpiv(itotv)=nkk;
        end
    end
end