    % 本程序用于建立力场和相场单元矩阵并装载至整体矩阵
    function [onetk,trhsk,oneph,trhsp]=golabkmat(element,totalstress,totalph,totalhm,tnijmart,tdetmart, ...
        totbmart,tndxmart,xk,dmat,gc,xlc,elmssor,mpiv,mpiv2,nnode,dnum,sdegfree,snode,nndex,nndex2)
        % 初始化劲度和荷载矩阵
        trhsk=zeros(1,sdegfree);       % 力场右端项矩阵
        trhsp=zeros(1,snode);          % 相场右端项矩阵
        onetk=zeros(1,length(nndex));  % 力场刚度矩阵
        oneph=zeros(1,length(nndex2)); % 相场刚度矩阵
        for i=1:size(element,1)
            % 单元所含节点矩阵
            elemnode2=element(i,:);       

            % 提取生成单元劲度矩阵所需的信息  3个总应力 3个总应变 1个phi值 1个hm值 
            estres=totalstress(elemnode2,:);
            ephi=totalph(elemnode2);
            ehm=totalhm(elemnode2,1); 
            tdetjac=tdetmart(nnode*(i-1)+1:nnode*(i-1)+4);
            tbmat=totbmart(12*(i-1)+1:12*(i-1)+12,:);
            tndx=tndxmart(8*(i-1)+1:8*(i-1)+8,:);            
            
            % 计算单元劲度矩阵
            [delemk,rhsk,delemphi,rhsp]=elemart ...
                (tdetjac,tbmat,tndx,ephi,ehm,estres,dmat,tnijmart,gc,xlc,xk,nnode);     


            % 将delemk和trhsk装载到整体矩阵中
            dssor=elmssor(i,:);
            for j1=1:dnum
                if dssor(j1)~=0
                    for j2=1:dnum
                        if dssor(j2)~=0 && dssor(j1)>=dssor(j2)
                            if dssor(j1)==1
                                onetk(1)=onetk(1)+delemk(j1,j2);
                            else
                                for p1=mpiv(dssor(j1)-1)+1:mpiv(dssor(j1))
                                    if dssor(j2)==nndex(p1)
                                        onetk(p1)=onetk(p1)+delemk(j1,j2);
                                    end                             
                                end
                            end
                        end
                    end
                    trhsk(dssor(j1))=trhsk(dssor(j1))+rhsk(j1);
                end  
            end

            % 将delemphi和oneph装载到整体矩阵中 
            pssor=element(i,:);
            for k1=1:nnode
                for k2=1:nnode
                    if pssor(k1)>=pssor(k2)
                        if pssor(k1)==1
                             oneph(1)=oneph(1)+delemphi(k1,k2);
                        else
                             for p2=mpiv2(pssor(k1)-1)+1:mpiv2(pssor(k1))
                                 if pssor(k2)==nndex2(p2)
                                     oneph(p2)=oneph(p2)+delemphi(k1,k2);
                                 end                             
                             end
                        end
                    end
                end 
                trhsp(pssor(k1))=trhsp(pssor(k1))+rhsp(k1);
            end 
        end 
    end
