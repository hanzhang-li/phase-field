% 该段程序用于确定同一节点出现在几个单元之中
function [linknode]=nodelink(snode,element)
linknode=zeros(1,snode);
[m,n]=size(element);
for i=1:m
    for j=1:n
         linknode(element(i,j))=linknode(element(i,j))+1;
    end
end
