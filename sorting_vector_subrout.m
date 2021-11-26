function [neword,neweval,newevect] = sorting_vector_subrout(size,oldeval,oldord,oldevect)
for i=1:size %copying into new
    neword(i)=oldord(i);
    neweval(i)=oldeval(i);
end
newevect=oldevect;
tempv=zeros(1,size);
for i=1:size %sorting the new
    for j=i:size
        %temp1=0;temp2=0; da cancellare
        if neweval(j)>neweval(i)
            temp1=neweval(j);
            neweval(j)=neweval(i);
            neweval(i)=temp1;
            
            temp2=neword(j);
            neword(j)=neword(i);
            neword(i)=temp2;    
            
            tempv(1,:)=newevect(:,i);
            newevect(:,i)=newevect(:,j);
            newevect(:,j)=tempv(1,:);
        end
    end
end
