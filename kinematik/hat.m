function vh=hat(v)
if (length(v)==6)
   vh=zeros(4,4);
   vh(1,2)=-v(6);
   vh(1,3)=v(5);
   vh(2,1)=v(6);
   vh(2,3)=-v(4);
   vh(3,1)=-v(5);
   vh(3,2)=v(4);
   vh(1:3,4)=v(1:3);
else
   vh=zeros(3,3);
   vh(1,2)=-v(3);
   vh(1,3)=v(2);
   vh(2,1)=v(3);
   vh(2,3)=-v(1);
   vh(3,1)=-v(2);
   vh(3,2)=v(1);
end

