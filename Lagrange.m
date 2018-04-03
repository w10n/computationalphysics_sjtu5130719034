function y=Lagrange(x,pointx,pointy)
   n = size(pointx,2);
   y = 0;
   for i=1:n
       p = 1;
       for j=1:n
           if (i~=j)
               p=p.*(x-pointx(j))/(pointx(i)-pointx(j));
           end
       end
       y = y+p*pointy(i);
   end
   y = simplify(y);
end