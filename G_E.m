function x = G_E(A,b)
ab = [A,b];
[r,c] = size(ab);
for k=1:r-1
    [rm,im]=max(abs(ab(k:r,k)));
    im = im +k-1;
    if(ab(im,k)~=0)
        if(im~=k)
            ab([k im],:)=ab([im k],:);
        end
    end
    for i=k+1:r
        ab(i,k:c)=ab(i,k:c)-ab(i,k)/ab(k,k)*ab(k,k:c);
    end
end
x=zeros(r,1);
x(r)=ab(r,c)/ab(r,r);
for i=r-1:-1:1
    x(i)=(ab(i,c)-ab(i,i+1:r)*x(i+1:r))/ab(i,i);
end