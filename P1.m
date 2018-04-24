close;
clear;
clc;
%%q(1)
U = {'1','2','3','4'};
No = [1, 2, 3, 4];
i = [2, 3, 2, 4, 3];
j = [1, 1, 3, 3, 4];
n = 4;
p=0.85;
delta = 3.75*10^(-11);
G = sparse(i,j,1,n,n);
full(G);
display(full(G));
c = full(sum(G));
c1 = full(sum(G'));
display(c);
%display(c1);
k = find(c~=0);
D = sparse(k,k,1./c(k),n,n);
%display(full(D))
e = ones(n,1);
I = speye(n,n);
z = ((1-p)*(c~=0) + (c==0))/n;
%display(z);
%display(e);
A = p*G*D + e*z;
display(A);
x = (I - p*G*D)\e;
%x = (I - A)\e;
x = x/sum(x);
display(x);
% construct two matrixes to present the result
%outputhead=['      url','      page-rank ','       in','         out'];
%PR = [No', x, c1', c'];
%display(PR);
%PR1 = sortrows(PR,-2);
%pagerank(U,G);
%for i=1:n
% output(i,:)=PR1(i,:);
%end
%format long
%format shortE
%disp(outputhead)
%disp(output)
PR1 =  P_R(No', x, c1', c',n);
bar(No,x);
title('Page-Rank');
%%q(2)
torr = 10^(-13);
%Delta = 1;
x1 = ones(n,1)/n;
x0 = zeros(n,1)/n;
sigma = 10^(-11);
%while  norm(x1-x0)>=sigma
while   sqrt(sum((x1-x0).^2))>=sigma
    x0 = x1;
    x1 = A*x0;
    %x1 =  G*x0 + e*(z*x0);
    %display(x2); 
end
%display(x1);
x1 = x1/sum(x1);
display(x1);
PR2 =  P_R(No', x1, c1', c',n);
%%q(3)
x_n = G_E(I-p*G*D,e);
x_n = x_n/sum(x_n);
display(x_n)
PR2 =  P_R(No', x_n, c1', c',n);
