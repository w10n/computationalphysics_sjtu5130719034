%% P1 

clear all
close all
clc
format long
%%----Question1.(a)----%%
a=5;
g = @(u) 1./2.*(u+a./u);
for u=1:1:4
    [y1(u),l1,x]=FixedPoint(g,u);
    n=1:1:(l1+1);
    figure(1)
    str=['a=5,initial=',num2str(u),',xn iterate with n'];
    subplot(2,2,u),plot(n,x),title(str),xlabel('iteration times n'),ylabel('xn');
end
u=1:1:200;
for n=1:200
    [y2(n),l2]=FixedPoint(g,u(n));
end
figure(2)
plot(u,y2),title('a=5,the final xn changes with the initial x '),xlabel('the initial x'),ylabel('the final xn')
%%----Question1.(b)----%%
a=5;
[y3,l3,x3]=FixedPoint(g,2);
rel_error=abs((y3-sqrt(5))/sqrt(5));
fprintf('The result of interation is %d when a=5 and initial x=2\n',y3);
fprintf('The relative error of y3 compared with sqrt(5) is %d\n',rel_error);
%%----Question1.(c)----%%
%
%% P2
%%----Question1.(a)----%%
clc,clear;
result=zerosofNo2(10,1);

% the result is listed in the matrix "result" 
% first line is the even func results, second line is the odd func results.

fprintf('the bound-state energies for even wave function at V=10 is %5.5f\n',result(1,:))
fprintf('the bound-state energies for even wave function at V=10 is %5.5f\n',result(2,1))

%%----Question1.(b)----%%
% examine ground state energy with various Vs 
% use for loop to examine 60 Vs, where V=2^(t/10)

for t=1:60
    V=2^(t/10);
    result=zerosofNo2(V,1);
    ground1(t)=min(result(find(result~=0)));
    depth(t)=V;
end
figure(1);
plot(depth,ground1);
xlabel('well depth');
ylabel('ground energy');
title('ground energy - well depth @ width=1');

%%----Question1.(c)----%%
% the same as (b),choose the same equation, set a=5


for t=1:60
    a=0.1*2^(t/10);
    result=zerosofNo2(100,a);
    ground2(t)=min(result(find(result~=0)));
    width(t)=a;
end
figure(2);
plot(width,ground2);
xlabel('well width');
ylabel('ground energy');
title('ground energy - well width @ depth=100');
%% P3

clear all
close all
clc
format long
fprintf('Question3:\n')
e2=14.4;V0=1.09*1000;r0=0.330;
g=@(x) -e2./x+V0.*exp(-x./r0);%V(r)
gm=@(x) e2./x-V0.*exp(-x./r0);%-V(r) in order to use fminsearch
f=@(x) e2./x.^2-V0.*exp(-x./r0)./r0; % the derivative of V(r)
df=@(x) -2*e2/x.^3+V0.*exp(-x./r0)./r0^2;%the second derivative of V(r)
x=0:0.01:1;
y1=feval(g,x);
y2=feval(f,x);
subplot(2,1,1),plot(x,y1); xlabel('r'); ylabel('V(r)'); grid on;title('Q3:Please take the point near the max extreme value in the first plot');
subplot(2,1,2),plot(x,y2); xlabel('r'); ylabel('dV(r)/r'); grid on;
Max=30;  % maximum trials
tolerr=1e-7;
n=1;delta=1;xnew=zeros(1,Max);
[xnew(1),ynew(1)]=ginput();
[r1,fval]=fminsearch(gm,xnew(1));
while delta>tolerr && n<Max 
    %fprintf('n=%d, x=%20.10f\n', n,xnew(n));
    n=n+1;
    xnew(n)=xnew(n-1)-feval(f,xnew(n-1))/feval(df,xnew(n-1));
    delta=abs(xnew(n)-xnew(n-1));    
end
if n<Max && feval(gm,xnew(n))<0
    fprintf('Using my method:x=%3.10f, dV(r)/r=%3.10f, V(r)=%3.10f\n', xnew(n),feval(f,xnew(n)),feval(g,xnew(n)));
else
    disp('Cannot find the extreme we want!');
end
fprintf('Using built-in fminsearch:x=%3.10f, dV(r)/r=%3.10f, V(r)=%3.10f\n',r1,feval(f,r1),-fval)
fprintf('The relative error of req is %f\n',abs(r1-xnew(n))/r1)

%% P5

clc
clear
close
format long
%Problema see in the report
%(b)
str1='please choose a ponit near the root';
f=@(x) 5*x^7+2*x-1;
fp=@(x) 35*x^6+2;
fpp=@(x) 210*x^5;
g=@(x) 1/x^3-10;
gp=@(x) -3/x^4;
gpp=@(x) 12/x^5;

figure(1)
ezplot(f,[0,1])
xlabel('x');ylabel('y');
title('f(x)=5*x^7+2*x-1');
text(0.1,3,str1)
xy1=ginput();
x1=xy1(1,1);

figure(2)
ezplot(g,[0.1,1])
xlabel('x');ylabel('y');
title('g(x)=1/x^3-10');
text(0.3,10,str1)
xy2=ginput();
x2=xy2(1,1);

[rh1 rhh1]=Halley(f,fp,fpp,x1,50,10^-10);
[rn1 rnn1]=Newton2(f,fp,x1,10^-10,50);
[rh2 rhh2]=Halley(g,gp,gpp,x2,50,10^-10);
[rn2 rnn2]=Newton2(g,gp,x2,10^-10,50);

fprintf('the root of f(x) by Halley''s method is %20.10f\n',rh1)
fprintf('the root of g(x) by Halley''s method is %20.10f\n',rh2)
z=1:50;
figure(3)
subplot(2,1,1)
plot(z,rhh1,'r',z,rnn1,'b')
xlabel('steps');ylabel('x')
title('a comparison of Newton and Halley method of f(x)')
subplot(2,1,2)
plot(z,rhh2,'r',z,rnn2,'b')
xlabel('steps');ylabel('x')
title('a comparison of Newton and Halley method of g(x)')