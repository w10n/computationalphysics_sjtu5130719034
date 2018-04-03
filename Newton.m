function yi=Newton(x,y,xi)
%Newton��ֵ����������һϵ�в�ֵ�ĵ�(x,y)���õ���x=xi���ģ�ţ�ٲ�ֵ�����ֵyi
n=length(x);
m=length(y);
if m~=n
 error('x,y�ĳ��Ȳ�һ��,����������!');
 return
end
A=zeros(n); %������̱�
A(:,1)=y; %���̱��һ��Ϊy
for j=2:n %jΪ�б�
    for i=1:(n-j+1)%iΪ�б�
    A(i,j)=(A(i+1,j-1)-A(i,j-1))/(x(i+j-1)-x(i)); %������̱�
    end
end
%���ݲ��̱�,���Ӧ��ţ�ٲ�ֵ����ʽ��x=xi����ֵyi
N(1)=A(1,1);
for j=2:n
    T=1;
    for i=1:j-1
        T=T*(xi-x(i));
    end
N(j)=A(1,j)*T;
end
yi=sum(N); %��x=xi����ţ�ٲ�ֵ����ʽ���õ���yi��ֵ
end 