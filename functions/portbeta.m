function beta=portbeta(portReturn,maketReturn)
%code by ariszheng@gmail.com 2012-5-17
%Э����������
temp_cov=cov(portReturn,maketReturn);
%������г���Э����/�г��ķ���
beta=temp_cov(1,2)/temp_cov(2,2);