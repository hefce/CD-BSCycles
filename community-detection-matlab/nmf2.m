% In the name of God
function nmf2(kernelcb)

clc
n=128;
m=4;
A=zeros(n);
f=[0,10,20];
p=0.7;
q=0.1;
%for k=2:length(f);
 %   for i=f(k-1)+1:f(k)
  %    for j=i+1:f(k)
   %        x=rand(1);
    %       if x<p
     %         A(i,j)=1;
      %        A(j,i)=1;
       %    end
      % end
   % end
%end
 %for k=2:length(f)
  %  for i=f(k)+1:n
   %     for j=f(k-1)+1:f(k)
    %        y=rand(1);
     %       if y<q          
      %        A(i,j)=1;
       %       A(j,i)=1;
		%	end 
       % end 
   % end 
 %end 
link=load('net.txt');
for i=1:length(link);
    A(link(i,1),link(i,2))=1;
    A(link(i,2),link(i,1))=1;
end
   %A=[0 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0;
    %1 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0;
    %1 1 0 1 1 1 0 0 0 0 0 0 0 0 0 0;
    %0 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0;
    %1 1 1 1 0 1 0 0 0 0 0 0 0 0 0 0;
    %1 0 1 1 1 0 1 0 0 0 0 0 0 0 0 0;
    %0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0;
    %0 0 0 0 0 0 1 0 1 1 0 1 0 0 0 0;
    %0 0 0 0 0 0 0 1 0 1 1 1 0 0 0 0;
    %0 0 0 0 0 0 0 1 1 0 1 1 0 0 0 0;
    %0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 0;
    %0 0 0 0 0 0 0 1 1 1 1 0 1 0 1 1;
    %0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 1;
    %0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1;
    %0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1;
    %0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0];
%A=[0 1 1 0 0;
 %  1 0 1 0 0;
 %  1 1 0 1 1;
 %  0 0 1 0 1;
 %  0 0 1 1 0];
%deg=sum(A);
D=diag(sum(A));
L=A-D;
%for i=1:n
%   for j=i+1:n
%       if(L(i,j)==1)
%         x1=rand(1);
%         if x1<0.5
%          L(i,j)=0;
%         else
%          L(j,i)=0;
%         end
%       end
%   end
%end
%v=abs(A-(deg'*deg./(sum(deg))));
%for i=1:n
%for j=1:n
%L(i,j)=L(i,j)/(sqrt(L(i,i).*L(j,j)));
% end
%end
%v =abs(corrcoef(L));
v=kernelcb;
[hhh,ggg]=eig(v);
for i=1:n
    if i<10
        s(i)=-1;
        else
        s(i)=1;
     end
end
figure(3)
for i=1:n
    e(i)=ggg(i,i);
    plot(real(e),imag(e),'b.')
    hold on
end
DDD=sum(sum(v(:,:)));
M=s*v*s';
M/DDD
%v=abs((L));
%v=L.*L;
W=abs(randn(n,m));
H=abs(randn(m,n));
for t=1:200
wh=W*H;
H1=W'*v;
H2=W'*W*H;
H=H.*H1./H2;
W1=v*H';
W2=W*H*H';
W=W.*W1./W2;
z=(v-wh).^2;
E=sum(sum(z));
%/(sum(sum(v.^2)))*(sum(sum(wh.^2)));
end
for i=1:m
w(:,i)=W(:,i)./(sqrt(sum(W(:,i).^2)));
end
w;
x=zeros(n,2);
for i=1:n
    x(i,1)=i;
    x(i,2)=find(w(i,:)==max(w(i,:)));
end
figure(1)
plot(x(:,1),x(:,2),'r.');
figure(2)
 N=60:n;
plot(N,w(60:n,1),'r.-')
hold on
plot(N,w(60:n,2),'b.-')
plot(N,w(60:n,3),'k.-')
plot(N,w(60:n,4),'g.-')
%plot(N,w(:,4),'g.-')
for i=1:n
    for j=1:n
        if x(i,2)<x(j,2)
           y=x(i,2);
           x(i,2)=x(j,2);
           x(j,2)=y;
           y1=x(i,1);
           x(i,1)=x(j,1);
           x(j,1)=y1;
           
        end
    end
    x(i,3)=i;
end
x
rrr=0;
 for i=1:n
     for j=1:i
         if A(i,j)==1
         rrr=rrr+1;
         link1(rrr,1)=i;
         link1(rrr,2)=j;
         end
     end
 end
fid=fopen('net.txt','wt');
for i=1:length(link1)
  fprintf(fid,'%d %d \n',link1(i,1),link1(i,2));
end
fclose(fid);
siz=zeros(1,3);
for i=1:n
    for j=1:m
        if x(i,2)==j
           siz(j)=siz(j)+1;
        end
    end
end
M=0;
for l=1:m
   for i=1+(i-1)*siz,siz+(i-1)*(n-siz)
      for j=1+(i-1)*siz,siz+(i-1)*(n-siz)
         M=M+v(x(i,1),x(j,1))
      end 
   end 
end 
end



