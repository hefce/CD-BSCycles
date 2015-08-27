function [Q]=modularity(A, Ci)

N=length(A);                            %number of vertices
n_perm = randperm(N);                   %DB: randomly permute order of nodes
A = A(n_perm,n_perm);                   %DB: use permuted matrix for subsequent analysis
Ki=sum(A,1);                            %in-degree
Ko=sum(A,2);                            %out-degree
m=sum(Ki);                              %number of edges
b=A-(Ko*Ki).'/m;
B=b+b';                                %directed modularity matrix

cn=1;                                   %number of communities
U=[1 0];                                %array of unexamined communites

ind=1:N;
Bg=B;
Ng=N;

 s=zeros(N,max(Ci));
 for i=1:N
     s(i,Ci(i))=1;
 end;




                     %compute modularity
Q=s'*B*s;
Q=trace(Q);
