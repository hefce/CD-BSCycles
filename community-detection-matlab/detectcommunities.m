
%Detecting Communities based on spectral clustering plus a similarity
% matrix captured from cycles

Graph=load('network.dat');

%disp('number of nodes' num2str(size(Graph(1))));
disp('loading :)');

size(Graph)
Adjmat=zeros(max(Graph(:,2)));

for i=1:size(Graph,1)
    Adjmat(Graph(i,1),Graph(i,2))=1;
end;


% adj=CIJctx;
%     for i=1:length(adj)
%         for j=i+1:length(adj)
%             if adj(i,j) == 1
%                 adj(i,j)=0;
%             end
%             if adj(i,j) >0
%                 adj(i,j)=1;
%             end
%            
%         
%         end%     end
% Spmat=dijk(adj');




disp('..');
% Adjmat(find(Adjmat==0))=2;
% Adjmat(find(Adjmat==1))=0;
% Adjmat(find(Adjmat==2))=1;

Spmat=dijk(Adjmat');
Cycles=Spmat+Spmat';
%   Cycles=Cycles+eye(size(Cycles,1));
%Cycles(find(Cycles==Inf))=;
% kernelcb=1./Cycles;
kernelcb=exp(-Cycles);
%kernelcb=kernelcb-eye(size(Cycles,1));
disp('shortest cycles :)');
disp('..');
% Simmat = Simmat-eye(size(Simmat,1));
Communities=load('community.dat');
%numcommunity=2;
numcommunity=max(Communities(:,2))
% numcommunity=max(Communities(:))

[assignment]=cluster_spectral_general(kernelcb,numcommunity,'ang_gen','kmeans');

NMI=nmi(Communities(:,2),assignment')

c1=Communities(:,2);
c2=assignment';

Mem1=c1;
Mem2=c2;

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end
C=Cont;

n=sum(sum(C));
nis=sum(sum(C,2).^2);		%sum of squares of sums of rows
njs=sum(sum(C,1).^2);		%sum of squares of sums of columns

t1=nchoosek(n,2);		%total number of pairs of entities
t2=sum(sum(C.^2));	%sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;		%no. agreements
D=  -t2+t3;		%no. disagreements

if t1==nc
   AR=0;			%avoid division by zero; if k=1, define Rand = 0
else
   AR=(A-nc)/(t1-nc);		%adjusted Rand - Hubert & Arabie 1985
end

RI=A/t1;			%Rand 1971		%Probability of agreement
MI=D/t1;			%Mirkin 1970	%p(disagreement)
HI=(A-D)/t1;	%Hubert 1977	%p(agree)-p(disagree)

disp('detection :)  NMI[Adjusted R, Uadjusted R, Mirkin, Hubert]= ');
[AR,RI,MI,HI] 
% 
% 
% n=length(c1); % number cass  
% ng1=max(c1);
% ng2=max(c2);
% dn=n*(n-1)/2; %number of possible pairwise comparisons of cases(a+b+c+d)
% 
% y1=pdist(c1(:),'ham');
% y2=pdist(c2(:),'ham');
% %size(y1)
% %size(y2)
% ad=sum((y1'==y2')); % number of pairwise concordances (matches (a) and mismatches(d))(a+d)
% bc=sum((y1'~=y2')); % number of pairwise discordances(b+c)
% %Rand Index
% ri=ad/dn;
% 
% a=sum((y1'==0).*(y2'==0));
% b=sum(y1'==0)-a;
% c=sum(y2'==0)-a;
% 
% d=ad-a;
% 
% %Jaccard Index
% jac=a/(dn-d)
% 
% 
% % ffff=assignment;
% % 
% % 
% % nnnn=length(ffff);
% % 
% % 
% % kkkk=max(ffff);
% % Z=zeros(nnnn,kkkk);
% % for i=1:nnnn
% %     Z(i,ffff(i))=1;
% % end;
% % 
% % 
% % 
% % D=diag(sum(kernelcb'));
% % L=D^(-1/2)*kernelcb*D^(-1/2);
% % reo=Z'*L*Z;
% % kk1=trace(reo);
% % kk2=sum(sum(reo));
% % disp('Qcb');Qcb=kk1/kk2
% 
% 
% 
%  mycolors = [
%              0 255 255
%              
%                     ];
%                 mycolors=mycolors./256;
% 
% set(gca, 'ColorOrder', mycolors);
% hold all
% 
% h = zeros(1,3);
% 
% 
% x=[4,6,8,10,12,14];
% %h(1) = plot(x,[0.41,0.46,0.50,0.51,0.52] ,'x-','LineWidth',4);
% 
% 
% 
% h(2) = plot(x,[0.2,0.7,0.95,1,1,1], 'x--','LineWidth',4);
% 
% 
% 
% %h(3) = plot(x,[1,1,1,1,1], 'x--','LineWidth',4);
% 
% 
% 
% legend({'Our Method'},'Location','NorthWest');
% title('Average Degree of Nodes-LFR Directed Networks');
% ylabel('Jaccard Index: J');
% xlabel('Average Degree( K-avg )');
% 
% A=Adjmat;
% 
% %MODULARITY_DIR     Optimal community structure and modularity
% %
% %   Ci = modularity_dir(W);
% %   [Ci Q] = modularity_dir(W);
% %
% %   The optimal community structure is a subdivision of the network into
% %   nonoverlapping groups of nodes in a way that maximizes the number of
% %   within-group edges, and minimizes the number of between-group edges. 
% %   The modularity is a statistic that quantifies the degree to which the
% %   network may be subdivided into such clearly delineated groups. 
% %
% %   Input:      W,      directed (weighted or binary) connection matrix.
% %
% %   Outputs:    Ci,     optimal community structure
% %               Q,      maximized modularity
% %
% %   Note: Ci and Q may vary from run to run, due to heuristics in the 
% %   algorithm. Consequently, it may be worth to compare multiple runs.
% %   Also see Good et al. (2010) Phys. Rev. E 81:046106.
% %
% %   Reference: Leicht and Newman (2008) Phys Rev Lett 100:118703.
% %
% %
% %   2008-2010
% %   Mika Rubinov, UNSW
% %   Jonathan Power, WUSTL
% %   Dani Bassett, UCSB
% 
% 
% %   Modification History:
% %   Jul 2008: Original (Mika Rubinov)
% %   Oct 2008: Positive eigenvalues are now insufficient for division (Jonathan Power)
% %   Dec 2008: Fine-tuning is now consistent with Newman's description (Jonathan Power)
% %   Dec 2008: Fine-tuning is now vectorized (Mika Rubinov)
% %   Sep 2010: Node identities are now permuted (Dani Bassett)
% 
% N=length(A);                            %number of vertices
% n_perm = randperm(N);                   %DB: randomly permute order of nodes
% A = A(n_perm,n_perm);                   %DB: use permuted matrix for subsequent analysis
% Ki=sum(A,1);                            %in-degree
% Ko=sum(A,2);                            %out-degree
% m=sum(Ki);                              %number of edges
% b=A-(Ko*Ki).'/m;
% B=b+b.';                                %directed modularity matrix
% Ci=ones(N,1);                           %community indices
% cn=1;                                   %number of communities
% U=[1 0];                                %array of unexamined communites
% 
% ind=1:N;
% Bg=B;
% Ng=N;
% 
% while U(1)                              %examine community U(1)
%     [V D]=eig(Bg);
%     [d1 i1]=max(diag(D));               %most positive eigenvalue of Bg
%     v1=V(:,i1);                         %corresponding eigenvector
% 
%     S=ones(Ng,1);
%     S(v1<0)=-1;
%     q=S.'*Bg*S;                         %contribution to modularity
% 
%     if q>1e-10                          %contribution positive: U(1) is divisible
%         qmax=q;                         %maximal contribution to modularity
%         Bg(logical(eye(Ng)))=0;         %Bg is modified, to enable fine-tuning
%         indg=ones(Ng,1);                %array of unmoved indices
%         Sit=S;
%         while any(indg);                %iterative fine-tuning
%             Qit=qmax-4*Sit.*(Bg*Sit);   %this line is equivalent to:
%             qmax=max(Qit.*indg);        %for i=1:Ng
%             imax=(Qit==qmax);           %       Sit(i)=-Sit(i);
%             Sit(imax)=-Sit(imax);       %       Qit(i)=Sit.'*Bg*Sit;
%             indg(imax)=nan;             %       Sit(i)=-Sit(i);
%             if qmax>q;                  %end
%                 q=qmax;
%                 S=Sit;
%             end
%         end
% 
%         if abs(sum(S))==Ng              %unsuccessful splitting of U(1)
%             U(1)=[];
%         else
%             cn=cn+1;
%             Ci(ind(S==1))=U(1);         %split old U(1) into new U(1) and into cn
%             Ci(ind(S==-1))=cn;
%             U=[cn U];
%         end
%     else                                %contribution nonpositive: U(1) is indivisible
%         U(1)=[];
%     end
% 
%     ind=find(Ci==U(1));                 %indices of unexamined community U(1)
%     bg=B(ind,ind);
%     Bg=bg-diag(sum(bg));                %modularity matrix for U(1)
%     Ng=length(ind);                     %number of vertices in U(1)
% end
% 
% s=Ci(:,ones(1,N));                      %compute modularity
% Q=~(s-s.').*B/(2*m);
% Qnewmanleitch=sum(Q(:))
% 
% Ci_corrected = zeros(N,1);              % DB: initialize Ci_corrected
% Ci_corrected(n_perm) = Ci;              % DB: return order of nodes to the order used at the input stage.
% Ci = Ci_corrected;
% lengthclusters=max(Ci)% DB: output corrected community assignments
% 


















