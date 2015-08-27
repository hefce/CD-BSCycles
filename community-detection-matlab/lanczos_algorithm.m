%dvenugopalarao%
% Date:09-02-2012,By D.Seshu Kumar.

% ----------- Lanczos ALgorithm  for Symmetric and Unsymmetric RK ----

function [PI,V,D]=lanczos_algorithm(RK)
%% RK is correlation matrix 
% PI is tridiagonal trnasformation matrix
% V is eigenvectors of PI
% D is diagonal matrix containing eigenvalues of PI
% The following books are referred for Lanczos Algorithm and Unsymmetric
% Lanczos Algorithm: B.N. Parlett, The Symmetric Eigenvalue Problem,
% Prentice Hall, Englewood Cli?s, NJ, 1980. Golub, G.H., Van Loan, C.F.
% Matrix Computations, 3rd Edition, The Johns Hopkins University Press,1996.
% Authors: 
% Seshu Kumar Damarla, Ph.D Student,Department of Chemical Engineering, National Institute of Technology Rourkela, India
% Prof.Madhusree Kundu, Associate Professor,Department of Chemical Engineering, National Institute of Technology Rourkela, India
%%
r0=(ones(size((RK),1),1))/sqrt(length(RK));
beta0=norm(r0);
beta=zeros(length(RK),1);
r=zeros(length(RK),length(RK));
alpha=zeros(length(RK),1);
% PI=zeros(length(RK),length(RK));
q=zeros(length(RK),length(RK));         % for storage of Lanczos vectors
% U=zeros(length(RK),1);
% --------------- Initial data for Lanczos Unsymmetric method -------------

p1=(ones(size((RK),1),1))/sqrt(length(RK));
P=zeros(length(RK),length(RK));
Q1=(ones(size((RK),1),1))/sqrt(length(RK));
Q=zeros(length(RK),length(RK));
S=zeros(length(RK),length(RK));
p0=0;
Q0=0;
t=0;
gama=zeros(length(RK),1);
% -------------------- ----------------------------------------------------
if (RK==RK') 
    disp('-----------------------------')
    disp(':correlation matrix is symmetric')
    disp('-----------------------------')
    c=1;
else
    disp('-----------------------------')
    disp(':correlation matrix is unsymmetric ')
    disp('-----------------------------')
    c=0;
end
b1=1;
q0=0;
% PI(1,:)=[];
% PI(:,1)=[];
% while(tracepi/traceRK<=e)
if c==1
    disp('-----------------------------')
    disp(':Symmetric Lanczos Algorithm')
    disp('-----------------------------')
    % ----------------Symmetric Lanczos Algorithm -------------------
    for k=1:1:length(RK)                
        if k==1                             % Calculations of symmetric matrix begin
            q(:,k)=r0/beta0;
        else
            q(:,k)=r(:,k-1)/beta(k-1);
        end
            U=RK*q(:,k);
            if k==1
                r(:,k)=U-q0*beta0;
            else
                r(:,k)=U-q(:,k-1)*beta(k-1);
            end
            alpha(k)=q(:,k)'*r(:,k);
            r(:,k)=r(:,k)-alpha(k)*q(:,k);
            sume=0;
    for k1=k:1:1                                     % Reorthogonalization
        sum1=sume+q(:,k1)*(q(:,k1)'*r(:,k));
        sume=sum1;
    end
    r(:,k)=r(:,k)-sum1;
    beta(k)=normest(r(:,k));                           % 2-norm 
%     k=k+1;
    end
% ---------------- End of calculations of symmetric matrix -----------

% ------------------ Unsymmetric Algorithm ------------------------
elseif c==0 
    disp('-----------------------------')
    disp(':Unsymmetric Lanczos Algorithm')
    disp('-----------------------------')
    s0=p1;r0=Q1;
    while (t<=length(RK))      % Unsymmetric Algorithm : (sum(r0)~=0)||(sum(s0)~=0)||(s0'*r0~=0)||
        if t==0
            beta0=normest(r0);
            gama0=(s0'*r0)/beta0;
            Q(:,t+1)=r0/beta0;
            P(:,t+1)=s0/gama0;
        elseif t>0
            beta(t)=normest(r(t));
            gama(t)=(S(:,t)'*r(:,t))/beta(t);
            Q(:,t+1)=r(:,t)/beta(t);
            P(:,t+1)=S(:,t)/gama(t);
        end
        t=t+1;
        alpha(t)=P(:,t)'*RK*Q(:,t);
        if t==1
        r(:,t)=(RK-alpha(t)*eye(length(RK)))*Q(:,t)-gama0*Q0;
        S(:,t)=(RK-alpha(t)*eye(length(RK)))'*P(:,t)-beta0*p0;
        elseif t>1
            r(:,t)=(RK-alpha(t)*eye(length(RK)))*Q(:,t)-gama(t-1)*Q(:,t-1);
            S(:,t)=(RK-alpha(t)*eye(length(RK)))'*P(:,t)-beta(t-1)*P(:,t-1);
        end
        r0=r(:,t);
        s0=S(:,t);
        if sum(r0)==0 || sum(s0)==0 || sum(r0)==0 && sum(s0)==0
            disp('---------------------------------')
            disp('1:either r0 or so or both is zero')
            disp(':breakdown occured')
            disp('---------------------------------')
            c=2;
            break;
        elseif sum(r0)~=0 && sum(s0)~=0 && (s0'*r0==0)
            disp('-------------------------------------')
            disp(':breakdown occured')
            c=2;
            disp(':there is no information about eigen values of tridiagonal matrix')
            disp(':use Lanczos look ahead algorithm to cure breakdown')
            disp('-----------------------------------------------------')
            break;
        end
    end
    if s0'*r0==0 || sum(r0)==0 || sum(s0)==0
%         cl=2;
        disp('-----------------------------')
        disp(':A Look ahead Lanczos Algorithm')
        disp('---- LAL ---------------------')
         % Initialization
%         Pl=zeros(length(RK),length(RK));
%         Ql=zeros(length(RK),length(RK));
%         Beta=zeros(length(RK),1);
        Gama=zeros(length(RK),1);
        B=zeros(length(RK),1);
        Alpha=zeros(length(RK),1);
%         u=1;
        z=zeros(length(RK),1);
        z(1)=1;
%         L=1;
        L=1+size((p0),1);
        R=zeros(length(RK),length(RK)+1);
        Sl=zeros(length(RK),length(RK)+1)';
        R(:,L)=(ones(length(RK),1)); %/sqrt(length(RK));
        Sl(L,:)=(ones(length(RK),1))'; %/sqrt(length(RK));
        ql=zeros(length(RK),length(RK));
        pl=zeros(length(RK),length(RK));
        Rn=zeros(length(RK),1);
        Sn=zeros(length(RK),1);
        w=Sl(L,:)*R(:,L);
        tol=0.0001;
% ------------------------ Look Ahead Lanczos Algorithm -----------------
        for I=1:1:length(RK)
                       % s': variable s is not designated as s'
            if I==1                      % but the values of s are transposed(row.and col.are interchanged)
            R(:,L+1)=RK*R(:,L)-q0*z(I)*w; % This mehtod is maintained throughout LAL algorithm  
            Sl(L+1,:)=(Sl(L,:))*RK-w*p0';               % without transpose
            else                                        % this consideration also holds for P
                R(:,L+1)=RK*R(:,L)-ql(I-1)*z(I)*w;
                Sl(L+1,:)=Sl(L,:)*RK-w*pl(I-1);
            end
            theta=Sl(L,:)*R(:,L+1);                % needed inner products
            wb=Sl(L+1,:)*R(:,L+1);
%             Rn(L)=norm()
            if L==2
            pi1=w/((norm(R(:,L)))*norm(Sl(L,:)));
            else
                pi1=w/(Rn(I-1)*Sn(I-1));
            end
            pi2=0;
            if theta==0
                if abs(pi1)<tol && pi2<tol      % test for failure
                    disp('-----------------------------')
                    dip(':both angles are less than tol')
                    disp(':algorithm will terminate without having any information about Lanczos vectors')
                    disp('------------------------------------------------------------------------------')
                    break;
                end
            else
                tau1=w/theta;
                tau2=wb/theta;
                if L==2
                rnb1=sqrt(((norm(R(:,L+1)))^2)-(2*tau2*(R(:,L)'*R(:,L+1)))+(tau2^2*norm(R(:,L))^2));
                snb=sqrt((norm(Sl(L,:)')^2)-2*tau1*(Sl(L,:)*Sl(L+1,:)')+(tau1^2)*norm(Sl(L+1,:)')^2);
                psi1=theta/(norm(R(:,L))*norm(Sl(L+1,:)));
                
                else
                    rnb1=sqrt(((norm(R(:,L+1)))^2)-(2*tau2*(R(:,L)'*R(:,L+1)))+(tau2^2*Rn(I-1)^2));
                    snb=sqrt((Sn(I-1)^2)-2*tau1*(Sl(L,:)*Sl(L+1,:)')+(tau1^2)*norm(Sl(L+1,:)')^2);
                    psi1=theta/(Rn(I-1)*norm(Sl(L+1,:)));
                end
                psi2=(w*tau2-theta)/(rnb1*snb);
                pi2=min(abs(psi1),abs(psi2));
            end
            biasf=0;         % single step is assumed
            if (abs(pi1))>=((biasf)*pi2);
                disp('-----------------------------------------------------------------')
                disp(':LAL reduces to standard Lanczos algorithm from this point onwards')
                disp('------- single step ----------------------------------------------')
                if L==2
                beta1=norm(R(:,L))*sqrt(pi1);
                else
                    beta1=Rn(I-1)*sqrt(pi1);
                end
                G1=w/beta1;
                ql(:,I)=R(:,L)/beta1;
                pl(I,:)=(Sl(L,:))/G1;
                Gama(I)=z(I)*G1;
                B(I)=beta1;
                R(:,L+1)=R(:,L+1)/beta1;                    % form new residuals
                Sl(L+1,:)=Sl(L+1,:)'/G1;
                Alpha(I)=1/tau1;
                R(:,L+1)=R(:,L+1)-ql(:,I)*Alpha(I);     % biorthogonalization
                Sl(L+1,:)=Sl(L+1,:)-Alpha(I)*pl(I,:);
                if L==2
                Rn(I)=sqrt((norm(R(:,L+1))^2)-(2*Alpha(I)*R(:,L)'*R(:,L+1))+(Alpha(I)^2)*norm(R(:,L))^2)/beta1;
                Sn(I)=sqrt(norm(Sl(L+1,:)')^2-2*Alpha(I)*Sl(L,:)*Sl(L+1,:)'+(Alpha(I)^2)*norm(Sl(L,:)')^2)/G1;
                else
                    Rn(I)=sqrt((norm(R(:,L+1))^2)-(2*Alpha(I)*R(:,L)'*R(:,L+1))+(Alpha(I)^2)*Rn(I-1)^2)/beta1;
                    Sn(I)=sqrt(norm(Sl(L+1,:)')^2-2*Alpha(I)*Sl(L,:)*Sl(L+1,:)'+(Alpha(I)^2)*Sn(I-1)^2)/G1;
                end
                w=(wb/w)-Alpha(I)^2;
                z(I+1)=[1];
            end                 % end of calculations for single step
            L=L+1;
        end
    end
% ------------------ Unsymmetric Lanczos Algorithm ------------------------             
end
e=0.995;
tracePI=0;
traceRK=1;
f=1;
% PI=zeros(f,f);
while ((tracePI/traceRK)<=e)          % termination, i.e.if eigen values of PI  
    PI=zeros(f,f);
    for w=1:1:f                     % are approximately equal to eig(RK).
        for l=1:1:f
            if w==l
                if c==2
                    PI(w,l)=Alpha(w);
                else
                PI(w,l)=alpha(w);
                end
            elseif w==(l+1) && l==(w-1)
                if c==2 
                    PI(w,l)=B(l);
                else
                    PI(w,l)=beta(l);
                end
            elseif l==w+1 && w==l-1 
                if c==1
                PI(w,l)=beta(w);
                elseif c==0 
                    PI(w,l)=gama(w);
                elseif c==2 
                    PI(w,l)=Gama(w);
                end
            end
        end
    end
    [vrk,drk]=eig(RK);
    [vpi,dpi]=eig(PI);
    tracePI=sum(diag(dpi));                  % sum of eigen values
    traceRK=sum(diag(drk));
    f=f+1;
    if f==(length(RK)+1)
        disp('-----------------------------')
        disp(':PI is not approximation of RK')
        disp(':no.of iterations are equal to dimensions of RK')
        disp('-----------------------------------------------')
        break;
    end
    
end
[V,D]=eig(PI);

