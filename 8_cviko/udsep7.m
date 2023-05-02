function [A D shat Shat1 signal signal1 M avSIR SIR A0 ch]=udsep7(x,L,r,numiter,ARorder,Atrue,Strue)
%
% Blind separation for undetermined mixtures of
% piecewise stationary sources
%
% Coded by Petr Tichavsky, May 2010. 
%
% inputs:
% x ..........  input signal d x N, d is the dimension
% L ..........  number of blocks
% r ..........  number of sources (rank), r < d^2/2 & r < L
% numiter ....  number of iterations
% ARorder ....  order of model for max. SDR beamformer
% Atrue ......  true mixing matrix (due to SIR and SDR computation)
% Strue ......  true sources (due to SIR and SDR computation)
% outputs:
% A  .........  estimated mixing matrix d x r
% D  .........  estimated variances of the sources in each block, L x r
% A0..........  initial estimate of the mixing matrix by sobium
% shat .......  reconstructed sources r x N by the maximum SIR beamformer
% Shat1 ......  reconstructed sources r x N by the maximum SDR beamformer
% signal .....  energy decomposition matrix for shat for exact SIR
% computation
% signal1 ....  energy decomposition matrix for Shat1 for exact SIR
% computation
% SIR ........  estimated SIR of each source in shat in each block, L x r

if nargin<4
   numiter=100;
end   
[d N]=size(x);
Nb=floor(N/L);
M=zeros(d,d*L);
M0=zeros(d,d);
for l=1:L
    ll=(l-1)*Nb;
    l2=(l-1)*d;
    M(:,l2+1:l2+d)=x(:,ll+1:ll+Nb)*x(:,ll+1:ll+Nb)'/Nb;%hessian(x(:,ll+1:ll+Nb),zeros(d,1));
    M0=M0+M(:,l2+1:l2+d);
end
[H E]=eig(M0);
W_est=diag(1./sqrt(abs(diag(E)/L)))*H';
for l=1:L   %%%%%%%%%%%%%%%%%%%%   pre-whitening
    ll=(l-1)*d;
    M(:,ll+1:ll+d)=W_est*M(:,ll+1:ll+d)*W_est';
end
%[A A0 iter]=GN2(M,r);  %%% performs the tensor factorization
[A A0 D iter]=LM1sb(M,r,numiter);
ch=chybajina(M,A,D);
%if ch>0.5
%   'chyba je velka'
%end  
A=W_est\A;
A0=W_est\A0;
%A=Atrue;
am=sqrt(sum(A.^2));
A=A*diag(1./am); D=D*diag(am.*am); 
SIR=zeros(L,r);
signal=zeros(r,r);
shat=zeros(r,N);      %%% computing the optimum beamformer and SIR in each
for l=1:L             %%% block separately
%    ll=(l-1)*d; 
    ll2=(l-1)*Nb;
%    Ma=M(:,ll+1:ll+d);
    Ma2=A*diag(D(l,:))*A';
    W=(Ma2\A).*(ones(d,1)*D(l,:)); %%% demixing matrix in the l-th block    
%    Rm=x(:,ll2+1:ll2+Nb)*x(:,ll2+1:ll2+Nb)'/Nb;
%    W=(Ma2\A)*diag(sqrt(D(l,:)'./diag((A'/Ma2)*Rm*(Ma2\A)))); 
    shat(:,ll2+1:ll2+Nb)=W'*x(:,ll2+1:ll2+Nb);
    if nargin>5
        G=(W'*Atrue).^2; signal=signal+G*diag(sum(Strue(:,ll2+1:ll2+Nb).^2,2));
    end
    num=sum(W.*A);
    SIR(l,:)=num./(1-num);
    if min(SIR(l,:))<0
       disp('one SIR is negative');
    end    
end 
avSIR=10*log10(mean(SIR));
SIR=10*log10(SIR);
%semilogy(iter,'g')


% % % % minimum SDR beamformer
%A=Atrue;
if ~exist('ARorder','var'), ARorder=3; end
Nb2=Nb/10;  % separation block-by-block of the length Nb2
I=eye(d);
signal1=zeros(r,r);
Shat1=zeros(r,N);      
dd=zeros(r,ARorder);
x=[x zeros(d,Nb+ARorder-1)];
RR=zeros(ARorder*d);
if nargin>5, Strue=[Strue zeros(r,ARorder-1)];end
for l=1:N/Nb2              
    ll2=(l-1)*Nb2;
    xx=x(:,max([ll2-Nb/2+1 1]):min([ll2+Nb/2-1 N])); % the estimation of covariance is done on block of length Nb
   % xx=x(:,ll2+1:ll2+Nb2);
    Rm=xx*xx'/Nb;    
    %C=Rm+mean(abs(Rm(:)))*I;
    %C=Rm+0.01*I;
    C=I;
    AAh=khatrirao(A,C\A);
    AAAh=(AAh'*AAh)\AAh';
    dd(:,1)=AAAh*reshape(C\Rm,d^2,1);  
    dd(:,1)=abs(dd(:,1));
    dd(dd<0.001)=0.001; % the variance must be positive
    for i=2:ARorder
        Rm=xx(:,1:end-i+1)*xx(:,i:end)';
        Rm=(Rm+Rm')/(Nb-i+1)/2;
        dd(:,i)=AAAh*reshape(C\Rm,d^2,1);        
    end
    if ARorder>1, dd=armodel(dd')'; end % the correction of AR model due to stability
    for i=1:ARorder
        for j=1:ARorder
            RR((i-1)*d+1:i*d,(j-1)*d+1:j*d)=A*diag(dd(:,abs(i-j)+1))*A';
        end
    end
    for i=1:r
       W=RR\kron(dd(i,:)',A(:,i)); 
       aux=zeros(r,Nb2);
       for j=1:ARorder
           Shat1(i,ll2+1:ll2+Nb2)=Shat1(i,ll2+1:ll2+Nb2)+W((j-1)*d+1:j*d)'*x(:,ll2+j:ll2+Nb2+j-1);
           if nargin>5 % signal energy decomposition due to exact SIR computation
               G=W((j-1)*d+1:j*d)'*Atrue; 
               aux=aux+diag(G)*Strue(:,ll2+j:ll2+Nb2+j-1);           
           end
       end
       if nargin>5 % signal energy decomposition due to exact SIR computation
           signal1(i,:)=signal1(i,:)+sum(aux.^2,2)';
       end
    end
end 


% % % % minimum SDR a la BARBI-II
% ARorder=1;
% I=eye(d);
% signal1=zeros(r,r);
% Shat1=zeros(r,N);      
% DD=zeros(r,ARorder);%dd=zeros(r,ARorder);
% x=[x zeros(d,ARorder-1)];
% RR=zeros(ARorder*d);
% AA=khatrirao(A,A);
% AAA=(AA'*AA)\AA';
% for i=1:ARorder
%     Rm=x(:,1:end-i+1)*x(:,i:end)';
%     Rm=(Rm+Rm')/N/2;
%     DD(:,i)=AAA*reshape(Rm,d^2,1);        
% end
% DD=diag(1./DD(:,1))*DD;
% if ARorder>1, DD=armodel(DD')';end
% Nb2=Nb/10;   % separace po blocich delky Nb2 
% if nargin>4, Strue=[Strue zeros(r,ARorder-1)];end
% for l=1:N/Nb2             
%     ll2=(l-1)*Nb2;
%     xx=x(:,max([ll2-Nb/2+1 1]):min([ll2+Nb/2-1 N])); % odhad statistik je proveden na bloku delky Nb
%     Rm=xx*xx'/length(xx);
%     C=Rm+mean(abs(Rm(:)))*I; % regularizace
%     AAh=khatrirao(A,C\A);
%     AAAh=(AAh'*AAh)\AAh';
%     dd1=AAAh*reshape(C\Rm,d^2,1);  
%     dd1=abs(dd1); dd1(dd1<0.001)=0.001; % rozptyl by mel byt odrazeny od nuly
%    % dd1=D(l,:)'; % puvodni odhad
%     dd=diag(dd1)*DD;
%     for i=1:ARorder
%         for j=1:ARorder
%             RR((i-1)*d+1:i*d,(j-1)*d+1:j*d)=A*diag(dd(:,abs(i-j)+1))*A';
%         end
%     end
%     for i=1:r
%        W=RR\kron(dd(i,:)',A(:,i)); 
%        aux=zeros(r,Nb2);
%        for j=1:ARorder
%            Shat1(i,ll2+1:ll2+Nb2)=Shat1(i,ll2+1:ll2+Nb2)+W((j-1)*d+1:j*d)'*x(:,ll2+j:ll2+Nb2+j-1);
%            if nargin>4
%                G=W((j-1)*d+1:j*d)'*Atrue; 
%                aux=aux+diag(G)*Strue(:,ll2+j:ll2+Nb2+j-1);           
%            end
%        end
%        if nargin>4
%            signal1(i,:)=signal1(i,:)+sum(aux.^2,2)';
%        end
%     end
% end 


% % cubic interpolation of variance
% Shat1=zeros(r,N);
% signal1=zeros(r,r);
% DD=zeros(r,N);
% for i=1:r
%     DD(i,:)=pchip([1 floor(Nb/2):Nb:N N],[D(1,i); D(:,i); D(end,i)],1:N);
% end
% for i=1:N
%     dd=DD(:,i);
%     Rm=A*diag(dd)*A';
%     W=Rm\A*diag(dd);
%     Shat1(:,i)=W'*x(:,i);
%     if nargin>4
%         G=(W'*Atrue).^2; signal1=signal1+G*diag(Strue(:,i).^2);
%     end
% end

 
% l1 sparsity optimization beamformer (using l1magic package)
% Shat1=zeros(r,N);
% signal1=zeros(r,r);
% Apinv=pinv(A);
% for i=1:N
%     Shat1(:,i) = l1eq_pd(Apinv*x(:,i), A, [], x(:,i), 1e-3);
% %    G=(W'*Atrue).^2; signal=signal+G*diag(Strue(:,i+floor(Nb/2)).^2);
% end


%evaluation
if nargin>5
    Strue=Strue(:,1:end-ARorder+1);
    [sir sdr order]=evalSIR(signal,Strue,shat);
    [sir1 sdr1]=evalSIR(signal1,Strue,Shat1);
    disp('          SIR [dB]')
    disp('maxSIRbeamfrm maxSDRbeamfrm')
    disp(10*log10([sir sir1]))
    disp('          SDR [dB]')
    disp('maxSIRbeamfrm maxSDRbeamfrm')
    disp(10*log10([sdr sdr1]))
    disp('Order of component')
    disp(order')
end

end % end of UDSEP7


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A A0 D iter2]=LM1sb(X,r,numit)
%
% computes parallel factor analysis
%
if nargin<3
   numit=100;
end  
d2=exp(100/numit*log(0.5));
%[A D]=zacsym(X,r);
%D=dejC(X,A,A);
[Ia Ibc]=size(X);
Id=floor(Ibc/Ia);
Y=reshape(X,Ia,Ia,Id);
[A A2 D1]=sobium2core(Y,r);
A0=A;
[D]=dejD(X,A,X);
if min(D(:))<0
 %  'INITIALIZATION FAILED'
   D=abs(D);
end
it=0; nu=2; tau=1; %len=r*(Ia+Id); incr=10;
alpha=1;
iter2=chybajina(X,A,D);
%incr=max(D(:));
incre=max(D')/10;
%C=X+incr*kron(ones(1,Id),eye(Ia));
C=X+kron(incre,eye(Ia));
for d=1:Id
    inp=(d-1)*Ia; 
    C(:,inp+1:inp+Ia)=inv(C(:,inp+1:inp+Ia));
end     
iter=chybajina2(X,A,D,C,alpha); % iter2=chybasym(X,A,D);
%[F R]=moje6(X,A,D,C,alpha);
%mu=tau*max(diag(F));
mu=tau*maxdi(A,D,C,alpha);
%F=2*F;
%h=waitbar(0,'iterating');
while it<numit %& incr>1e-6
   it=it+1;
   %if it > 1
   %[F R]=moje6(X,A,D,C,alpha); % viz udsep4
   %F=2*F;
   [dA dD gA gD]=moje9(X,A,D,C,alpha,mu);
   %end
   %d=(F+mu*eye(len))\R;
   err=chybajina2(X,A,D,C,alpha);
   A1=A; D1=D; 
   A=A1+dA; D=D1+dD;
   am=sqrt(sum(A.^2));
   A=A*diag(1./am); D=D*diag(am.*am); 
   if min(D(:))<0
      err2=1e30;
   else   
   err2=chybajina2(X,A,D,C,alpha); % iter2=[iter2 chybasym(X,A,D)];
   end
   lin=sum(sum(dA.*(gA+mu*dA)))+sum(sum(dD.*(gD+mu*dD)));
   %rho=(err-err2)/(d'*(R+mu*d));
   rho=(err-err2)/lin;
   if err2>err
      mu=mu*nu; nu=2*nu; err2=err;
      A=A1; D=D1;
   else
      nu=2; 
      mu=mu*max([1/3 1-(2*rho-1)^3]);
   end   
   if rem(it,5)==0 
      incre=incre*d2; 
      alpha=alpha*d2;
      C=X+kron(incre,eye(Ia));
      for d=1:Id
          inp=(d-1)*Ia; 
          C(:,inp+1:inp+Ia)=inv(C(:,inp+1:inp+Ia));
      end    
   end   
%   incr=(err-err2)/err;
   iter=[iter err2];
   iter2=[iter2 chybajina(X,A,D)];
  % waitbar(it/numit,h)
end   
%close(h)
end % of LM1sb
%iter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err=chybajina(X,A,D,C)
%
% computes an error of approximation of X by sum
% of outer products of columns of A,B and C
%
[Id r]=size(D);
[Ia r]=size(A);
if nargin<4
   for d=1:Id
       inp=(d-1)*Ia; 
       C(:,inp+1:inp+Ia)=inv(X(:,inp+1:inp+Ia));
   end 
end   
err=0;
for d=1:Id 
     inp=(d-1)*Ia;
     Aux=C(:,inp+1:inp+Ia)*(X(:,inp+1:inp+Ia)-A*diag(D(d,:))*A');
     err=err+sum(sum(Aux'.*Aux));
end
end  % of chybajina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err=chybajina2(X,A,D,C,alpha)
%
% computes an error of approximation of X by sum
% of outer products of columns of A,B and C
%
[Id r]=size(D);
[Ia r]=size(A);
err=chybajina(X,A,D,C)-alpha*sum(sum(log(D)))-alpha*Id*sum(log(sum(A.^2)));
end % of chybajina2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dA dD R4a R4d]=moje9(X,A,D,C,alpha,mu)
%
% Computes the Hessian and the gradient of the criterion with the
% barrier
%
[Ia r]=size(A);
[Id r]=size(D);
V=zeros(r,Ia*Id);
W=zeros(r,r*Id);
V2=zeros(Id,Ia*r);
W2=zeros(Id,r*r);
Z2=zeros(Id,Ia*Ia);
H=zeros(Ia,Ia*Id);
for d=1:Id
     inp=(d-1)*Ia; inp1=(d-1)*r;
 %    Za=inv(C(:,inp+1:inp+Ia));
     V(:,inp+1:inp+Ia)=A'*C(:,inp+1:inp+Ia);
     W(:,inp1+1:inp1+r)=A'*C(:,inp+1:inp+Ia)*A;
     H(:,inp+1:inp+Ia)=C(:,inp+1:inp+Ia)*(X(:,inp+1:inp+Ia)-A*diag(D(d,:))*A')*C(:,inp+1:inp+Ia);
     V2(d,:)=reshape(diag(D(d,:))*V(:,inp+1:inp+Ia),1,Ia*r);
     W2(d,:)=reshape(diag(D(d,:))*W(:,inp1+1:inp1+r)*diag(D(d,:)),1,r*r);
     Z2(d,:)=reshape(C(:,inp+1:inp+Ia),1,Ia^2);
 %    Z2(d,:)=reshape(inv(C(:,inp+1:inp+Ia)),1,Ia*Ia);
end   
roh1=V2'*V2;
roh2=Z2'*W2;
Z=roh1; 
for i=1:Ia
    ind=(i-1)*r;
    ina=(i-1)*Ia;
    for j=1:Ia
        indj=(j-1)*r;
        Z(indj+1:indj+r,ind+1:ind+r)=2*(roh1(indj+1:indj+r,ind+1:ind+r)'+...
        reshape(roh2(ina+j,:),r,r)');
    end
end
for i=1:r
    Z(i:r:Ia*r,i:r:Ia*r)=Z(i:r:Ia*r,i:r:Ia*r)-alpha*Id*(eye(Ia)-2*A(:,i)*A(:,i)');
end  %
R4a=zeros(Ia,r);
R4d=zeros(Id,r);
for i=1:r
  R4a(:,i)=alpha*Id*A(:,i);
  for k=1:Id
      inp=(k-1)*Ia;
      R4a(:,i)=R4a(:,i)+2*D(k,i)*H(:,inp+1:inp+Ia)*A(:,i);
  end        
  Aux=reshape(A(:,i)'*H,Ia,Id);
  R4d(:,i)=Aux'*A(:,i)+alpha./D(:,i);
end  
%
Dm=zeros(r,r*Id); DmG=zeros(Id*r,1); DmC=zeros(Ia*r,Id*r);
FG3=zeros(Ia*r,Id*r);
for i=1:Id
    ini=(i-1)*r; ini2=(i-1)*Ia;
    Dm(:,ini+1:ini+r)=inv(W(:,ini+1:ini+r).^2+diag(0.5*(alpha./D(i,:).^2+mu)));
    DmG(ini+1:ini+r,1)=Dm(:,ini+1:ini+r)*R4d(i,:)';
    FG3A=zeros(r*Ia,r);
    for j=1:Ia
        inj=(j-1)*r;
    %    FG3A(inj+1:inj+r,:)=diag(D(i,:))*W(:,ini+1:ini+r)*diag(V(:,ini2+j));
    %     Aux=diag(D(i,:))*W(:,ini+1:ini+r)*diag(V(:,ini2+j));
        FG3A(inj+1:inj+r,:)=(V(:,ini2+j)*D(i,:))'.*W(:,ini+1:ini+r);
    end    
    DmC(:,ini+1:ini+r)=2*FG3A*Dm(:,ini+1:ini+r);
    FG3(:,ini+1:ini+r)=2*FG3A;
%    Z=Z-DmC(:,ini+1:ini+r)*FG3A';
end
Z=Z-DmC*FG3';
ZBz=(Z+mu/2*eye(Ia*r))\(reshape(R4a',Ia*r,1)-DmC*reshape(R4d',r*Id,1));
dA=0.5*reshape(ZBz,r,Ia)';
dD=0.5*reshape(DmG-DmC'*ZBz,r,Id)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D]=dejD(X,A,C)
%
% computes columns of  D in parallel factor analysis X=[A A D];
% using special weight matrices.
%
[Ia r]=size(A);
[Ia Ibc]=size(X);
Id=floor(Ibc/Ia);
D=zeros(Id,r);
for k=1:Id
    ind=(k-1)*Ia;
    M=C(:,ind+1:ind+Ia)\A;
    D(k,:)=(A'*M).^2\diag(M'*X(:,ind+1:ind+Ia)*M);
end
end % of dejD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val=maxdi(A,D,C,alpha)
%
% Computes the maximum diagonal element of the Hessian 
%
[Ia r]=size(A);
[Id r]=size(D);
V=zeros(r,Ia*Id);
W=zeros(r,r*Id);
%V2=zeros(Id,Ia*r);
W2=zeros(Id,r);
roh1b=zeros(1,Ia*r);
Z2=zeros(Id,Ia);
%H=zeros(Ia,Ia*Id);
for d=1:Id
     inp=(d-1)*Ia; inp1=(d-1)*r;
     V(:,inp+1:inp+Ia)=A'*C(:,inp+1:inp+Ia);
     W(:,inp1+1:inp1+r)=A'*C(:,inp+1:inp+Ia)*A;
  %   V2(d,:)=reshape(diag(D(d,:))*V(:,inp+1:inp+Ia),1,Ia*r);
     roh1b=roh1b+reshape(diag(D(d,:))*V(:,inp+1:inp+Ia),1,Ia*r).^2;
     W2(d,:)=diag(diag(D(d,:))*W(:,inp1+1:inp1+r)*diag(D(d,:)))';
     Z2(d,:)=diag(C(:,inp+1:inp+Ia))';
end   
%roh1=V2'*V2;
%roh1=sum(V2.^2)';
%roh2=Z2'*W2;
%roh2=Z2(:,1:Ia+1:end)'*W2;
roh2=Z2'*W2;
val=-1e20;
for i=1:r
  roh=roh1b(1,i:r:r*Ia)'+roh2(:,i);
  dFA=2*roh-Id*alpha*(1-2*A(:,i).^2); 
  dFD=W(i,i:r:r*Id).^2'+0.5*alpha./D(:,i).^2;
  val=max([val; dFA; dFD]);
end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rs=armodel(R,rmax)
%
% to compute AR coefficients of the sources given covariance functions 
% but if the zeros have magnitude > rmax, the zeros are pushed back.
%
[M,d]=size(R);  AR=R;
if nargin<2, rmax=0.99; end
for id=1:d
    AR(:,id)=[1; -toeplitz(R(1:M-1,id),R(1:M-1,id)')\R(2:M,id)];
    v=roots(AR(:,id)); %%% mimicks the matlab function "polystab"
    vs=0.5*(sign(abs(v)-1)+1);
    v=(1-vs).*v+vs./conj(v);
    vmax=max(abs(v));
    if vmax>rmax
       v=v*rmax/vmax;
    end   
    AR(:,id)=real(poly(v)'); %%% reconstructs back the covariance function
end 
Rs=ar2r(AR);
Rs=Rs*diag(R(1,:)./Rs(1,:));
end %%%%%%%%%%%%%%%%%%%%%%%  of armodel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ r ] = ar2r( a )
%%%%%
%%%%% Computes covariance function of AR processes from 
%%%%% the autoregressive coefficients using an inverse Schur algorithm 
%%%%% and an inverse Levinson algorithm (for one column it is equivalent to  
%%%%%      "rlevinson.m" in matlab)
% 
  if (size(a,1)==1)
      a=a'; % chci to jako sloupce
  end
  
  [p m] = size(a);    % pocet vektoru koef.AR modelu
  alfa = a;
  K=zeros(p,m);
  p = p-1;
  for n=p:-1:1
      K(n,:) = -a(n+1,:);
      for k=1:n-1
          alfa(k+1,:) = (a(k+1,:)+K(n,:).*a(n-k+1,:))./(1-K(n,:).^2);
      end
      a=alfa;
  end
%  
  r  = 1./prod(1-K.^2);
  f=[r(1,:); zeros(p,m)];
  b=f;
  for k=1:p 
      for n=k:-1:1
          f(n,:)=f(n+1,:)+K(n,:).*b(k-n+1,:);
          b(k-n+1,:)=-K(n,:).*f(n+1,:)+(1-K(n,:).^2).*b(k-n+1,:);
      end
      b(k+1,:)=f(1,:);
      r=[r; f(1,:)]; 
  end       
end %%%%%%%%%%%%%%%%%%%%%%%%%%%  of ar2r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B,C,niter,time,Fit]=sobium2core(X,R,tol,max_iter)
%-------------------------------------------------------------------------- 
% CANDECOMP/PARAFAC BY SIMULTANEOUS MATRIX DIAGONALIZATION (THIRD-ORDER CASE)
%-------------------------------------------------------------------------
% Description:          
% The CANDECOMP/PARAFAC problem is reformulated as a square simultaneous matrix 
% diagonalization problem, which is solved using cp_extQZ. 
% The SOBIUM2 uniqueness conditions have to be satisfied. 
% In particular, it is necessary that R <= K and R*(R-1) <= I*(I-1)*J*(J-1)/2 (complex data).
% The role of I,J,K can of course be interchanged.
%------------------------------------------------------------------
% INPUTS: 
% X: tensor of size (IxJxK)
% R: number of rank-1 components
%------------------------------------------------------------------
% Optional parameters:
% - tol (default 1e-4): threshold value to stop the algorithm
% - max_iter (default 200): max number of iterations 
%------------------------------------------------------------------
% OUTPUTS:
% - Estimated loading matrices: A (IxR), B(JxR), C(KxR)
% - niter: number of iterations to converge 
% - time: time to converge 
% - Fit: Frobenius Norm of the residual tensor
%-------------------------------------------------------------------
%
% Refs: L. De Lathauwer, J. Castaing, ``Blind Identification of
% Underdetermined Mixtures by Simultaneous Matrix Diagonalization'',
% IEEE Trans. Signal Processing, Vol. 56, No. 3, March 2008, pp. 1096-1105.
% 
% Feed-back: lieven.delathauwer@kuleuven-kortrijk.be
%
% Noncommercial use only
%
% Version October 1, 2008
% Copyright Lieven De Lathauwer
%

%---------------------------------------------------------------------
% Check inputs and define default parameters
%---------------------------------------------------------------------
[I1, J1, K1]=size(X);
[size_sort,perm_vec]=sort([I1 J1 K1],'ascend');               
       
%--- Permute tensor if necessary such that the long dimension is the third one
[size_sort,perm_vec]=sort([I1 J1 K1],'ascend'); 
X=permute(X,perm_vec);
[I,J,K]=size(X);       % K is now the longest dimension
   
%----------------
if nargin==1
   tol=1e-4;max_iter=200;
   Y = reshape(X, I*J, K);      
   singval=svd(Y);
   disp('Please enter the rank R of the tensor by inspection of the following singular values');
   singval
   R=input('Which rank do you choose?  ');
end 
while R==[]
  R = input(['Please enter the rank R of the tensor','\n']);
end
% if R*(R-1)>I*(I-1)*J*(J-1)/2
%     error('Error: R has to bounded by R*(R-1) <= I*(I-1)*J*(J-1)/2')	
% end
if R>K
    error('Error: R should not be greater than K')	
end
%-----------------
if nargin==2
    tol=1e-4;max_iter=200;
end
%------------------
if nargin==3
    max_iter=200;
end
%-----------------
if isreal(X); rc=0; else;rc=1; end

param=[max_iter,tol,rc];  

%-------------------------------------------------------------------
%  Now compute the PARAFAC decomposition of X
%-------------------------------------------------------------------
tic

        %%%%%%%%%%% step 1 %%%%%%%%%%
        %disp('step 1: best rank-R approximation of (JI x K) matrix representation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        Y = reshape(X, I*J, K);      % J varying more slowly             
        %[U_y,S_y,V_y]	= svds(Y,R);   % compute only the R first singular vectors and values
        [U_y,S_y,V_y]	= svd(Y,0);
        U_y=U_y(:,1:R);
        S_y=S_y(1:R,1:R);
        V_y=V_y(:,1:R);
        E_y = U_y*S_y;     
        Y=E_y*V_y';
        
        %%%%%%%%%%% step 2 %%%%%%%%%%
        %disp('step 2: Computation of P')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Create matrix P=[Phi2(Es,Et)] with (s, t) in [1, R] 
        P=zeros(size(E_y,1)^2,R*(R+1)/2);
        egal_ind=0;
        diff_ind=0;

        for t=1:R
            Et=unvec(E_y(:,t),J);
            for s=1:t
                Es=unvec(E_y(:,s),J);
                Pst=Phi3(Es, Et);
                if s==t
                    egal_ind=egal_ind+1;
                    P(:,egal_ind) = Pst;    % R terms
                else
                    diff_ind=diff_ind+1;
                    P(:,R+diff_ind)  = 2*Pst;   % R*(R-1)/2 terms
                end
            end
        end
        %%%%%%%%%%% step 3 %%%%%%%%%%
        %disp('step 3: Estimation of the set of R matrices Wr to be diagonalized')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Solve the equation sum_{s,t} (Phi(s,t) W(s, t))=0 with W symmetric

        [U_p, S_p, V_p]	= svd(P,0);
        W_inter		= V_p(:, end:-1:end-R+1) ; %% W_inter belongs to the kernel of P, size (R(R+1)/2xR)

        W=zeros(R, R, R);
        for i3=1:R    
            W_diag(:,:, i3)=diag(W_inter(1:R, i3));  % Fill diagonal: the first R columns of P corresponds to s==t
            r=R+1;
            for i2=1:R                               % Fill upper triangular part
                for i1=1:i2-1
                    W(i1,i2,i3)=W_inter(r,i3);
                    r=r+1;
                end     
            end    
        end
        W=W+permute(W, [2,1,3])+W_diag;   % Build the symmetric matrix
        
        
        %%%%%%%%%%% step 4 %%%%%%%%%%%%%%%%%%%%%%%%%
        %disp('step 4: Simultaneous Diagonalization')
        % Use extended QZ iteration 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %[niter,F_est]	= cp_extQZ(W, param); %% Solve Wi = F Di F.'
        iF_est=uwedge(reshape(W,R,R*R));
        %%%%%%%%%%% step 5 %%%%%%%%%%
        %disp('step 5: Estimation of A, B, C') by using Khatri-Rao
        % product properties, i.e. each column results from the Kronecker
        % product of 2 vectors so it is the vectorized representation of
        % a rank-1 matrix --> take the principal left and right sing. vectors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	        
        %C	= conj(V_y)*inv(F_est).';
        %P_est	= E_y*F_est;   % P_est is close to Khatri_Rao(B, A)
        %C	= conj(V_y)*iF_est.';
        P_est	= E_y*inv(iF_est);   % P_est is close to Khatri_Rao(B, A)

        for r=1:R
   %         Pr_est(:,:,r) 		= reshape(P_est(:,r), I, J);   
   %         [Ar_est,Dr_est,Br_est]	= svd(Pr_est(:,:,r),0);  % best rank-1 approx
            [Ar_est,Dr_est,Br_est]	= svd(reshape(P_est(:,r), I, J),0);  
            A(:,r)		= Ar_est(:,1);
            B(:,r)		= conj(Br_est(:,1));
        end

        C = Y.'/kat_rao(B,A).';
      
        Fit=norm(reshape(X, I*J, K)-kat_rao(B,A)*C.','fro');
      
  
time=toc;

%---------------------------------------------------------------
% Permute the order of the loading matrices according
%---------------------------------------------------------------
L=cell(1,3);L{1}=A;L{2}=B;L{3}=C;
L(:,perm_vec)=L;
A=L{1};B=L{2};C=L{3};

end % of sobium2core 
%*************** END OF CORE FUNCTION *****************************
%*******************************************************************
      
%*******************************************************************        
function C = kat_rao(A,B)
%-------------------------------------------------------------------
% Calculate the Khatri-Rao product of 2 matrices
%-------------------------------------------------------------------
        I=size(A,1);R1=size(A,2);
        J=size(B,1);R2=size(B,2);
        if R1~=R2
            error('Input matrices must have the same number of columns for Khatri-Rao product')
        end
        
        C=zeros(I*J,R1);
        for j=1:R1
            C(:,j)=reshape(B(:,j)*A(:,j).',I*J,1);
        end
        
        % Note: this is much faster than the following loop, because kron is avoided
%         for j=1:R1
%              C(:,j)=kron(A(:,j),B(:,j));     % column-wise Kronecker product
%         end
end % of kat_rao 
       
%*************************************************************************
function S=Phi3(X, Y)
        % Input:		X(IxJ)	Y(KxL)
        %
        % Output:	S= vector representation of Phi(X, Y) defined by
        %		Phi(X, Y)(i,j,k,l) = xij ykl + yij xkl - xkj yil - ykj xil
        %
        [I, J]=size(X);
        [K, L]=size(Y);

        Y2=permute(Y,[3 4 1 2]);
        tens_X=repmat(X,  [1 1 K L]);
        tens_Y=repmat(Y2, [I J 1 1]);
        A=tens_X.*tens_Y;   %A(i,j,k,l)=X(i,j)*Y(k,l)
        
        X2=permute(X,[3 4 1 2]);
        tens_Y=repmat(Y,  [1 1 I J]);
        tens_X=repmat(X2, [K L 1 1]);
        B=tens_X.*tens_Y;   %B(i,j,k,l)=Y(i,j)*X(k,l)

        P1=A+B;
        P2=permute(P1, [3 2 1 4]);
        D=P1-P2;
        S=reshape(D, numel(X)^2,1,1,1);
end % of Phi3
%*************************************************************************
function [W_est Ms]=uwedge(M,W_est0)
%
% Uniformly Weighted Exhaustive Diagonalization using Gauss itEration (U-WEDGE)
%
% Coded by Petr Tichavsky, March 2008, updated July 2008
%
% Please cite
%
%  P. Tichavsky, A. Yeredor and J. Nielsen,
%     "A Fast Approximate Joint Diagonalization Algorithm
%      Using a Criterion with a Block Diagonal Weight Matrix",
%      ICASSP 2008, Las Vegas
%
% Input: M .... the matrices to be diagonalized, stored as [M1 M2 ... ML]
%        West0 ... initial estimate of the demixing matrix, if available
%
% Output: W_est .... estimated demixing matrix
%                    such that W_est * M_k * W_est' are roughly diagonal
%         Ms .... diagonalized matrices composed of W_est*M_k*W_est'
%         crit ... stores values of the diagonalization criterion at each
%                  iteration
%
[d Md]=size(M);
L=floor(Md/d);
Md=L*d;
iter=0;
eps=1e-7;
improve=10;
if nargin<2
   [H E]=eig(M(:,1:d));
   W_est=diag(1./sqrt(abs(diag(E))))*H';
else
   W_est=W_est0;
end 
Ms=M; 
Rs=zeros(d,L);
for k=1:L
      ini=(k-1)*d;
      M(:,ini+1:ini+d)=0.5*(M(:,ini+1:ini+d)+M(:,ini+1:ini+d)');
      Ms(:,ini+1:ini+d)=W_est*M(:,ini+1:ini+d)*W_est';
      Rs(:,k)=diag(Ms(:,ini+1:ini+d));
end
crit=sum(Ms(:).^2)-sum(Rs(:).^2);
C1=zeros(d,d);
while improve>eps && iter<15
  B=Rs*Rs';
  for id=1:d
      C1(:,id)=sum(Ms(:,id:d:Md).*Rs,2);
  end
  D0=B.*B'-diag(B)*diag(B)';
  A0=eye(d)+(C1.*B-diag(diag(B))*C1')./(D0+eye(d));
  W_est=A0\W_est;
  Raux=W_est*M(:,1:d)*W_est';
  aux=1./sqrt(abs(diag(Raux)));
  W_est=diag(aux)*W_est;  % normalize the result
  for k=1:L
     ini=(k-1)*d;
     Ms(:,ini+1:ini+d) = W_est*M(:,ini+1:ini+d)*W_est';
     Rs(:,k)=diag(Ms(:,ini+1:ini+d));
  end
  critic=sum(Ms(:).^2)-sum(Rs(:).^2);
  improve=abs(critic-crit);
  crit=critic;
  iter=iter+1;
end 
end  % of uwedge
%*******************************************************************
function V = unvec(v,n)
        % function V = unvec(v,n)
        %
        % form matrix V from column vector v, where V has n columns.
        % If n is not specified, V is made square

        [m,N] = size(v);
        if (nargin == 1), 
           n = sqrt(m);
        end

        m = m/n;
        V = zeros(m,n);
        for i = 0:n-1,
            V(:,i+1) = v(i*m+1:i*m+m);
        end
end % of unvec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SIR SDR p]=evalSIR(signal,S,shat)
% exact SIR and SDR evaluation
% signal .... signal energy decomposition matrix
% S ......... true sources
% shat ...... sources' estimates
signal=diag(1./sqrt(mean(signal,2)))*signal;

d=size(signal,1);
SIR=zeros(d,1);
SDR=zeros(d,1);
p=zeros(d,1);
for i=1:d
    sir=signal(:,i)./sum(signal(:,[1:i-1 i+1:d]),2);
    [sir, ind]=max(sir);
    SIR(i)=sir;
    p(i)=ind;    
    CC=mean(S(i,:).*shat(ind,:))/mean(shat(ind,:).^2);
    SDR(i)=mean(S(i,:).^2)/mean((S(i,:)-CC*shat(ind,:)).^2);
end
end