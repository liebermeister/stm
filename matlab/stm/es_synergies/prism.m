function [Z,Conflict,Direct]=prism(ECM,Alpha,Gamma,PrismR)
%
% function [Z,Conflict,Direct]=prism(ECM,Alpha,Gamma,PrismR)
% 
% Cluster a two-color network into monochromatically interacting clusters
% 
% ECM: A symmetric square matrix describing the interaction network
% ECM(i,j): Interaction between nodes i and j: 
% (0) No interaction; (1) Positive ; (-1) Negative ;
%
% Alpha: Associative (0) versus Direct (1) contributions to cluster affinities. 
% Set between 0 and 1. 
%
% Gamma: Raise associative correlation to a power for increasing the 
%        importance of high similarity (default is 1 - no enhancement). 
%
% PrismR: Run in Prism (0, default) or in Prism-R (1) mode. 
%
% The function returns the clustering in Z (see help linkage)
% Also returns:
%    Conflict: the number of conflicts generated in each clustering step
%    where, Qmodule = sum((Conflict>0).*Conflict) 
%    Direct: the direct postive and negative connections between the clusters 
%
% 
% "Modular epistasis in yeast metabolism"
% Daniel Segrè, Alexander DeLuna, George M. Church, Roy Kishony
% 
%
% For questions or comments please contact:
%
% Roy Kishony
% Bauer Center for Genomics Research
% Harvard University
% 7 Divinity Ave
% Cambridge, MA 02138
% E-mail: rkishony@cgr.harvard.edu
%
%
% R.K. 4/24/04
%


if nargin<3
    Gamma=1;
end
if nargin<4
    PrismR=0;
end

Y=zeros(size(ECM,1),size(ECM,2),2) ;
Y(:,:,1) = (ECM== 1) ;
Y(:,:,2) = (ECM==-1) ;

warning off MATLAB:divideByZero

m=size(Y,1) ;
Z = zeros(m-1,3); 
Conflict = zeros(m-1,1) ;
Direct = zeros(m-1,2) ;
N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; 
R = 1:n;

for s=1:(n-1)
   om=ones(1,m) ; o12=[1 2] ; o21=[2 1] ;
   
   X_dir = zeros(size(Y)) ;                        
   X_dir(:,:,1) = Y(:,:,1) ./ (N(R)'*N(R)) ;
   X_dir(:,:,2) = Y(:,:,2) ./ (N(R)'*N(R)) ;
   
   if s==1 
      S_ass = zeros(size(Y)) ;          
      S_ass(:,:,1) = (X_dir(:,:,1)^2 + (1-X_dir(:,:,1))^2) / (2*m) ;
      S_ass(:,:,2) = (X_dir(:,:,2)^2 + (1-X_dir(:,:,2))^2) / (2*m) ;   
      S_ass_max=sum(S_ass,3) ;
   end
   
   X = (1-Alpha) * S_ass_max.^Gamma + Alpha * sum(X_dir,3) ;
   X=1.01-X ;

   Y_bln = Y>0 ;

   Y_exc = Y_bln(:,:,o12) .* ~Y_bln(:,:,o21) ;     
   
   S_con = zeros(size(Y)) ;                        
   S_con(:,:,1) = Y_exc(:,:,1) * Y_exc(:,:,2) ;
   S_con(:,:,2) = Y_exc(:,:,2) * Y_exc(:,:,1) ;
   
   Y_exc_d = zeros(size(Y)) ;                      
   dia=diag(Y_exc(:,:,1)) ; diap=dia' ; Y_exc_d(:,:,1)=dia(:,om) + diap(om,:);
   dia=diag(Y_exc(:,:,2)) ; diap=dia' ; Y_exc_d(:,:,2)=dia(:,om) + diap(om,:);
   
   N_con = zeros(size(X)) ;                        
   N_con = N_con - Y_bln(:,:,1) .* Y_bln(:,:,2) ;  
   N_con = N_con + sum( S_con - Y_exc .* Y_exc_d(:,:,o21) , 3 ) ; 
   
   z_up = find(triu(ones(m,m),1)) ;
   N_con_min=min(N_con(z_up)) ;
   z_min_con=z_up(find(N_con(z_up)==N_con_min)) ;
   
   if ~PrismR
       [v,k]=min(X(z_min_con)) ; k=z_min_con(k) ; 
   else
       k=floor(rand*length(z_min_con))+1 ;k=z_min_con(k) ; 
       v=X(k) ;
   end
   i=mod(k-1,m)+1 ;
   j=(k-i)/m + 1 ;
   Z(s,:) = [R(i) R(j) v] ; 
   Conflict(s) = N_con(k) ;
   Direct(s,:) = [X_dir(i,j,1), X_dir(i,j,2)] ;
      
   new_Y = Y(i,:,:) + Y(j,:,:) ;
   new_Y(1,i,:)=new_Y(1,i,:) + Y(j,j,:) + Y(i,j) ;
   Y(i,:,:) = new_Y ; Y(:,i,:) = new_Y ;
   Y(j,:,:) = [] ; Y(:,j,:) = [] ;
   
   S_ass_max_new = max(S_ass_max(i,:),S_ass_max(j,:)) ;
   S_ass_max(i,:) = S_ass_max_new ;
   S_ass_max(:,i) = S_ass_max_new' ;
   S_ass_max(j,:) = [] ; S_ass_max(:,j) = [] ;

   N(n+s) = N(R(i)) + N(R(j));
   R(i) = n+s;
   R(j:(m-1))=R((j+1):m);  R(m)=[];
   m = m-1; 
end
