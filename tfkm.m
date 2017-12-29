% tfkm carries out trimmed factorial k-means algorithm, introduced in
%   
% Farnè, M. and Vouldis, A. T. (2017), 'Business models of the Banks in the Euro Area', 
% No. 2070. ECB Working Paper.
%  
% This method is a robustified version of the factorial k-means algorithm developed
% in Vichi, M., & Kiers, H. A. (2001). Factorial k-means analysis for two-way data. 
% Computational Statistics & Data Analysis, 37(1), 49-64.
%
% The INPUT arguments are: 
% z: n (objects) times p (variables) data matrix.
% alpha: allowed proportion of outliers. alpha=0 means that the original
% factorial k-means is applied.
% r_all: vector of tested number of factors.
% g_all: vector of te numbers of clusters.
% cov_yes=1 if the input is the covariance matrix, 0 if the input is the
% correlation matrix.
% r_G: 0 for a rational initialization of loadings, memberships and
% centroids (default), 1 for a random initialization (see Farnè, M. and
% Vouldis, A. T. (2017) for the details).
% N: number of replicates (defaults to 100).

% The OUTPUT arguments are:
% fac_scores: the matrix of factor scores for each object.
% Member_mod: the vector of group membership (outliers=0).
% group_num: the number of observations in each cluster (excluding
% outliers).
% out_tfkm: outlier identifiers.
% n_out_tfkm: number of identified outliers.
% Member: the vector of group membership (without excluding outliers).
% group_comp: the number of observations in each cluster (without excluding
% outliers).
% Hartigan: the matrix of Hartigan statistics for each rank-number of
% groups combination.
% i1: optimal rank according to the maximization of Hartigan statistics.
% i2: optimal number of groups according to the maximization of Hartigan
% statistics.
% model_select: 1 if rank and n_groups are selected maximizing Hartigan statistics, 0 otherwise. 
% Available only if g_all has length larger than one. In case its length is one while 
% the length of r_all is larger than 1, the first rank of r_all is automatically selected.
% Please note that the length of g_all must be such that Hartigan
% statistics may be computed for all pair (r,n_groups) with r<=n_groups-1.
% If we want to consider such a pair (r,n_groups), n_groups+1 must be
% contained in g_all.

function[output_all]=tfkm(z,alpha,r_all,g_all,cov_yes,r_G,N)

if~exist('r_G','var')
    r_G=0;
end

if~exist('N','var')
    N=100;
end

%if ~exist('alpha','var')
%    alpha=0.1;
%end

if alpha==0%%~exist('robust','var') &&
    robust=0;
end

if ne(alpha,0)==1%%~exist('robust','var') && 
    robust=1;
end

n=size(z,1);%%number of objects
p=size(z,2);%%number of variables
if cov_yes==1
    C=cov(z)%%covariance matrix
else
C=corr(z);
end%%correlation matrix
r_length=length(r_all);
g_length=length(g_all);

if length(g_all)==1
    model_select=0
else
    model_select=1
end

    Hartigan=0;

if model_select==1
    
    if r_all(1) > g_all(g_length)-1
        % the rank can not exceed the number of cluster minus one
        output_all=strcat('The rank is larger then the number of clusters minus one')
        return;
        %quit;
    end

% the algorithm starts
for i1=r_all(1):r_all(r_length)%%rank
    for i2=g_all(1):g_all(g_length)%%number of groups

        if i1 <=i2-1
        
% Initialize rank
[U_r, D_r]=svds(C,i1);

% Initialisation: 
R=zeros(p,p);
V=zeros(p,p);
E=zeros(p,p);
v=zeros(p,1);
Random=zeros(p,p,N);  
r=i1;
p=length(C);

clear Random E V

for number=1:N
    K=rand(p)*eye(p);    
    % Here a Gram-Schmidt algorithm is used to derive the orthogonal base of the space defined by K
    E(:,1)=K(:,1);  
        
    for j=2:p
    for i=1:(j-1)
        R(i,j)=dot(K(:,j),E(:,i))/(norm(E(:,i))^2); 
    end
    for i=1:(j-1)
    V(:,i)=R(i,j)*E(:,i);  
    end
   
        for h=1:p
        v(h,1)=sum(V(h,1:(j-1)));  
        end
       E(:,j)=K(:,j)-v;    
    end

    for i=1:p
    E(:,i)=E(:,i)/norm(E(:,i)); 
    end

    Random(:,:,number)=E*U_r;
end

% Initialize n_groups

n_groups=i2;
group=1:n_groups;
RawData=z;  
% in case we want to run the clustering in a subset of z, this can be easily done here !!
n_Raw=size(RawData,1);

%% Random group initializer

%r_G=0;
clear G Random_G centr maxRaw minRaw rangeRaw dist_G t_hotelling_pre F_pre pace

%% The clustering initialisation loop
for number=1:N
   
for i=1:p
maxRaw(i)=max(RawData(:,i));      
minRaw(i)=min(RawData(:,i));
rangeRaw(i)=maxRaw(i)-minRaw(:,i);
end

vvar=var(z');
mmean=mean(z');

if r_G==1
for j=1:n_Raw
    casual(j)=randi(n_groups); 
    G(j,casual(j))=1;
end
end

if r_G==0
    fac_pre=z*Random(:,:,number);
for j=1:n_Raw
     dist_G(:,j)=(fac_pre(j,:)-mean(fac_pre))';
     t_hotelling_pre(:,j)=n_Raw*dist_G(:,j)'*inv(cov(fac_pre))*dist_G(:,j);
end   
    
  q_G=quantile(t_hotelling_pre,2*n_groups-1);

  % The distances of objects from the various centroids is calculated
  for j=1:n_Raw
  for m=1:n_groups
  F_pre(j,m)=norm(t_hotelling_pre(j)-q_G(2*m-1));
  end
  end
  
% Here the assignment of objects to clusters is done, based on the t-hotelling variance measure

for j=1:n_Raw
for m=1:n_groups
    if min(F_pre(j,:))==F_pre(j,m)
        G(j,m)=1;
    else
        G(j,m)=0;
    end
end
end
end
    
Random_G(:,:,number)=G;

clear G
end

%%

clear YY GG UU_GG F IQR Member

n_it=0;
%robust=1;

for number=1:N%%number of initializers
G=Random_G(:,:,number);
U_G=Random(:,:,number);
Y_bar=inv(G'*G)*G'*RawData*U_G;
ob_pre=norm(RawData*U_G-G*Y_bar)^2;   
diff=ob_pre;
n_it=0;

% The main loop of clustering starts. It runs while diff>0
while ((diff))>0

    % For each object j, for each group m, the loss for the object is calculated
    for j=1:n
        
        for m=1:n_groups
            count=zeros(1,n_groups);
            count(m)=1; 
            F(j,m)=norm(z(j,:)*U_G-count*Y_bar)^2; 
        end
        
        % each object is assigned to clusters according to the minimum distance
        for m=1:n_groups
            if min(F(j,:))==F(j,m)
                G(j,m)=1;
            else
                G(j,m)=0;
            end
        end

    end

clear Member
for j=1:n
    for m=1:n_groups
    if G(j,m)==1
       Member(j)=m;
    end
    end
end

% Performs Step 1 in (Vichi and Kiers, 2001, p. 57), so in case the loop ends, 
% the centroids and groupings are updated

centroids=Y_bar;
fac_scores=(z)*U_G;
cen_true=G*centroids;

if robust==1
for m=1:n_groups
for j=1:n
    diff_scores(:,j)=fac_scores(j,:)'-cen_true(Member(j),:)';
    t_hotelling(j)=n*diff_scores(:,j)'*inv(cov(fac_scores))*diff_scores(:,j);
    t_test(j)=(n-r)/(r*(n-1))*t_hotelling(j);
end
end

% the 100(1-alpha)-th percentile of the test distribution is calculated
thr_clust_max=quantile(t_test,1-alpha);

clear flag
boundmax=thr_clust_max;
%  The objects exceeding the threshold are flagged as outliers
for j=1:n
    if t_test(j)>boundmax 
       flag(j)=1;
    else flag(j)=0;
    end
end

% outliers are removed from the groups
for j=1:n
    if flag(j)==1
       G(j,:)=zeros(1,n_groups);
    end
end
end

if rank(G)==size(G,2)   
    % if rank(G) equals the number of columns in G
[U_G, D_G]=svds(z'*(G*inv(G'*G)*G'-eye(n))*z,i1);
% extraction of eigenvectors in Step 2 
% see ten Berge, J. M. (1993). Least squares optimization in multivariate
% analysis. Leiden, The Netherlands: DSWO Press, Leiden University.

Y_bar=inv(G'*G)*G'*z*U_G;
if n_it>1
ob_pre=ob_post;
end
ob_post= norm(z*U_G-G*Y_bar)^2;
diff=ob_post-ob_pre;
diff_pre=diff;
n_it=n_it+1;
if diff<0
YY(:,:,number)=Y_bar;
GG(:,:,number)=G;
UU_GG(:,:,number)=U_G;
clear obob
obob(number)=ob_post;
diffdiff(number)=diff;
N_it(number)=n_it;
end
else diff=-1;           
    % exits immediately, without saving nothing. goes on with the next initialiser
end
end
if diff==0
   obob(number)=0;
   N_it(number)=0;
end

end

% We choose the solution among the set (=N) of initial starts

if sum(N_it)~=0
gg=find(obob==min(obob(N_it~=0)));
find(obob==min(obob(N_it~=0)));
end

n_it=0;
%robust=1;
%thr_clust=1.5;

% the procedure is repeated only for the best solution (gg) in order to
% have the variables ready for plots and other calculations below!
for number=gg(1):gg(1)
G=Random_G(:,:,number);
U_G=Random(:,:,number);
Y_bar=inv(G'*G)*G'*z*U_G;
ob_pre=norm(z*U_G-G*Y_bar)^2;
diff=ob_pre;
n_it=0;
while ((diff))>0
for j=1:n
for m=1:n_groups
    count=zeros(1,n_groups);
    count(m)=1;
    F(j,m)=norm(z(j,:)*U_G-count*Y_bar)^2; 
end

for m=1:n_groups
    if min(F(j,:))==F(j,m)
        G(j,m)=1;
    else
        G(j,m)=0;
    end
end

end

clear Member
for j=1:n
    for m=1:n_groups
    if G(j,m)==1
       Member(j)=m;
    end
    end
end

centroids=Y_bar;
fac_scores=(z)*U_G;
cen_true=G*centroids;

if robust==1
clear diff_scores t_hotelling
for m=1:n_groups
for j=1:n
    diff_scores(:,j)=fac_scores(j,:)'-cen_true(Member(j),:)';
    t_hotelling(j)=n*diff_scores(:,j)'*inv(cov(fac_scores))*diff_scores(:,j);
    t_test(j)=(n-r)/(r*(n-1))*t_hotelling(j);
end
end

IQR=quantile(t_test,0.75)-quantile(t_test,0.25);
thr_clust_max=quantile(t_test,1-alpha);

clear flag
boundmax=thr_clust_max;
for j=1:n
    if t_test(j)>boundmax
       flag(j)=1;
    else flag(j)=0;
    end
end

for j=1:n
    if flag(j)==1
       G(j,:)=zeros(1,n_groups);
    end
end
end

if rank(G)==size(G,2)
[U_G, D_G]=svds(z'*(G*inv(G'*G)*G'-eye(n))*z,i1);
Y_bar=inv(G'*G)*G'*z*U_G;
if n_it>1
ob_pre=ob_post;
end
ob_post= norm(z*U_G-G*Y_bar)^2;
diff=ob_post-ob_pre;
diff_pre=diff;
n_it=n_it+1;
if diff<0
YY(:,:,number)=Y_bar;
GG(:,:,number)=G;
UU_GG(:,:,number)=U_G;
obob(number)=ob_post;
diffdiff(number)=diff;
end
else diff=-1;
end
end
end

% Hartigan statistics see Hartigan, J. A. (1975). Clustering algorithms (Vol. 209). 
% New York: Wiley.

    for i=1:n
        fac_good=z*U_G;
            diffscores(i)=(norm(fac_good(i,:)'-Y_bar(Member(i),:)'))^2;
                               
    end

H(i1,i2)=sum(reshape(diffscores,[],1));   


i1
i2

clearvars -except gg r_G N robust i1 i2 model_select z r_length g_length r_all g_all alpha C z p n H Hartigan;
        end
    end
end

for k=g_all(1):g_all(g_length-1)
    for r=r_all(1):r_all(r_length)
    if r <= k-1
% the rank can not exceed the number of cluster minus one
       Hartigan(r,k)=((H(r,k)/H(r,k+1))-1)*(p-k-1);
    end 
    end
end
Hartigan

        if size(Hartigan,1)==1
        output_all=strcat('It is impossible to compute Hartigan statistics. Please enlarge the vector g_all or set two specific values for r and n_groups s.t. r<=n_groups-1')
        return;
        end

[min2 ii1]=max(Hartigan)
[min1 i2]=max(min2)
i1=ii1(i2)
  
else
    
    i1=r_all(1)
    i2=g_all
    
    if i1>i2-1
        % the rank can not exceed the number of cluster minus one
        output_all=strcat('The rank is larger then the number of clusters minus one')
        return;
    end
end


% OPTIMAL SOLUTION IMPLEMENTATION
% Initialize rank
[U_r, D_r]=svds(C,i1);         
R=zeros(p,p);
V=zeros(p,p);
E=zeros(p,p);
v=zeros(p,1);
Random=zeros(p,p,N);  
r=i1;
p=length(C);

clear Random E V
for number=1:N
    K=rand(p)*eye(p); 
    rank(K);     
    % Here a Gram-Schmidt algorithm is used to derive the orthogonal base of the space defined by K
    E(:,1)=K(:,1);  
        
    for j=2:p
    for i=1:(j-1)
        R(i,j)=dot(K(:,j),E(:,i))/(norm(E(:,i))^2); 
    end
    for i=1:(j-1)
    V(:,i)=R(i,j)*E(:,i);  
    end
   
        for h=1:p
        v(h,1)=sum(V(h,1:(j-1)));  
        end
       E(:,j)=K(:,j)-v;    
    end
                
     rank(E);
    E'*E;

    for i=1:p
    E(:,i)=E(:,i)/norm(E(:,i)); 
    end
    rank(E);
    E;
    E'*E;  

    Random(:,:,number)=E*U_r;
end

% Initialize n_groups
n_groups=i2;
group=1:n_groups;
RawData=z;
n_Raw=size(RawData,1);
%r_G=0;
clear G Random_G centr maxRaw minRaw rangeRaw dist_G t_hotelling_pre F_pre pace

% The clustering initialisation loop
for number=1:N
   
for i=1:p
maxRaw(i)=max(RawData(:,i));      
minRaw(i)=min(RawData(:,i));
rangeRaw(i)=maxRaw(i)-minRaw(:,i);
end

vvar=var(z');
mmean=mean(z');

if r_G==1
for j=1:n_Raw
    casual(j)=randi(n_groups); 
    G(j,casual(j))=1;
end
end

if r_G==0
    fac_pre=z*Random(:,:,number);
for j=1:n_Raw
     dist_G(:,j)=(fac_pre(j,:)-mean(fac_pre))';
     t_hotelling_pre(:,j)=n_Raw*dist_G(:,j)'*inv(cov(fac_pre))*dist_G(:,j);
end   
    
  q_G=quantile(t_hotelling_pre,2*n_groups-1);

  % The distances of objects from the various centroids is calculated                   
  for j=1:n_Raw
  for m=1:n_groups
  F_pre(j,m)=norm(t_hotelling_pre(j)-q_G(2*m-1));
  end
  end

% Here the assignment of objects to clusters is done, based on the t-hotelling variance measure

for j=1:n_Raw
for m=1:n_groups
    if min(F_pre(j,:))==F_pre(j,m)
        G(j,m)=1;
    else
        G(j,m)=0;
    end
end
end
end
    
Random_G(:,:,number)=G;

clear G
end

clear YY GG UU_GG F IQR Member
n_it=0;
%robust=1;

for number=1:N
G=Random_G(:,:,number);   
U_G=Random(:,:,number);
Y_bar=inv(G'*G)*G'*RawData*U_G;   
ob_pre=norm(RawData*U_G-G*Y_bar)^2;      
diff=ob_pre;
n_it=0;

% The main loop of clustering starts. It runs while diff>0
while ((diff))>0

    for j=1:n
        
        for m=1:n_groups
            count=zeros(1,n_groups);
            count(m)=1; 
            F(j,m)=norm(z(j,:)*U_G-count*Y_bar)^2; 
        end
        
        % each object is assigned to clusters according to the minimum distance
        for m=1:n_groups
            if min(F(j,:))==F(j,m)
                G(j,m)=1;
            else
                G(j,m)=0;
            end
        end

    end

clear Member
for j=1:n
    for m=1:n_groups
    if G(j,m)==1
       Member(j)=m;
    end
    end
end

% Performs Step 1 in (Vichi and Kiers, 2001, p. 57), so in case the loop ends, 
% the centroids and groupings are updated

centroids=Y_bar;
fac_scores=(z)*U_G;
cen_true=G*centroids;

if robust==1
for m=1:n_groups
for j=1:n
    diff_scores(:,j)=fac_scores(j,:)'-cen_true(Member(j),:)';
    t_hotelling(j)=n*diff_scores(:,j)'*inv(cov(fac_scores))*diff_scores(:,j);
    t_test(j)=(n-r)/(r*(n-1))*t_hotelling(j);
end
end

thr_clust_max=quantile(t_test,1-alpha);

clear flag

boundmax=thr_clust_max;

% The objects exceeding the threshold, are flagged as outliers
for j=1:n
    if t_test(j)>boundmax
       flag(j)=1;
    else flag(j)=0;
    end
end

% outliers are removed from the groups
for j=1:n
    if flag(j)==1
       G(j,:)=zeros(1,n_groups);
    end
end
end

if rank(G)==size(G,2)   
[U_G, D_G]=svds(z'*(G*inv(G'*G)*G'-eye(n))*z,i1);
% extraction of eigenvectors in Step 2 
% see ten Berge, J. M. (1993). Least squares optimization in multivariate
% analysis. Leiden, The Netherlands: DSWO Press, Leiden University.

Y_bar=inv(G'*G)*G'*z*U_G;
if n_it>1
ob_pre=ob_post;
end
ob_post= norm(z*U_G-G*Y_bar)^2;
diff=ob_post-ob_pre;
diff_pre=diff;
n_it=n_it+1;
if diff<0
YY(:,:,number)=Y_bar;
GG(:,:,number)=G;
UU_GG(:,:,number)=U_G;
clear obob
obob(number)=ob_post;
diffdiff(number)=diff;
N_it(number)=n_it;
end
else diff=-1;         
end
end
if diff==0
   obob(number)=0;
   N_it(number)=0;
end

end

% We choose the solution among the set (=N) of initial starts.
% The solution with the minimum loss function is chosen.


if sum(N_it)~=0
gg=find(obob==min(obob(N_it~=0)));
find(obob==min(obob(N_it~=0)));
end

n_it=0;
%robust=1;

for number=gg(1):gg(1)
G=Random_G(:,:,number);
U_G=Random(:,:,number);
Y_bar=inv(G'*G)*G'*z*U_G;
ob_pre=norm(z*U_G-G*Y_bar)^2;
diff=ob_pre;
n_it=0;
while ((diff))>0
for j=1:n
for m=1:n_groups
    count=zeros(1,n_groups);
    count(m)=1;
    F(j,m)=norm(z(j,:)*U_G-count*Y_bar)^2; 
end

for m=1:n_groups
    if min(F(j,:))==F(j,m)
        G(j,m)=1;
    else
        G(j,m)=0;
    end
end

end

clear Member
for j=1:n
    for m=1:n_groups
    if G(j,m)==1
       Member(j)=m;
    end
    end
end

centroids=Y_bar;
fac_scores=(z)*U_G;
cen_true=G*centroids;

if robust==1
clear diff_scores t_hotelling
for m=1:n_groups
for j=1:n
    diff_scores(:,j)=fac_scores(j,:)'-cen_true(Member(j),:)';
    t_hotelling(j)=n*diff_scores(:,j)'*inv(cov(fac_scores))*diff_scores(:,j);
    t_test(j)=(n-r)/(r*(n-1))*t_hotelling(j);
end
end

IQR=quantile(t_test,0.75)-quantile(t_test,0.25);
thr_clust_max=quantile(t_test,1-alpha);

clear flag
boundmax=thr_clust_max;
for j=1:n
    if t_test(j)>boundmax
       flag(j)=1;
    else flag(j)=0;
    end
end

for j=1:n
    if flag(j)==1
       G(j,:)=zeros(1,n_groups);
    end
end
end

if rank(G)==size(G,2)
[U_G, D_G]=svds(z'*(G*inv(G'*G)*G'-eye(n))*z,i1);
Y_bar=inv(G'*G)*G'*z*U_G;
if n_it>1
ob_pre=ob_post;
end
ob_post= norm(z*U_G-G*Y_bar)^2;
diff=ob_post-ob_pre;
diff_pre=diff;
n_it=n_it+1;
if diff<0
YY(:,:,number)=Y_bar;
GG(:,:,number)=G;
UU_GG(:,:,number)=U_G;
obob(number)=ob_post;
diffdiff(number)=diff;
end
else diff=-1;
end
end
end

%%OUTPUTS

if robust==1
Member;
Member_mod=Member;
for j=1:n
    if flag(j)==1
       Member_mod(j)=0;
    end
end
Member_mod;

out_tfkm=find(Member_mod==0);
n_out_tfkm=length(out_tfkm);

clear group_num % gives number of members of each group, excluding outliers!
for m=1:n_groups
group_num(m)=sum(flag(Member_mod==m));
end
group_num;
end

% here the composition also including outliers
clear group_comp
for m=1:n_groups
group_comp(m)=length(find(Member==m));
end
sum(group_comp);
group_comp;

if robust==0
    out_tfkm=0;
    n_out_tfkm=[];
    group_num=group_comp;
    Member_mod=Member;
end

   r_yes=i1;
   g_yes=i2;

output_all={fac_scores,'fac_scores',Member_mod,'Member_mod',group_num,'group_num',out_tfkm,'out_tfkm',n_out_tfkm,'n_out_tfkm',Member,'Member',group_comp,'group_comp',Hartigan,'Hartigan',r_yes,'rank',g_yes,'cluster_number',model_select,'model_select'}
    