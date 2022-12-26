function [H,Hi] = Main_NFCCE(MuxNetwork,NOC,alpha,factorization)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
X=MuxNetwork;%Multiplex network
k=NOC;%number of clusters
a=alpha;
Hi=cell(1,size(X,2));
Si=cell(1,size(X,2));

n1=size(X{1,1},1);
if strcmp(factorization, 'SNMF')
    for i=1:size(X,2)
     Hi{1,i}=Symnmfrule_Esraa(X{1,i},k);
    end
    A_avg=zeros(n1,n1);
    for i=1:size(X,2)
      Ai_mod=X{1,i}+(a/2)*(Hi{1,i}*Hi{1,i}');
      A_avg=A_avg+Ai_mod;
       clear Ai_mod     
    end
    H=Symnmfrule_Esraa(A_avg,k);
    
    
elseif strcmp(factorization, 'PNMF')
    for i=1:size(X,2)
     Hi{1,i}=pnmfeu_Esraa(X{1,i},k);
    end
    A_avg=zeros(n1,n1);
    
    for i=1:size(X,2)
      Ai_mod=(X{1,i}*X{1,i}')+a*(Hi{1,i}*Hi{1,i}');
      A_avg=A_avg+Ai_mod;
       clear Ai_mod     
    end
    H=pnmfeu_Esraa(A_avg,k);
    
    
elseif strcmp(factorization, 'SNMTF')
    for i=1:size(X,2)
    [~,Si{1,i},Hi{1,i},~,~,~]=orthnmfrule_Esraa(X{1,i},k);
    end
    H = SNMTF_Esraa(X,k,Si,Hi,a);
    
    
    
    
end

% [~,Labels]=max(H,[],2);



end

