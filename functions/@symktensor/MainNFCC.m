function [QualityMetrics]=MainNFCC(MuxNetwork,W,GT_CubeMux,MinNOC,MaxNOC,NOCBy,factorization)
%Input:
% MuxNetwork: is a structure variable of the multiplex network of L layers in the form
%            MuxNetwork.connections: Is a 1*L  cell variable that contains the adjacency matrices of the L layers 
%            MuxNetwork.GroundTruth: is a double variable that contains a nL*1 vector
%            of the nodes groundtruth

%W: Supra adjacency matrix
%GT_CubeMux: Ground truth membership matrix
%MinNOC: Minimum number of clusters 
%MaxNOC: maximum number of clusters
%NOCBy: 'dispersion' or 'modularity'
%factorization: 'SNMF' or 'PNMF' or 'SNMTF'
% W, GT_CubeMux,MinNOC,MaxNOC,NOCBy you need these if you are finding the
% number of clusters using modularity or dispersion coefficent 
% If the number of clusters is known you can set MinNOC=MaxNOC=the number
% of clusters and comment out the modularity line
X=MuxNetwork.connections;

Tic_NFCC=tic; 

if strcmp(NOCBy,'dispersion')

numSample=size(X{1,1},2);
rho=zeros(MaxNOC,1);

Labels_initial=cell(1,MaxNOC);
Y_initial=cell(1,MaxNOC);
for NOC=MinNOC:MaxNOC %round(sqrt(numSample))
    C=zeros(numSample,numSample);
    Y_i=zeros(numSample,NOC);
    
        [H,Hi]= Main_NFCCE(X,NOC,0.5,factorization);
      
        ind=getClusters(H');
        C=C+getRelationMatrix(ind);
        Y_i=H;
    
     rho(NOC)=dispersionCoefficient(C);
Y_initial{1,NOC}=Y_i;
[~,b]=max(getClusters(Y_i'));

Labels_initial{1,NOC}=b';

end
[val,k1]=max(rho);
k1=find(rho==val);
k1=k1(end);
Y_initial_suggested=Y_initial{1,k1};
Labels_n=Labels_initial{1,k1};
DNOC=k1;

elseif strcmp(NOCBy,'modularity')
    
  Labels_initial=cell(1,MaxNOC);
  Mod=zeros(1,MaxNOC);  
 for NOC=MinNOC:MaxNOC
      [H,Hi]= Main_NFCCE(X,NOC,0.5,factorization);
      [~,idx]=max(getClusters(H'));
      Labels_initial{1,NOC}=idx';
      Mod(1,NOC)= modularityMIT(repmat(Labels_initial{1,NOC},size(X,2),1),W);%uncomment
      %this
     
 end   
    
 [~,b]=max(Mod(2:end));%uncomment this
DNOC=b+1;%uncomment this
% DNOC=NOC;
Labels_n=Labels_initial{1,DNOC};   
    
    
end

Labels=repmat(Labels_n,size(MuxNetwork.GroundTruth,1)/size(Labels_n,1),1);

NFCC_cpuTimes=toc(Tic_NFCC); 

[VIn_NFCCE, MIn_NFCCE] = partition_distance(MuxNetwork.GroundTruth,Labels);
[CA_NFCCE, CR_NFCCE, CP_NFCCE, CF_NFCCE, ~] = computePerformance(Labels,GT_CubeMux);
purity_NFCCE = calc_purity(MuxNetwork.GroundTruth,Labels);

QualityMetrics.VIn=VIn_NFCCE;
QualityMetrics.MIn=MIn_NFCCE;
QualityMetrics.CA=CA_NFCCE;
QualityMetrics.CR=CR_NFCCE;
QualityMetrics.CP=CP_NFCCE;
QualityMetrics.CF=CF_NFCCE;
QualityMetrics.purity=purity_NFCCE;
QualityMetrics.DNOC=DNOC;
QualityMetrics.cpuTime=NFCC_cpuTimes;

end

