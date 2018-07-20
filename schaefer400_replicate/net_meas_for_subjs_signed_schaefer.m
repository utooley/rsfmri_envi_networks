%% FUNCTION SETUP
%% SETUP
%Running Locally
datadir=fullfile('~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_schaefer400/Schaefer400zNetworks/')
listdir='~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists'
outdir='~/Documents/bassett_lab/tooleyEnviNetworks/analyses'

%Running Locally Bassett
%datadir=fullfile('~/Documents/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
%listdir='~/Documents/tooleyEnviNetworks/subjectLists'
%outdir='~/Documents/tooleyEnviNetworks/analyses'


%Running on the cluster
datadir=fullfile('/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/restNetwork_schaefer400/Schaefer400zNetworks/')
listdir='/data/jag/bassett-lab/tooleyEnviNetworks/subjectLists'
outdir='/data/jag/bassett-lab/tooleyEnviNetworks/analyses/'
%Load each Pearson FC matrix
%ursula(:,:,1)=load(fullfile(datadir,'2632_GlasserPNC_network.txt'))

%for each subject
%read the subject list in without the header
subjlist=csvread(fullfile(listdir,'n1015_healthT1RestExclude.csv'),1, 0 )
%preallocate variables
avgclustco_both=zeros(length(subjlist));
avgweight=zeros(length(subjlist),1);
modul=zeros(length(subjlist));

for n=1:length(subjlist)
    sub=subjlist(n,2)
    %load Pearson correlation matrix
    file=fullfile(datadir,strcat(num2str(sub),'_Schaefer400_znetwork.txt'));
	subfcmat = load(file);
    %parcel 52 is already gone
    %replace the diagonal of 1's with 0's
    for x=1:359
        subfcmat(x,x)=0;
    end

%% Network Measures
%No thresholding, just use the z-transformed pearson correlations....

%average network strength (the mean of all network weights in the matrix that are not equal to
%0))
avgweight(n,1)=mean(subfcmat(subfcmat~=0));

%Average strength 
%Local assortivity
% [loc_assort_pos,loc_assort_neg]=local_assortativity_wu_sign(subfcmat);
% avglocass_pos(n,1)=sum(loc_assort_pos);
% avglocass_neg(n,1)=sum(loc_assort_neg);

%CLUSTERING COEFFICIENT
%Using Onnela et al. default clust coeff, which is calculated separately
%for pos and neg weights
%Documentation for BCT says must normalize prior to using the weighted
%clustering coefficient.

%Normalize weights to the range [0,1]
%normsubfcmat=weight_conversion(subfcmat, 'normalize');
% [C_pos,C_neg,Ctot_pos, Ctot_neg]=clustering_coef_wu_sign(subfcmat);
% avgclco_pos(n,1)=Ctot_pos;
% avgclco_neg(n,1)=Ctot_neg;

%Using Constantini & Perugini's generalization of the Zhang &
%Horvath formula (option2). This formula takes both positive &
%negative weights into account simultaneously, produces 1 value for each
%node. Then take the mean.
avgclustco_both(n,1)=mean(clustering_coef_wu_sign(subfcmat,3));

%Community Louvain outputs a measure of modularity and can take signed
%networks as input. Weighted the negative connections asymmetrically, Q* as
%recommended by Rubinov & Sporns
[M Q]=community_louvain(subfcmat, 1, [], 'negative_asym');
modul(n,1)=Q;


%% PRIOR TO THE MEASURES BELOW MUST CONVERT TO A CONNECTION_LENGTHS matrix
%%L and then a distance matrix D
% connlenmat=weight_conversion(subfcmat, 'lengths');
% distmat=distance_wei(connlenmat);
% 
% %characteristic path length
%  %Characteristic path length is defined here as the mean shortest
% %   path length between all pairs of nodes, for consistency with common
% %   usage.
% cpathleng(n,1)=charpath(distmat);

%betweenness centrality
%Node betweenness centrality is the fraction of all shortest paths in 
%   the network that contain a given node. Nodes with high values of 
%   betweenness centrality participate in a large number of shortest paths.
%betweenc=betweenness_wei(distmat);

end

%% Write outfiles
outfile=dataset(subjlist, avgweight(:,1), avgclustco_both(:,1), modul(:,1))
export(outfile,'File',fullfile(outdir,'n1015_sub_net_meas_schaefer_signed.csv'),'Delimiter',',')
% %% Compare with network null models?
% 
% %for n=1:length(subjlist)
%     %sub=subjlist(n,2)
%     %load Pearson correlation matrix
%     file=fullfile(datadir,strcat(num2str(sub),'_GlasserPNC_znetwork.txt'));
% 	subfcmat = load(file);
%     %parcel 52 is already gone
%     %replace the diagonal of 1's with 0's
%     for x=1:359
%         subfcmat(x,x)=0;
%     end
% %make 100 null models and take the mean clustering coeff across them
% %take a given subject's network, randomize it, 
% %this function DOES NOT preserve the strength function in weighted networks
% %should change the number of iterations?
% null_1=zeros(359,359,100);
% avgweight=zeros(100,1);
% avgclustco_both=zeros(100,1);
% for i=1:100
%     null_1(:,:,i)=randmio_und_signed(subfcmat, 20);
%     avgweight(i,1)=mean(null_1(null_1~=0));
%     avgclustco_both(i,1)=mean(clustering_coef_wu_sign(null_1,3));
% end
% %add each subject's average across 100 null models for clust co and avg
% %weight
% avgweight_null1=mean(avgweight(:,1))
% avgclustco_both_null1=mean(avgclustco_both(:,1))
% 
% %make 100 null models and take the mean clustering coeff across them
% %this function does preserve both the strength and degree distribution
% null_2=zeros(359,359,100);
% avgweight=zeros(100,1);
% avgclustco_both=zeros(100,1);
% for i=1:100
%     null_2(:,:,i)=null_model_und_sign(subfcmat, 5, 0.3);
%     avgweight(i,1)=mean(null_2(null_2~=0));
%     avgclustco_both(i,1)=mean(clustering_coef_wu_sign(null_2,3));
% end
% %add each subject's average across 100 null models for clust co and avg
% %weight
% avgweight_null2=mean(avgweight(:,1))
% avgclustco_both_null2=mean(avgclustco_both(:,1))
% 
% 
% % %run for net=null_1 and then net=null_2
% % net=null_1;
% % %calculate net meas on null net models
% % avgweight(n,1)=mean(net(net~=0));
% % %Local assortivity
% % [loc_assort_pos,loc_assort_neg]=local_assortativity_wu_sign(net);
% % avglocass_pos(n,1)=sum(loc_assort_pos);
% % avglocass_neg(n,1)=sum(loc_assort_neg);
% % %CLUSTERING COEFFICIENT
% % %Using Onnela et al. default clust coeff, sep for pos and neg weights
% % [C_pos,C_neg,Ctot_pos, Ctot_neg]=clustering_coef_wu_sign(net);
% % avgclco_pos(n,1)=Ctot_pos;
% % avgclco_neg(n,1)=Ctot_neg;
% % %Using Constantini & Perugini's generalization for both neg and pos weights
% % avgclustco_both(n,1)=mean(clustering_coef_wu_sign(net,3));
% % %Community Louvain outputs a measure of modularity and can take signed nets
% % [M Q]=community_louvain(net, 1, [], 'negative_asym');
% % modul(n,1)=Q;
% % %PRIOR TO THE MEASURES BELOW MUST CONVERT TO A CONNECTION_LENGTHS matrix
% % %%L and then a distance matrix D
% % connlenmat=weight_conversion(net, 'lengths');
% % distmat=distance_wei(connlenmat);
% % %characteristic path length
% % cpathleng(n,1)=charpath(distmat);
% %end
% 
% %change this path for exporting null_1 or null_2
% outfile=padcat(avgweight_null1, avgweight_null2, avgclustco_both_null1, avgclustco_both_null2)
% dlmwrite(fullfile(outdir,strcat(num2str(sub), '_null_models_100x.txt')), outfile)

% %change this path for exporting null_1 or null_2
% outfile=dataset(subjlist, avgweight, avglocass_pos, avglocass_neg, avgclco_pos, avgclco_neg, avgclustco_both, modul, cpathleng)
% export(outfile,'File',fullfile(outdir,'n1012_net_meas_null_1_20rewir_signed.csv'),'Delimiter',',')


%%just in case
%csvwrite(fullfile(outdir, 'assort.csv'),assort)
%csvwrite(fullfile(outdir, 'avgclustc.csv'),avgclustc)
%csvwrite(fullfile(outdir, 'Eglob.csv'),Eglob)
%csvwrite(fullfile(outdir, 'avgeloc.csv'),avgeloc)
%csvwrite(fullfile(outdir, 'mod.csv'),modul)
%csvwrite(fullfile(outdir, 'cpathleng.csv'),cpathleng)

