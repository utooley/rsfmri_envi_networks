%% SETUP
%Running Locally
datadir=fullfile('~/Documents/projects/in_progress/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='~/Documents/projects/in_progress/tooleyEnviNetworks/subjectLists'
outdir='~/Documents/projects/in_progress/tooleyEnviNetworks/analyses'
parceldir='~/Documents/projects/in_progress/tooleyEnviNetworks/parcels'

%Running Locally Bassett
datadir=fullfile('~/Documents/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='~/Documents/tooleyEnviNetworks/subjectLists'
outdir='~/Documents/tooleyEnviNetworks/analyses'

%Running on the cluster
datadir=fullfile('/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='/data/jag/bassett-lab/tooleyEnviNetworks/subjectLists'
outdir='/data/jag/bassett-lab/tooleyEnviNetworks/analyses'
%Load each Pearson FC matrix
%ursula(:,:,1)=load(fullfile(datadir,'2632_GlasserPNC_network.txt'))

%read the subject list in without the header
subjlist=csvread(fullfile(listdir,'n1012_healthT1RestExclude_parcels.csv'),1, 0 )
%read in the mapping of Glasser to Yeo networks
yeo_mapping=csvread(fullfile(parceldir, 'Glasser_to_Yeo.csv'),1,1 )
%can't insert a row because it will throw things off with 0's.
%remove the extra row from Yeo, it's 103
yeo_mapping(103,:)=[]

%for each subject
for n=1:length(subjlist)
    sub=subjlist(n,2);
    %load Pearson correlation matrix
    file=fullfile(datadir,strcat(num2str(sub),'_GlasserPNC_znetwork.txt'));
	subfcmat = load(file);
    %parcel 52 is already gone
    %replace the diagonal of 1's with 0's
    for x=1:359
        subfcmat(x,x)=0;
    end
    
    %%% Network Measures
% %No thresholding, just use the z-transformed pearson correlations....

%average network strength (the mean of all network weights in the matrix that are not equal to
%0))
avgweight(n,1)=mean(subfcmat(subfcmat~=0));

%CLUSTERING COEFFICIENT
%Documentation for BCT says must normalize prior to using the weighted
%clustering coefficient, per Mikhail it will not change things when using
%Constantini formula

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

%Network/system segregation per Chan et al. 2018    
%Calculate system segregation
[S W B]=segregation(subfcmat, yeo_mapping);
sys_segreg(n,1)=S;

%Threshold networks for positive-only weights, as Chan 2018 did, and then
%calculate system segregation
subfcmatthresh=threshold_absolute(subfcmat, 0);
[S W B]=segregation(subfcmatthresh, yeo_mapping);
sys_segreg_posonly(n,1)=S;
end

%% Write outfiles
outfile=dataset(subjlist, avgweight, avgclustco_both, modul)
export(outfile,'File',fullfile(outdir,'n1012_sub_net_meas_signed.csv'),'Delimiter',',')

%% Write system segregation outfile
outfile=dataset(subjlist, avgweight, avgclustco_both, modul, sys_segreg, sys_segreg_posonly)
export(outfile,'File',fullfile(outdir,'n1012_sub_net_meas_signed_w_segregation.csv'),'Delimiter',',')

%%just in case
%csvwrite(fullfile(outdir, 'avgclustc.csv'),avgclustc)
%csvwrite(fullfile(outdir, 'mod.csv'),modul)

