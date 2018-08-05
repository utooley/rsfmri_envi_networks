%% SETUP
%Running Locally
datadir=fullfile('~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists'
outdir='~/Documents/bassett_lab/tooleyEnviNetworks/analyses'

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

%for each subject
%read the subject list in without the header
subjlist=csvread(fullfile(listdir,'n1012_healthT1RestExclude_parcels.csv'),1, 0 )

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

end

%% Write outfiles
outfile=dataset(subjlist, avgweight, avgclustco_both, modul)
export(outfile,'File',fullfile(outdir,'n1012_sub_net_meas_signed.csv'),'Delimiter',',')

%%just in case
%csvwrite(fullfile(outdir, 'avgclustc.csv'),avgclustc)
%csvwrite(fullfile(outdir, 'mod.csv'),modul)

