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

%% Compare with network null models?

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
%take a given subject's network, randomize it, 
%this function DOES NOT preserve the strength function in weighted networks
%should change the number of iterations?
for i=1:100
    null_1(:,:,i)=randmio_und_signed(subfcmat, 20);
end
%this function does preserve both the strength and degree distribution
for i=1:100
    null_2(:,:,i)=null_model_und_sign(subfcmat, 5, 0.3);
end

%run for net=null_1 and then net=null_2 separately
net=null_2;
%calculate net meas on null net models
avgweight(n,1)=mean(net(net~=0));
%CLUSTERING COEFFICIENT
%Using Constantini & Perugini's generalization for both neg and pos weights
avgclustco_both(n,1)=mean(clustering_coef_wu_sign(net,3));
%Community Louvain outputs a measure of modularity and can take signed nets
%[M Q]=community_louvain(net, 1, [], 'negative_asym');
%modul(n,1)=Q;
end

%change this path for exporting null_1 or null_2
outfile=dataset(subjlist, avgweight, avgclustco_both)
export(outfile,'File',fullfile(outdir,'n1012_net_meas_null_2_normed3clustcoefs_signed.csv'),'Delimiter',',')
