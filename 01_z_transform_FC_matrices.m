%Running Locally Personal
%datadir=fullfile('~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCTimeseries/')
%listdir='~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists'
%outdir='~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCWaveletNets/'

%Running Locally Bassett
datadir=fullfile('~/Documents/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCNetworks/')
listdir='~/Documents/tooleyEnviNetworks/subjectLists'
outdir='~/Documents/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/'

%Running on the cluster
datadir=fullfile('/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCTimeseries/')
listdir='/data/jag/bassett-lab/tooleyEnviNetworks/subjectLists'
outdir='/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/'

%get the subject list,excluding those who have NAs
subjlist=csvread(fullfile(listdir,'n1012_healthT1RestExclude_parcels.csv'),1, 0 )

for n=1:length(subjlist)
    sub=subjlist(n,2)
    file=fullfile(datadir,strcat(num2str(sub),'_GlasserPNC_network.txt'))
	subfcmat = load(file);
    %elimate parcel 52 (parcel index 103), delete row 103
    subfcmat=removerows(subfcmat, 'ind', [103]);
    %remove column 103
    subfcmat(:,103)=[];
    %replace the diagonal of 1's with 0's
    for x=1:359
        subfcmat(x,x)=0;
    end
    %create an empty z-matrx
    zfc=[]
    for i=1:359
        %cycle through each column of the FC matrix and do a fisher r-to-z
        %for each value
        zfc(:,i)=fisherz(subfcmat(:,i));
    end
    outfile=fullfile(outdir, strcat(num2str(sub),'_GlasserPNC_znetwork.txt'))
    csvwrite(outfile, zfc)
end