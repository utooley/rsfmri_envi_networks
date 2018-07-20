%Running Locally Personal
datadir=fullfile('~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_schaefer400/Schaefer400Networks')
listdir='~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists'
outdir='~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_schaefer400/Schaefer400zNetworks'

%Running Locally Bassett
datadir=fullfile('~/Documents/tooleyEnviNetworks/data/rest/restNetwork_schaefer400/Schaefer400Networks')
listdir='~/Documents/tooleyEnviNetworks/subjectLists'
outdir='~/Documents/tooleyEnviNetworks/data/rest/restNetwork_schaefer400/Schaeferz400Networks'

%Running on the cluster
datadir=fullfile('/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/restNetwork_schaefer400/Schaefer400Networks')
listdir='/data/jag/bassett-lab/tooleyEnviNetworks/subjectLists'
outdir='/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/restNetwork_schaefer400/Schaefer400zNetworks'

%get the subject list, including everyone, as no one has NAs
subjlist=csvread(fullfile(listdir,'n1015_healthT1RestExclude.csv'),1, 0 )

for n=1:length(subjlist)
    sub=subjlist(n,2)
    file=fullfile(datadir,strcat(num2str(sub),'_Schaefer400_network.txt'))
	subfcmat = load(file);
    %replace the diagonal of 1's with 0's
    for x=1:400
        subfcmat(x,x)=0;
    end
    %create an empty z-matrx
    zfc=[]
    for i=1:400
        %cycle through each column of the FC matrix and do a fisher r-to-z
        %for each value
        zfc(:,i)=fisherz(subfcmat(:,i));
    end
    outfile=fullfile(outdir, strcat(num2str(sub),'_Schaefer400_znetwork.txt'))
    csvwrite(outfile, zfc)
end