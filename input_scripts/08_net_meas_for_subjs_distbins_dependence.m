% %% SETUP
% %Running Locally
datadir=fullfile('~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists'
outdir='~/Documents/bassett_lab/tooleyEnviNetworks/analyses'
% 
% %Running Locally Bassett
% datadir=fullfile('~/Documents/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
% listdir='~/Documents/tooleyEnviNetworks/subjectLists'
% outdir='~/Documents/tooleyEnviNetworks/analyses'

%Running on the cluster
datadir=fullfile('/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='/data/jag/bassett-lab/tooleyEnviNetworks/subjectLists'
outdir='/data/jag/bassett-lab/tooleyEnviNetworks/analyses/null_models_subjects_dist_bins'

%for each subject
%read the subject list in without the header
subjlist=csvread(fullfile(listdir,'n1012_healthT1RestExclude_parcels.csv'),1, 0 )

%% Distance Dependence of the Effect Distance Bins
%load Euclidean distances matrix of the Glasser nodes
distmatdir='~/Documents/bassett_lab/tooleyEnviNetworks/parcels/Glasser';
distmat=load(fullfile(distmatdir,'GlasserEuclideanDistanceMatrix.txt'));

%take parcel 52 out of the distance matrix
distmat=removerows(distmat, 'ind', [103]);
    %remove column 103
distmat(:,103)=[];
%threshold the distance matrix at several distances
for p=0:9
    my_field=strcat('threshold', num2str(p), 'to', num2str(p+1), 'longest')
    matrices.(my_field)=threshold_proportional_bins(distmat,(0.1 * p), ((0.1*p)+0.10));
    %convert to binary for each threshold
    mean(matrices.(my_field)(matrices.(my_field)~=0)) %check average length for each bin
    matrices.(my_field)=weight_conversion(matrices.(my_field), 'binarize');
end
for n=1:length(subjlist)
     sub=subjlist(n,2)
    %load Pearson correlation matrix
    file=fullfile(datadir,strcat(num2str(sub),'_GlasserPNC_znetwork.txt'));
	subfcmat = load(file);
    %parcel 52 is already gone
    %replace the diagonal of 1's with 0's
    for x=1:359
        subfcmat(x,x)=0;
    end
    %multipy subfcmat by binarized distance mat
    mymat=matrices.threshold0to1longest.*subfcmat;
    %calculate clustco on the resulting weights
    avgclustco_both_0to1_longest(n,1)=mean(clustering_coef_wu_sign(mymat,3));
    avgweight_0to1_longest(n,1)=mean(mymat(mymat~=0));
    mymat=matrices.threshold1to2longest.*subfcmat;
    avgclustco_both_1to2_longest(n,1)=mean(clustering_coef_wu_sign(mymat,3));
    avgweight_1to2_longest(n,1)=mean(mymat(mymat~=0));
    mymat=matrices.threshold2to3longest.*subfcmat;
    avgclustco_both_2to3_longest(n,1)=mean(clustering_coef_wu_sign(mymat,3));
    avgweight_2to3_longest(n,1)=mean(mymat(mymat~=0));
    mymat=matrices.threshold3to4longest.*subfcmat;
    avgclustco_both_3to4_longest(n,1)=mean(clustering_coef_wu_sign(mymat,3));
    avgweight_3to4_longest(n,1)=mean(mymat(mymat~=0));
     mymat=matrices.threshold4to5longest.*subfcmat;
    avgclustco_both_4to5_longest(n,1)=mean(clustering_coef_wu_sign(mymat,3));
    avgweight_4to5_longest(n,1)=mean(mymat(mymat~=0));
     mymat=matrices.threshold5to6longest.*subfcmat;
    avgclustco_both_5to6_longest(n,1)=mean(clustering_coef_wu_sign(mymat,3));
    avgweight_5to6_longest(n,1)=mean(mymat(mymat~=0));
     mymat=matrices.threshold6to7longest.*subfcmat;
    avgclustco_both_6to7_longest(n,1)=mean(clustering_coef_wu_sign(mymat,3));
    avgweight_6to7_longest(n,1)=mean(mymat(mymat~=0));
     mymat=matrices.threshold7to8longest.*subfcmat;
    avgclustco_both_7to8_longest(n,1)=mean(clustering_coef_wu_sign(mymat,3));
    avgweight_7to8_longest(n,1)=mean(mymat(mymat~=0));
    mymat=matrices.threshold8to9longest.*subfcmat;
    avgclustco_both_8to9_longest(n,1)=mean(clustering_coef_wu_sign(mymat,3));
    avgweight_8to9_longest(n,1)=mean(mymat(mymat~=0));
    mymat=matrices.threshold9to10longest.*subfcmat;
    avgclustco_both_9to10_longest(n,1)=mean(clustering_coef_wu_sign(mymat,3));
    avgweight_9to10_longest(n,1)=mean(mymat(mymat~=0));
end
    %% Write outfiles
outfile=dataset(subjlist,avgclustco_both_0to1_longest, avgweight_0to1_longest, avgclustco_both_1to2_longest, avgweight_1to2_longest, avgclustco_both_2to3_longest, avgweight_2to3_longest, avgclustco_both_3to4_longest, avgweight_3to4_longest, avgclustco_both_4to5_longest, avgweight_4to5_longest, avgclustco_both_5to6_longest, avgweight_5to6_longest, avgclustco_both_6to7_longest, avgweight_6to7_longest, avgclustco_both_7to8_longest, avgweight_7to8_longest, avgclustco_both_8to9_longest, avgweight_8to9_longest, avgclustco_both_9to10_longest, avgweight_9to10_longest)
export(outfile,'File',fullfile(outdir,'n1012_sub_net_meas_signed_distance_bins.csv'),'Delimiter',',')

%% Edge weight Dependence of the Effect
for n=1:length(subjlist)
    for p=0:9
    sub=subjlist(n,2);
    %load Pearson correlation matrix
    file=fullfile(datadir,strcat(num2str(sub),'_GlasserPNC_znetwork.txt'));
	subfcmat = load(file);
    %parcel 52 is already gone
    %replace the diagonal of 1's with 0's
    for x=1:359
        subfcmat(x,x)=0;
    end
    %threshold for bins 1-10% edge weight, 20-30%, etc.
    %wrote in function for binning thresholds
    threshmat=threshold_proportional_bins(subfcmat, (0.1 * p), ((0.1*p)+0.10));
    %check mean
    %mean(threshmat(threshmat~=0))
    my_field=strcat('threshold', num2str(p), 'to', num2str(p+1));
    %calculate clustco on thresholded matrix
    avgweight.(my_field)(n,1)=mean(threshmat(threshmat~=0));
    avgclustco_both.(my_field)(n,1)=mean(clustering_coef_wu_sign(threshmat,3));
  
    end
end
%%write outfile
outfile=dataset(subjlist, avgweight.threshold0to1, avgclustco_both.threshold0to1, avgweight.threshold1to2, avgclustco_both.threshold1to2 ,avgweight.threshold2to3, avgclustco_both.threshold2to3, avgweight.threshold3to4,avgclustco_both.threshold3to4,avgweight.threshold4to5, avgclustco_both.threshold4to5, avgweight.threshold5to6, avgclustco_both.threshold5to6, avgweight.threshold6to7, avgclustco_both.threshold6to7, avgweight.threshold7to8, avgclustco_both.threshold7to8, avgweight.threshold8to9, avgclustco_both.threshold8to9, avgweight.threshold9to10, avgclustco_both.threshold9to10)
export(outfile,'File',fullfile(outdir,'n1012_sub_net_meas_signed_weight_bins.csv'),'Delimiter',',')
