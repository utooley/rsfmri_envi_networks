%% FUNCTION SETUP
function [avgclustco_both_0to1_longestnull, avgweight_0to1_longestnull, avgclustco_both_1to2_longestnull, avgweight_1to2_longestnull, avgclustco_both_2to3_longestnull, avgweight_2to3_longestnull, avgclustco_both_3to4_longestnull, avgweight_3to4_longestnull, avgclustco_both_4to5_longestnull, avgweight_4to5_longestnull, avgclustco_both_5to6_longestnull, avgweight_5to6_longestnull, avgclustco_both_6to7_longestnull, avgweight_6to7_longestnull, avgclustco_both_7to8_longestnull, avgweight_7to8_longestnull, avgclustco_both_8to9_longestnull, avgweight_8to9_longestnull, avgclustco_both_9to10_longestnull, avgweight_9to10_longestnull] = net_meas_for_subjs_signed_dist_bins(sub)
% %% SETUP
% %Running Locally
% datadir=fullfile('~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
% listdir='~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists'
% outdir='~/Documents/bassett_lab/tooleyEnviNetworks/analyses'
% 
% %Running Locally Bassett
% datadir=fullfile('~/Documents/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
% listdir='~/Documents/tooleyEnviNetworks/subjectLists'
% outdir='~/Documents/tooleyEnviNetworks/analyses'

%Running on the cluster
datadir=fullfile('/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='/data/jag/bassett-lab/tooleyEnviNetworks/subjectLists'
outdir='/data/jag/bassett-lab/tooleyEnviNetworks/analyses/null_models_subjects_dist_bins'

%% Null Models for Distance Dependence of the Effect
% %load Euclidean distances matrix of the Glasser nodes
% distmatdir='~/Documents/tooleyEnviNetworks/parcels/Glasser';
distmatdir='~/Documents/bassett_lab/tooleyEnviNetworks/parcels/Glasser';
distmatdir='/data/jag/bassett-lab/tooleyEnviNetworks/parcels/Glasser/distances';
distmat=load(fullfile(distmatdir,'GlasserEuclideanDistanceMatrix.txt'));

%subjlist=csvread(fullfile(listdir,'n1012_healthT1RestExclude_parcels.csv'),1, 0 )
%take parcel 52 out of the distance matrix
distmat=removerows(distmat, 'ind', [103]);
    %remove column 103
distmat(:,103)=[];
%threshold the distance matrix at several distances
for p=0:9
    my_field=strcat('threshold', num2str(p), 'to', num2str(p+1), 'longest');
    matrices.(my_field)=threshold_proportional_bins(distmat,(0.1 * p), ((0.1*p)+0.10))
    %convert to binary for each threshold
    matrices.(my_field)=weight_conversion(matrices.(my_field), 'binarize')
end
%preallocate variables helps with parfor loop?
% avgclustco_both_0to1_longestnull_subs=zeros(1012,1)
% avgweight_0to1_longestnull_subs=zeros(1012,1)
% avgclustco_both_1to2_longestnull_subs=zeros(1012,1)
% avgweight_1to2_longestnull_subs=zeros(1012,1)
% avgclustco_both_2to3_longestnull_subs=zeros(1012,1)
% avgweight_2to3_longestnull_subs=zeros(1012,1)
% avgclustco_both_3to4_longestnull_subs=zeros(1012,1)
% avgweight_3to4_longestnull_subs=zeros(1012,1)
% avgclustco_both_4to5_longestnull_subs=zeros(1012,1)
% avgweight_4to5_longestnull_subs=zeros(1012,1)
% avgclustco_both_5to6_longestnull_subs=zeros(1012,1)
% avgweight_5to6_longestnull_subs=zeros(1012,1)
% avgclustco_both_6to7_longestnull_subs=zeros(1012,1)
% avgweight_6to7_longestnull_subs=zeros(1012,1)
% avgclustco_both_7to8_longestnull_subs=zeros(1012,1)
% avgweight_7to8_longestnull_subs=zeros(1012,1)
% avgclustco_both_8to9_longestnull_subs=zeros(1012,1)
% avgweight_8to9_longestnull_subs=zeros(1012,1)
% avgclustco_both_9to10_longestnull_subs=zeros(1012,1)
% avgweight_9to10_longestnull_subs=zeros(1012,1)
%for n=1:length(subjlist)
%datadir=fullfile('~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
file=fullfile(datadir,strcat(num2str(sub),'_GlasserPNC_znetwork.txt'));
subfcmat=load(file);
%parcel52isalreadygone
%replacethediagonalof1'swith0's
for x=1:359
subfcmat(x,x)=0;
end
%multipysubfcmatbybinarizeddistancemat
mymat=matrices.threshold0to1longest.*subfcmat;
mymat=randmio_und_signed(mymat,20);
%calculateclustcoontheresultingweights
avgclustco_both_0to1_longestnull=mean(clustering_coef_wu_sign(mymat,3));
avgweight_0to1_longestnull=mean(mymat(mymat~=0));
mymat=matrices.threshold1to2longest.*subfcmat;
mymat=randmio_und_signed(mymat,20);
avgclustco_both_1to2_longestnull=mean(clustering_coef_wu_sign(mymat,3));
avgweight_1to2_longestnull=mean(mymat(mymat~=0));
mymat=matrices.threshold2to3longest.*subfcmat;
mymat=randmio_und_signed(mymat,20);
avgclustco_both_2to3_longestnull=mean(clustering_coef_wu_sign(mymat,3));
avgweight_2to3_longestnull=mean(mymat(mymat~=0));
mymat=matrices.threshold3to4longest.*subfcmat;
mymat=randmio_und_signed(mymat,20);
avgclustco_both_3to4_longestnull=mean(clustering_coef_wu_sign(mymat,3));
avgweight_3to4_longestnull=mean(mymat(mymat~=0));
mymat=matrices.threshold4to5longest.*subfcmat;
mymat=randmio_und_signed(mymat,20);
avgclustco_both_4to5_longestnull=mean(clustering_coef_wu_sign(mymat,3));
avgweight_4to5_longestnull=mean(mymat(mymat~=0));
mymat=matrices.threshold5to6longest.*subfcmat;
mymat=randmio_und_signed(mymat,20);
avgclustco_both_5to6_longestnull=mean(clustering_coef_wu_sign(mymat,3));
avgweight_5to6_longestnull=mean(mymat(mymat~=0));
mymat=matrices.threshold6to7longest.*subfcmat;
mymat=randmio_und_signed(mymat,20);
avgclustco_both_6to7_longestnull=mean(clustering_coef_wu_sign(mymat,3));
avgweight_6to7_longestnull=mean(mymat(mymat~=0));
mymat=matrices.threshold7to8longest.*subfcmat;
mymat=randmio_und_signed(mymat,20);
avgclustco_both_7to8_longestnull=mean(clustering_coef_wu_sign(mymat,3));
avgweight_7to8_longestnull=mean(mymat(mymat~=0));
mymat=matrices.threshold8to9longest.*subfcmat;
mymat=randmio_und_signed(mymat,20);
avgclustco_both_8to9_longestnull=mean(clustering_coef_wu_sign(mymat,3));
avgweight_8to9_longestnull=mean(mymat(mymat~=0));
mymat=matrices.threshold9to10longest.*subfcmat;
mymat=randmio_und_signed(mymat,20);
avgclustco_both_9to10_longestnull=mean(clustering_coef_wu_sign(mymat,3));
avgweight_9to10_longestnull=mean(mymat(mymat~=0));

outfile=padcat(avgclustco_both_0to1_longestnull, avgweight_0to1_longestnull, avgclustco_both_1to2_longestnull, avgweight_1to2_longestnull, avgclustco_both_2to3_longestnull, avgweight_2to3_longestnull, avgclustco_both_3to4_longestnull, avgweight_3to4_longestnull, avgclustco_both_4to5_longestnull, avgweight_4to5_longestnull, avgclustco_both_5to6_longestnull, avgweight_5to6_longestnull, avgclustco_both_6to7_longestnull, avgweight_6to7_longestnull, avgclustco_both_7to8_longestnull, avgweight_7to8_longestnull, avgclustco_both_8to9_longestnull, avgweight_8to9_longestnull, avgclustco_both_9to10_longestnull, avgweight_9to10_longestnull)
dlmwrite(fullfile(outdir,strcat(num2str(sub), '_null_models_dist_bins.txt')), outfile)
%export(outfile,'File',fullfile(outdir,strcat(num2str(sub), '_null_models_dist_bins.txt')),'Delimiter',' ')
end