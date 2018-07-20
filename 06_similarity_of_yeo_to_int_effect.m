%% SETUP
parceldir='~/Documents/tooleyEnviNetworks/parcels'

%get the interaction effects for the effect of interest
node_est=dlmread('~/Dropbox/bassett_lab/clustco_paper/nodewise_betas_for_lms.csv', ',', 1,1)
estimates=node_est(:,2)
%get the glasser-yeo mapping
nodes_in_yeo=csvread(fullfile(parceldir, 'Glasser_to_Yeo.csv'),1,1 )
%can't insert a row because it will throw the estimate of the effect size
%off
%estimates=insertrows(estimates, 0, 102)
%remove the row from Yeo, it's 103
nodes_in_yeo(103,:)=[]
%look at average estimate of interaction effect in each system
master=[nodes_in_yeo estimates]
for i=1:7
m(i,1)=mean(master(master(:,1)==i,2))
end
%% COMPARING SIMILARITY OF PARTITIONS

%normalize the interaction effect estimates (between 0 and 1)
estimates=estimates-min(estimates)
estimates=estimates/max(estimates)

%zRand between yeo assignments and interaction effect
[realzrand SR SAR VI]= zrand(nodes_in_yeo,(round(estimates*(length(unique(nodes_in_yeo)) -1 ))));

%generate a distribution of zRand scores, by permuting the interaction
%estimates and calculating the similarity 10,000x.
for n=1:10000
    permestim=randsample(estimates, 359);
    [zr permSR permSAR VI]= zrand(nodes_in_yeo,(round(permestim*(length(unique(nodes_in_yeo))))));
    perm_zrand(n,1)=zr;
    perm_rand_sim_co(n,1)=permSR;
    perm_rand_adjust(n,1)=permSAR;
    perm_VI(n,1)=VI;
end

%make a histogram with with real z_rand value on it
hist(perm_zrand)
hold on
line([realzrand realzrand], [0 3500])
%make a histogram with with real rand similarity coefficient on it
hist(perm_rand_sim_co)
hold on
line([SR SR], [0 3500])
%make a histogram with with real adjusted rand similarity coefficient on it
hist(perm_rand_adjust)
hold on
line([SAR SAR], [0 3500])
%make a histogram with with real VI on it
hist(perm_VI)
hold on
line([VI VI], [0 3500])

%look at where the real zrand value falls in the permuted distribution
perm_zrand_dist=fitdist(perm_zrand,'Normal')
p=cdf(perm_zrand_dist, realzrand)
%look at where the real zrand adj value falls in the permuted distribution
perm_rand_adj_dist=fitdist(perm_rand_adjust,'Normal')
p=cdf(perm_rand_adj_dist, real_rand_adj)
%save the annotation contribution MI variables to read into R
MI_out=dataset(subjlist, similarity_diff, similarity_unpermuted, similarity_permuted)
export(MI_out,'File',fullfile(outdir,'wsbm_ct_annot_Pearson_k7_MI_annotation_contribution.csv'),'Delimiter',',')
