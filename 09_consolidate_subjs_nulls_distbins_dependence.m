%% SETUP
%Running Locally
datadir=fullfile('~/Documents/bassett_lab/tooleyEnviNetworks/analyses/null_models_subjects_dist_bins/')
listdir='~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists'
outdir='~/Documents/bassett_lab/tooleyEnviNetworks/analyses/null_models_subjects_dist_bins/'

subjlist=csvread(fullfile(listdir,'n1012_healthT1RestExclude_parcels.csv'),1, 0 )
null_models_sub=zeros(1012,20)
for n=1:length(subjlist)
    sub=subjlist(n,2);
    %load null' models dist bins
    file=fullfile(datadir,strcat(num2str(sub),'_null_models_dist_bins.txt'));
	null_model_coef=readtable(file);
    null_model_coef=table2array(null_model_coef)';
    null_models_sub(n,:)=null_model_coef;
    
end
header={'avgclustco_both_0to1_longestnull', 'avgweight_0to1_longestnull', 'avgclustco_both_1to2_longestnull', 'avgweight_1to2_longestnull', 'avgclustco_both_2to3_longestnull', 'avgweight_2to3_longestnull', 'avgclustco_both_3to4_longestnull', 'avgweight_3to4_longestnull', 'avgclustco_both_4to5_longestnull', 'avgweight_4to5_longestnull', 'avgclustco_both_5to6_longestnull', 'avgweight_5to6_longestnull', 'avgclustco_both_6to7_longestnull', 'avgweight_6to7_longestnull', 'avgclustco_both_7to8_longestnull', 'avgweight_7to8_longestnull', 'avgclustco_both_8to9_longestnull', 'avgweight_8to9_longestnull', 'avgclustco_both_9to10_longestnull', 'avgweight_9to10_longestnull'}
outfile=[header; num2cell(null_models_sub)]   
l=dataset(outfile)
export(l,'File', fullfile(outdir,'n1012_sub_null_models_dist_bins.csv'), 'Delimiter', ',')