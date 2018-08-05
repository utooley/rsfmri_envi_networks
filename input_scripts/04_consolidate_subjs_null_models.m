%% SETUP
%Running Locally
datadir=fullfile('~/Documents/bassett_lab/tooleyEnviNetworks/analyses/null_models_subjects/')
listdir='~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists'
outdir='~/Documents/bassett_lab/tooleyEnviNetworks/analyses/null_models_subjects/'

subjlist=csvread(fullfile(listdir,'n1012_healthT1RestExclude_parcels.csv'),1, 0 )
null_models_sub=zeros(1012,4)
for n=1:length(subjlist)
    sub=subjlist(n,2);
    %load null' models dist bins
    file=fullfile(datadir,strcat(num2str(sub),'_null_models_100x.txt'));
	null_model_coef=readtable(file,'ReadVariableNames', 0);
    null_model_coef=table2array(null_model_coef)';
    null_models_sub(n,:)=null_model_coef;
    
end
header={'avgweight_null1', 'avgweight_null2', 'avgclustco_both_null1', 'avgclustco_both_null2'}
outfile=[header; num2cell(null_models_sub)]
l=dataset(subjlist,null_models_sub)
export(l,'File', fullfile(outdir,'n1012_sub_null_models.csv'), 'Delimiter', ',')