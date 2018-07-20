%% SETUP
%Running Locally
%Use NON-Z-SCORED FC MATRICES
datadir=fullfile('~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists'
outdir='/Users/utooley/Dropbox (Personal)/bassett_lab/clustco_paper/'
parceldir='~/Documents/bassett_lab/tooleyEnviNetworks/parcels/Glasser/'
clustcodir='~/Dropbox/bassett_lab/clustco_paper/'

%Running Locally Bassett
datadir=fullfile('~/Documents/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCNetworks/')
listdir='~/Documents/tooleyEnviNetworks/subjectLists'
clustcodir='~/Dropbox/bassett_lab/clustco_paper/'
parceldir='~/Documents/tooleyEnviNetworks/parcels/Glasser/'
outdir='~/Dropbox/bassett_lab/clustco_paper/brains'
subinfodir='~/Documents/tooleyEnviNetworks/data/subjectData/'
qadir='~/Documents/tooleyEnviNetworks/data/rest/'

%Running on the cluster
datadir=fullfile('/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='/data/jag/bassett-lab/tooleyEnviNetworks/subjectLists'
outdir='/data/jag/bassett-lab/tooleyEnviNetworks/analyses'
%read the subject list in without the header
subjlist=csvread(fullfile(listdir,'n1012_healthT1RestExclude_parcels.csv'),1, 0 )

stackedMatrix = zeros(359, 359, 1012);
for n=1:length(subjlist);
    sub=subjlist(n,2);
    %load Pearson correlation matrix
    file=fullfile(datadir,strcat(num2str(sub),'_GlasserPNC_znetwork.txt'));
	subfcmat = load(file);
    %parcel 52 is already gone
    %replace the diagonal of 1's with 0's
    %for x=1:359
    %    subfcmat(x,x)=0;
    %end
    
    stackedMatrix(:,:, n)=subfcmat;
    
end

meanMatrix = mean(stackedMatrix,3); %doc mean for more info.
%make the diagonal 0's
    for x=1:359
        meanMatrix(x,x)=0;
    end
%remove missing ROI of 52
meanMatrix=removerows(meanMatrix, 'ind', [103]);
    %remove column 103
meanMatrix(:,103)=[];
%plot it
imagesc(meanMatrix);
colormap(jet);
colorbar;
%export illustrator compatible image

%export mean matrix
csvwrite(fullfile(outdir,'average_FC_mat_all_subjects.csv'), meanMatrix)

%% rewrite edges for each subject into a vector
size_vec=triu(ones(359,359),1)
for n=1:length(subjlist);
    sub=subjlist(n,2);
    %load Pearson correlation matrix
    file=fullfile(datadir,strcat(num2str(sub),'_GlasserPNC_znetwork.txt'));
	subfcmat = load(file);
    %cycle through edges from 1,1 to 1,2 to 2,2 to 1,3 to 2,3 tec. (vertically)
    %rewrite into a vector
    vector_edges(n,:)=subfcmat(size_vec==1);
end
outfile=dataset(subjlist, vector_edges)
export(outfile,'File',fullfile(outdir,'zedges_for_each_subj_64261.csv'),'Delimiter',',')

%export file
csvwrite(fullfile(outdir,'zedges_for_each_subj_64621.csv'), vector_edges)

%% write edge weight betas/etc back into a matrix
size_vec=triu(ones(359,359),1)
beta_edge_weight_mat=size_vec
beta_edge_weight(size_vec ==1)= vector_edges



n = numel(vector_edges)/2-1; % compute size of output matrix 
B = zeros(round(n)); % define B of required size
[ii jj] = ndgrid(1:n); % ii and jj are row and column indices respectively
B(ii>=jj) = A; % fill in values in column-major order
B = B.'; % transpose to get result

B = triu(ones(vector_edges));
B(B==1) = vector_edges;
triu(toeplitz(vector_edges)
  

%% do the regression for the age x SES effect on each edge
%take out the extra parcel if using the 360 x 360 matrix
stackedMatrix=removerows(stackedMatrix, 'ind', [103]);
    %remove column 103
stackedMatrix(:,103)=[];
%read in the covariates data
subjlist=readtable(fullfile(listdir,'n1012_healthT1RestExclude_parcels.csv'))
demo=readtable(fullfile(subinfodir,'n1601_demographics_go1_20161212.csv'))
enviro=readtable(fullfile(subinfodir, 'n1601_go1_environment_factor_scores_tymoore_20150909.csv'))
%Get QA values to include in analyses
restqa=readtable(fullfile(qadir, 'n1601_RestQAData_20170509.csv'));
master=join(subjlist,demo, 'Keys', 'scanid');
master=join(master,enviro, 'Keys', 'scanid');
master=join(master,restqa, 'Keys', 'scanid');
%set up covariates data
%make sure factor variables are factors and numeric are numeric
master.race2=categorical(master.race2, [1 2 3], {'White', 'Black', 'Other'});
master.sex=categorical(master.sex, [1 2], {'Male', 'Female'});
%cell2table
master.restRelMeanRMSMotion=str2double(master.restRelMeanRMSMotion);
%should center variables before looking at interactions
master.ageAtScan1cent=(master.ageAtScan1-mean(master.ageAtScan1));
master.ageAtScan1yrs=(master.ageAtScan1)/12;
master.ageAtScan1yrscent=(master.ageAtScan1yrs-mean(master.ageAtScan1yrs));
master.envSEScent=(master.envSES-mean(master.envSES))
%split SES on the median 
%master.envSEShigh=discretize(master.envSES, [-2.6 0.0178 1.5],'categorical', {'Low', 'High'});
master.envSEShigh=discretize(master.envSES, [-2.6 0.0178 1.5]);

%cycle through each edge and regress the weight on predictors, get the
%pvals for most significant edges in age x SES interaction
" ~ ageAtScan1cent+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh+ageAtScan1cent*envSEShigh"
predictors=[ones(size(master.ageAtScan1cent)) master.ageAtScan1cent master.sex master.race2 master.restRelMeanRMSMotion master.envSEShigh master.ageAtScan1cent.*master.envSEShigh]
lm=fitlm(master, 'restRelMeanRMSMotion~ageAtScan1cent+sex+race2+envSEShigh+ageAtScan1cent.*envSEShigh')


for edges=1:359
      Iset=nonzeros(triu( reshape(1:numel(subfcmat), size(subfcmat)) ));
      for i =Iset
          
      end
      
pvals=regress(stackedMatrix(edges, edges,:), predictors)

end
