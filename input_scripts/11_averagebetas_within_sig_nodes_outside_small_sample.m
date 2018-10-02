%% SETUP
%Running Locally Personal Computer
datadir=fullfile('~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists'
outdir='~/Documents/bassett_lab/tooleyEnviNetworks/analyses'
clustcodir='~/Dropbox (Personal)/bassett_lab/clustco_paper'

%Running Locally Bassett
datadir=fullfile('~/Documents/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='~/Documents/tooleyEnviNetworks/subjectLists'
outdir='~/Documents/tooleyEnviNetworks/analyses'
clustcodir='~/Dropbox (Personal)/bassett_lab/clustco_paper'
  
%Running on the cluster
datadir=fullfile('/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='/data/jag/bassett-lab/tooleyEnviNetworks/subjectLists'
outdir='/data/jag/bassett-lab/tooleyEnviNetworks/analyses'
%Load each Pearson FC matrix
%ursula(:,:,1)=load(fullfile(datadir,'2632_GlasserPNC_network.txt'))
%what are the significant nodes?
indexofcols=[103 167 272 313]
indexofnodes=linspace(1,359,359)
indexofnonsignodes=setdiff( indexofnodes, indexofcols)
%read in the betas for edge weights for age x SES interaction
edgebetas=csvread(fullfile(clustcodir,'edge_betas_agexses_int_scaled_small_sample.csv'),1, 1 );
%edgebetas(:,1)=[]
%make a matrix with only the 26 nodes of significance and their edgess
%pull out only the 26 nodes that show an age x SES effect colnum-wise.
    for l=1:length(indexofcols);
    i=indexofcols(l);
    signodesonly(:,l)=edgebetas(:,i);
    end
%pull out only the 26 nodes that show an age x SES effect rowise.
    for j=1:length(indexofcols);
    i=indexofcols(j);
    signodesonly2(j,:)=signodesonly(i,:);
    end
    
%look at the avg value within those 26 node
avgbetainside=mean(signodesonly2(signodesonly2~=0))
%look at average value outside the 26 nodes
%make a matrix with only the 26 nodes of significance and their edges
a=edgebetas([indexofcols],:)
avgbetabetween=mean(a(a~=0))

%what is the true difference between the effect on edges within and
%outside?
truediff1=avgbetainside-avgbetabetween
%make a matrix of everything outside the 26 nodes of significance)
    a=edgebetas(:, [indexofcols indexofnonsignodes]);
    nodesandedges=a([indexofcols indexofnonsignodes],:);
    nonsigedgesonly=nodesandedges([5:359],[5:359]);
%look at the avg value within nonsig node edges only
avgbetatorest=mean(nonsigedgesonly(nonsigedgesonly~=0))

%what is the true difference between the effect on edges within and
%those in the rest of the brain?
truediff2=avgbetainside-avgbetatorest
%%  PERMUTATION TESTING

%%% WITHIN SIG NODES AND BETWEEN THOSE AND REST OF BRAIN

%randomise the edge betas at the beginning, then do everything the same,
%build up a distribution.
edgebetas=csvread(fullfile(clustcodir,'edge_betas_agexses_int_scaled_small_sample.csv'),1, 1 );
permdiff=[]
for x=1:10000
    %only shuffle the rows (as that's what I'm filtering on below), otherwise will disrupt the number of
    %0's in each row, leading to possibly different numbers of betas
    permedgebetas=edgebetas(randperm(size(edgebetas,1)),:);
%make a matrix with only the 26 nodes of significance and their edgess
%pull out only the 26 nodes that show an age x SES effect colnum-wise.
    for l=1:length(indexofcols);
    i=indexofcols(l);
    signodesonly(:,l)=permedgebetas(:,i);
    end
%pull out only the 26 nodes that show an age x SES effect rowise.
    for j=1:length(indexofcols);
    i=indexofcols(j);
    signodesonly2(j,:)=signodesonly(i,:);
    end
    
%look at the avg value within those 26 node
avgbetainside=mean(signodesonly2(signodesonly2~=0));
%look at average value outside the 26 nodes
%make a matrix with only the 26 nodes of significance and their edges
a=edgebetas([indexofcols],:);
avgbetabetween=mean(a(a~=0));

%what is the true difference between the effect on edges within and
%outside?
permdiff(x)=avgbetainside-avgbetabetween;
end

truediff= 
%make a histogram and see where the real difference falls on it
hist(permdiff)
hold on
line([truediff1 truediff1], [0 3500])
%get a p-value
perm_edgesbeta_dist=fitdist(permdiff','Normal')
p=cdf(perm_edgesbeta_dist, truediff1)

%%% WITHIN SIG NODES VS. NODES IN REST OF BRAIN

%randomise the edge betas at the beginning, then do everything the same,
%build up a distribution.
edgebetas=csvread(fullfile(clustcodir,'edge_betas_agexses_int_scaled.csv'),1, 1 );
permdiff=[]
for x=1:10000
    %only shuffle the rows (as that's what I'm filtering on below), otherwise will disrupt the number of
    %0's in each row, leading to possibly different numbers of betas
    permedgebetas=edgebetas(randperm(size(edgebetas,1)),:);
%make a matrix with only the 26 nodes of significance and their edgess
%pull out only the 26 nodes that show an age x SES effect colnum-wise.
    for l=1:length(indexofcols);
    i=indexofcols(l);
    signodesonly(:,l)=permedgebetas(:,i);
    end
%pull out only the 26 nodes that show an age x SES effect rowise.
    for j=1:length(indexofcols);
    i=indexofcols(j);
    signodesonly2(j,:)=signodesonly(i,:);
    end
    %look at the avg value within those 26 node
    avgbetainside=mean(signodesonly2(signodesonly2~=0));
    %make a matrix of everything outside the 26 nodes of significance
    a=edgebetas(:, [indexofcols indexofnonsignodes]);
    nodesandedges=a([indexofcols indexofnonsignodes],:);
    nonsigedgesonly=nodesandedges([5:359],[5:359]);
    %look at the avg value within nonsig node edges only
    avgbetatorest=mean(nonsigedgesonly(nonsigedgesonly~=0));

%what is the true difference between the effect on edges within and
%those in the rest of the brain?
permdiff(x)=avgbetainside-avgbetatorest;

end

%make a histogram and see where the real difference falls on it
hist(permdiff)
hold on
line([truediff2 truediff2], [0 3500])
%get a p-value
perm_edgesbeta_dist=fitdist(permdiff','Normal')
p=cdf(perm_edgesbeta_dist, truediff2)