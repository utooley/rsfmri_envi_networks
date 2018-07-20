%% SETUP
%Running Locally
datadir=fullfile('~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists'
outdir='/Users/utooley/Dropbox (Personal)/bassett_lab/clustco_paper/brains'
parceldir='~/Documents/bassett_lab/tooleyEnviNetworks/parcels/Glasser/'
clustcodir='~/Dropbox/bassett_lab/clustco_paper/'

%Running Locally Bassett
datadir=fullfile('~/Documents/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='~/Documents/tooleyEnviNetworks/subjectLists'
clustcodir='~/Dropbox/bassett_lab/clustco_paper/'
parceldir='~/Documents/tooleyEnviNetworks/parcels/Glasser/'
outdir='~/Dropbox/bassett_lab/clustco_paper/brains'

%Running on the cluster
datadir=fullfile('/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/restNetwork_GlasserPNC/GlasserPNCzNetworks/')
listdir='/data/jag/bassett-lab/tooleyEnviNetworks/subjectLists'
outdir='/data/jag/bassett-lab/tooleyEnviNetworks/analyses'
%read the subject list in without the header
subjlist=csvread(fullfile(listdir,'n1012_healthT1RestExclude_parcels.csv'),1, 0 )

%% Calculate Clustering Coef at Each Node

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
    
    %Node-wise calculation of the clustering coefficient
    %creates a matrix where each row is a subject's node-wise clustco
    avgclustco_both(n,:)=clustering_coef_wu_sign(subfcmat,3);
    %Nodewise calculation of average weight
    %creates a matrix where each row is a subject's node-wise clustco
    avgstrength(n,:)=strengths_und(subfcmat);
    
end

%calculate column means for each node (average clust co across subjects)
node_clust_co_mean=mean(avgclustco_both);
%calculate column sd for each node (sd of clustering coefficient for a node
%across subjects)
node_clust_co_std=std(avgclustco_both);

%write outfile
outfile=dataset(subjlist, avgclustco_both)
export(outfile,'File',fullfile(outdir,'n1012_clust_co_nodewise_by_subj.csv'),'Delimiter',',')

%% plot the node clust co on the brain
%need to account for parcel 52 being missing!!! insert an empty row
%at row 103 which is 52, and make it gray in the
%colortable
node_clust_co_mean=insertrows(node_clust_co_mean', 0, 102)
%get the p-vals/betas for the effect of interest
node_pvals=dlmread('~/Dropbox (Personal)/bassett_lab/clustco_paper/nodewise_estimates_for_ct_effect_corrected.csv', ',', 1, 1)
pvals=node_pvals
pvals=node_pvals(:,2);
%pvals=insertrows(pvals, 1, 102);
%load annotation file
[num_vertices label colortable]=read_annotation(fullfile(parceldir,'lh.HCP-MMP1.annot'));
%[num_vertices label colortable]=read_annotation(fullfile(outdir,'mean_clustco_blue_rh.annot'));
%read in the glasser look-up table
mapping = readtable(fullfile(parceldir,'glasser_lookup_for_matlab.csv'),'ReadVariableNames',false);
%make empty vector for new labels
newlabel=label;
%for each region in the glasser atlas lh annotation
for n=1:numel(unique(colortable.struct_names));
    name=char(colortable.struct_names(n));
    %try
    %the row of mapping that matches the structure name
    ind=strmatch(name, mapping{:,2});
    region_num=mapping{ind,1};
    %find the value/clustco/community assignment for that region in the mapping
    comm_assign=pvals(region_num);
    %go into the colortable struct-names and find the colortable.table(:,5)
    %that matches the vertex ID.
    vertind=find(label == colortable.table(n,5));
    newlabel(vertind)=comm_assign;
    %end
end
%create a new empty struct named copy
copy = struct();
copy.numEntries=numel(unique(newlabel))+1; %num entries=num new labels, plus one for the unknown label needed at the beginning
copy.orig_tab='custom_tab_clustco_nodewise' %name the new colortable
copy.struct_names=cell(numel(unique(newlabel))+1,1); %%make enough structure names for num of labels
ranked=sort(unique(newlabel))%rank the labels (clustering coef by node) from low to high
copy.struct_names{1}='Unknown'  %make the first struct_name the empty unknown
copy.table(1,:)=[0 0 0 0 0];  %make the first colortable entry the empty unknown
for i=1:length(unique(newlabel))
    %make the struct_names the clust co, in order of magnitude
    %make sure this works
    copy.struct_names{i+1}=[char(num2str(ranked(i)))];
end
%make an empty colortable of 0's
copy.table=zeros(numel(unique(newlabel))+1,5)
% 
% % COLORMAP 1
% %create a vector that ranges evenly from 0-255 in increments to match the
% %length of the number of unique clust co values
% Yellow0to255 = round(linspace(0,255,numel(unique(newlabel)))');
% for i=1:length(unique(newlabel));
%     copy.table(i+1,:)=[255 Yellow0to255(i) 0 0 0]; %map that vector to the colortable
% end
% %get the final column of colortable.table by = R + G*2^8 + B*2^16 + flag*2^24
% copy.table(:,5)=copy.table(:,1) + (copy.table(:,2).*(2.^8)) + (copy.table(:,3).*(2.^16))

% COLORMAP 2
clear cmap
n = numel(unique(newlabel));                %// number of colors
% cmap(1,:) = [255 255 240];   %// color first row - off-white
% cmap(2,:) = [255 0 127];   %// color mid row - mid red
% cmap(3,:) = [156 0 76];   %// color mid row - dark red
%reho agexses previous
%cmap(1,:) = [255 255 240];   %// color first row - off-white
%cmap(2,:) = [51 255 255];   %// color mid row - light blue
%cmap(3,:) = [46 139 87];   %// color last row - dark green
%age effect clust co
cmap(1,:) = [51 0 102];   %// color first row - dark purple
cmap(2,:) = [153 51 255];   %// color mid row - mid purple
cmap(3,:) = [229 204 255];   %// color last row - light purple
%[255 255 102];   %// color first row - light yellow
%mean clust co
% cmap(1,:) = [0 0 255];   %// color first row - dark blue
% cmap(2,:) = [0 191 255];   %// color mid row - light blue
% cmap(3,:) = [255 255 240];   %// color last row - off-white
%age x ses
% cmap(1,:) = [255 140 0];   %// color first row - mid-orange
% cmap(2,:)=[255 255 0];   %// color mid row - mid-yellow
% cmap(3,:)=[255 255 102];   %// color last row - light yellow

[X,Y] = meshgrid([1:3],[1:numel(unique(newlabel))]);  %// mesh of indices, number of unique(

cmap = interp2(X([1,60,numel(unique(newlabel))],:),Y([1,60,numel(unique(newlabel))],:),cmap,X,Y); %// interpolate colormap
%  for i=1:5
%       cmap(i,:)=[255 0 0]%%//color all sig p-values dark red, look at ranked to know how far to index
%  end 
cmap(:,4)=zeros(numel(unique(newlabel)),1)
cmap(:,5)=zeros(numel(unique(newlabel)),1)
copy.table(2:numel(unique(newlabel))+1, :)=round(cmap)
copy.table(:,5)=copy.table(:,1) + (copy.table(:,2).*(2.^8)) + (copy.table(:,3).*(2.^16))
copy.table(length(copy.table), :)=[0 0 0 0 0] %make the last entry in the colortable, which will be mapped to subcort regions, 0
%copy.table((length(copy.table)-1), :)=[ 0           0         102           0     6684672] %make the missing parcel the least sig color
% RELABEL
%now transfer the RGB values in the fifth col back to the ranked clust co values for nodes
%want each value of the colortable (low to high) to be mapped to the value
%in the ranked list of clust co values
for n=1:(numel(unique(newlabel)));
    vertind=find(newlabel==ranked(n));
    %vertind=find(newlabel == uniqlabels(n));
    newlabel(vertind)=copy.table(n+1,5);
    %the row of mapping that matches the structure name
end
%write out the annotation file
filename=fullfile(outdir, 'ct_age_effect_linear_lh.annot')

write_annotation(filename, num_vertices, newlabel, copy)

%% CREATE AN IMAGE OF THE COLORBAR
cmap1=round(cmap)/255
figure
colormap(cmap1(1:140, 1:3))
set(1,'Color','w')
h=colorbar
set( h, 'YDir', 'reverse' );

set(h, 'YTick', [])

