% Function to create a code plot with cluster colouring
% each point is labelled with the input cluster mapping
% on to it.
% Guy Billings UCL, 2009

function [dominant_in_cluster]=cluster_plot(out_cluster,code_in,code_out,index_in,index_out)

dims=size(index_in);
clusts=dims(1);

col_tmp=[[1,1,0];[1,0,1];[0,1,1];[1,0,0];[0,1,0];[0,0,1];[0,0,0]];
cols=zeros(clusts,3);

for col_ch=1:clusts
    
    cols(col_ch,:)=col_tmp(randsample(7,1),:);
    cols(col_ch,randsample(3,1))=1-rand;
    
end    


this_clust=index_in(1,find(index_in(1,:)~=0));
this_clust=intersect(this_clust,index_out(out_cluster,find(index_out(out_cluster,:)~=0)));
maxx=max(code_in(this_clust));
maxy=max(code_out(this_clust));
gmin=min([maxx,maxy]);
plot(code_in(this_clust),code_out(this_clust),'color',cols(1,:),'linestyle','.');
hold on

% count through input clusters. For each cluster find the codes who
% intersect with the chosen output cluster. Each such intersection
% is plotted with a differing colour: Therefore all output codes in this
% cluster are labelled with their originating input cluster

last_card=0;
dominant_in_cluster=0;

for c=2:clusts
       
    this_clust=index_in(c,find(index_in(c,:)~=0));
    this_clust=intersect(this_clust,index_out(out_cluster,find(index_out(out_cluster,:)~=0)));
    cardin=size(this_clust);
    
    if cardin(2)>last_card
        dominant_in_cluster=c;   % record the input cluster with the largest contribution to this
        last_card=cardin(2);     % output cluster
    end    
    
    plot(code_in(this_clust),code_out(this_clust),'color',cols(c,:),'linestyle','.')
    title(['Plot of output cluster ' num2str(out_cluster)])
    maxx=max(code_in(this_clust));
    maxy=max(code_out(this_clust));
    this_min=min([maxx,maxy]);
    if this_min>gmin
        gmin=this_min;
    end    
    
end    

plot([1:gmin],[1:gmin],'k')

% Label the data sets with the originating cluster number
lstr='legend(''1''';
for c=2:clusts
    
   lstr=[lstr,',''',num2str(c),''''];
   
end
lstr=[lstr,')'];
eval(lstr);
    
hold off




