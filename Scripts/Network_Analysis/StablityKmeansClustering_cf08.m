%% K-means on confidence matrix of stability selection:


path='Results_subsample/analysis/'
n=0.8
prefix='testdata'


% input: prefix, path, n, kend
d=importdata(sprintf('%s/consensus_edges_merged_cf%g.txt',path,n));
sprintf('%s/consensus_edges_merged_cf%g.txt',path,n)
outdir=sprintf('%s/kmeansclustering_cf%g/',path,n);


if ~exist(outdir)
    mkdir(outdir)
end


data=d.data;
font=7
name=d.textdata(2:end,1);
cells=d.textdata(1,2:end)
cellabel=strrep(cells,'cluster','C')
SI=[];
smat=[];
kstart=1 %length(cells)+1
kend=30
distortion=zeros(kend,1);


for k=kstart:kend
    rand('seed',0);
    [rid,C,sumd,D]=kmeans(data,k,'replicates',100);
    distortion(k,1)=sum(sumd);
    [a,b]=hist(rid,unique(rid));
    a
    s=silhouette(data,rid); %[s,h]=silhouette(data,rid);
    SIk=mean(s)
    SI=[SI,SIk]
    smat=[smat,s];

    % output the cluster id:
    B=array2table(data,'VariableNames',d.textdata(1,2:end));
    Cl=table(rid,'VariableNames',{'cluster'});
    A = table(name,'VariableNames',{'Edge'});
    out = horzcat(A,B,Cl);
    writetable(out,sprintf('%s/%s_stability_kmeans_cf%g_k%d.txt',outdir,prefix,n,k),'Delimiter','\t','WriteRowNames',false)
    [ig,rorder]=sort(rid);
    rparts=[];
    oldcid=rid(rorder(1));
    for i=2:size(rorder,1)
        newcid=rid(rorder(i));
        if(oldcid~=newcid)
            rparts=[rparts i];
            oldcid=newcid;
        end
    end
    
    cmap=[ones(101,1),(1:-0.01:0)',(1:-0.01:0)'];

    f=figure;
    imagesc(data(rorder,:),[0 1]);
    for i=1:length(rparts)
        line([0 size(data,2)+2],[rparts(i) rparts(i)],'color','black','LineWidth',0.5);
    end    
    ylabel('Edges','FontSize',font);
    xlabel('Cells','FontSize',font);   
    title(sprintf('%s k=%d',strrep(prefix,'_','-'),k),'FontSize',font);   
    set(gca,'ytick',[1,rparts]+10,'yticklabel',1:k,'TickLength',[0 0]);
    set(gca,'xtick',1:length(cells),'xticklabel',cellabel,'FontSize',font,'TickLength',[0 0]);  %'XTickLabelRotation',90,
    %
    colormap(cmap)
    colorbar
    set(gcf,'PaperPosition',[ 0 0 4 5], 'PaperPositionMode','manual', 'PaperSize',[4 5]);%  W H
    saveas(gcf,sprintf('%s/%s_stability_kmeans_cf%g_k%d.pdf',outdir,prefix,n,k),'pdf');
end

skyblue=[0.3010, 0.7450, 0.9330]

SI(1)=0
[SIsorted,I] = sort(SI,'descend')
optK2=I(1)
f=figure;
subplot(2,1,1)
krange=I(1:kend-kstart+1)
bp=bar(SIsorted(1:kend-kstart+1)) 
ylim([0,1])
set(bp,'FaceColor',skyblue);
ylabel('Silhouette Index');
xlim([0,kend+1]);
set(gca,'xticklabel',krange,'xtick',1:kend-kstart+1);  

variance=distortion(1:end-1)-distortion(2:end);
distortion_percent=cumsum(variance)/(distortion(1)-distortion(end));
subplot(2,1,2)
plot(distortion_percent,'b*--');
ylabel('Percentage of variance');
set(gca,'xticklabel',kstart+1:kend,'xtick',1:kend-kstart+1);  
set(gcf,'PaperPosition',[ 0 0 10 10], 'PaperPositionMode','manual', 'PaperSize',[10 10]);%  W H 
saveas(gcf,sprintf('%s/%s_kmeans_silhouette.pdf',outdir,prefix),'pdf');


[r,~]=find(distortion_percent>0.95); % 95% k=22, 96%:25, 97%:27, 98%:31 99% k=35
optK=r(1,1)+1
optK2=I(1)
sprintf('optimal K based on SI is %d, K explained at least 95%% variance is %d',optK2,optK)

