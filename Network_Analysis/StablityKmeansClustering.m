%% K-means on confidence matrix of stability selection:

prefix='liger_sqrt_ncell50_k10_filterhumanbc_macs2'
path='/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/Results/scCVN_upto1transition/liger_sqrt_ncell50_k10_filterhumanbc_macs2/pg0.2_pm0.8_pr0.2_maxReg50_b4_bm4/subsample/analysis/'
n=5000

prefix='liger_sqrt_ncell50_k8_filterhumanbc_FBS_macs2'
path='/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/Results/scCVN_upto1transition/liger_sqrt_ncell50_k8_filterhumanbc_FBS_macs2/pg0.2_pm0.8_pr0.2_maxReg50_b4_bm4/subsample/analysis/'
n=4000

path='/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/Results/scCVN_upto1transition/liger_4cells_sqrt_19genesrm_varthre0.05_k10_macs2/pg0.2_pm0.8_pr0.2_maxReg50_b4_bm4/subsample/analysis/'
n=5000
prefix='Buenrostro'

% input: prefix, path, n, kend
d=importdata(sprintf('%s/consensus_edges_stability_selection_filteredlowexpression_top%d.txt',path,n));
sprintf('%s/consensus_edges_stability_selection_filteredlowexpression_top%d.txt',path,n)
outdir=sprintf('%s/kmeansclustering_lowexpression/',path);
if ~exist(outdir)
    mkdir(outdir)
end
data=d.data;
font=7
name=d.textdata(2:end,1);
cells=d.textdata(1,2:end)
cellabel=strrep(cells,'cluster','C') %strrep(cells,'score_','')
SI=[];
smat=[];
kstart=1 %length(cells)+1
kend=30
distortion=zeros(kend,1);
for k=kstart:kend
    k
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
    writetable(out,sprintf('%s/%s_stability_kmeans_top%d_k%d.txt',outdir,prefix,n,k),'Delimiter','\t','WriteRowNames',false)
    %save(sprintf('%s/kmeans_rid_top%d_k%d.mat',outdir,n,k), 'rid','-v7.3')
    
    % sort by row id 1-k:
    [ig,rorder]=sort(rid);
    %rparts=getparts(rid,rorder);
    %function parts=getparts(rid,rorder)
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
    %set(gca,'ytick',1:length(name),'yticklabel',name(rorder),'TickLength',[0.001, 0.01]);   
    set(gca,'ytick',[1,rparts]+10,'yticklabel',1:k,'TickLength',[0 0]);
    set(gca,'xtick',1:length(cells),'xticklabel',cellabel,'FontSize',font,'TickLength',[0 0]);  %'XTickLabelRotation',90,
    %
    colormap(cmap)
    colorbar
    set(gcf,'PaperPosition',[ 0 0 4 5], 'PaperPositionMode','manual', 'PaperSize',[4 5]);%  W H
    set(gcf,'PaperPosition',[ 0 0 4 10], 'PaperPositionMode','manual', 'PaperSize',[4 10]);%  W H
    %saveas(gcf,sprintf('%s_stability80_kmeans_k%d.pdf',prefix,k),'pdf');
    %saveas(gcf,sprintf('%s_stability80_kmeans_k%d.png',prefix,k),'png');

    saveas(gcf,sprintf('%s/%s_stability_kmeans_top%d_k%d-1.pdf',outdir,prefix,n,k),'pdf');
    rparts(2:end)-rparts(1:end-1)

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
%set(gcf,'PaperPosition',[ 0 0 8 4], 'PaperPositionMode','manual', 'PaperSize',[8 4]);%  W H 
%saveas(gcf,sprintf('%s/%s_kmeans_silhouette.pdf',path,prefix),'pdf');

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

