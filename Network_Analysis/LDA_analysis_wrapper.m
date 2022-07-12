addpath('/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/scripts/topicnet/')

k=10
cells=importdata('/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/data/CVNdata/liger_4cells_sqrt_19genesrm_varthre0.05_k10/celltype_order.txt')
indir='/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/Results/scCVN_upto1transition/liger_4cells_sqrt_19genesrm_varthre0.05_k10_macs2/pg0.2_pm0.8_pr0.2_maxReg50_b4_bm4/subsample/analysis/lda_TFcellbygene/'
cf='_full'; 
cf='_cf0.8'
cf='_top5k';
%cf='_cf0.8_binary'
prefix='_filteredlowexpression'
dat='Buenrostro'
binry=0
for k=[30,20,10]
    [mdl10]=LDA_analysis(indir,cells,k,prefix,cf,dat,binry)
end

cells=importdata('/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/data/liger_sqrt_ncell50_k10_filterhumanbc/celltype_order.txt')
indir='/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/Results/scCVN_upto1transition/liger_sqrt_ncell50_k10_filterhumanbc_macs2/pg0.2_pm0.8_pr0.2_maxReg50_b4_bm4/subsample/analysis/lda_TFcellbygene/'
cf='_full'; 
cf='_cf0.8';
cf='_top5k';
%cf='_cf0.8_binary'
prefix='_filteredlowexpression'
dat='A2S'
binry=0
for k=[30,20,10]
    [mdl10]=LDA_analysis(indir,cells,k,prefix,cf,dat,binry)
end

cells=importdata('/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/data/liger_sqrt_ncell50_k8_filterhumanbc_FBS/celltype_order.txt')
indir='/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/Results/scCVN_upto1transition/liger_sqrt_ncell50_k8_filterhumanbc_FBS_macs2/pg0.2_pm0.8_pr0.2_maxReg50_b4_bm4/subsample/analysis/lda_TFcellbygene/'
cf='_full'; 
cf='_cf0.8';
cf='_top4k';
prefix='_filteredlowexpression'
dat='FBS'
binry=0
for k=[30,20,10]
    [mdl10]=LDA_analysis(indir,cells,k,prefix,cf,dat,binry)
end



k=10
outdir=sprintf('%s/network%s/k%d/',indir,cf,k)
outdir=sprintf('%s/k%d/',indir,k)
load(sprintf('%s/lda_model_k%d_%s%s%s.mat',outdir,k,dat,prefix,cf))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/scripts/topicnet/LDA_analysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indir=sprintf('%s/network%s/',indir,cf)

d=[]
st=[]
en=[]
regulators=[]
idrm=[]
for i=1:length(cells)
    cells{i}
    dt=importdata(sprintf('%s/%s_consensus_edges%s_mat%s.txt',indir,cells{i},prefix,cf));
    st=[st,size(d,1)+1];
    d=[d;dt.data];  %document*words TF-cell*targets
    en=[en,size(d,1)];
    regulators=[regulators;dt.textdata(2:end,1)];
end
rowsum=sum(d,2);
idrms=find(rowsum==0);

%save(sprintf('%s/consensus_edges_data.mat',indir),'d','-v7.3');

n=size(dt.data,1); % tf
v=size(dt.data,2); % gene
regnames=dt.textdata(2:end,1);
gnames=dt.textdata(1,2:end)';


%load(sprintf('%s/consensus_edges_data.mat',indir));
%for k=[10,15,20]
for k=[10,20,30]
    rng('default')
    outdir=sprintf('%s/k%d/',indir,k)
    mkdir(outdir)
    if binry==1
        mdl10=fitlda(d,k);
    else
        mdl10=fitlda((ceil(d*100)),k);
    end
    save(sprintf('%s/lda_model_k%d_%s%s%s.mat',outdir,k,dat,prefix,cf),'mdl10','-v7.3');
end
load(sprintf('%s/lda_model_k%d_%s%s%s.mat',outdir,k,dat,prefix,cf))

for seed=0:100
    rng(seed)
    mdl10=fitlda((ceil(d*100)),k);
    if 15597715.043696<mdl10.FitInfo.NegativeLogLikelihood && mdl10.FitInfo.NegativeLogLikelihood<15597715.043698
        break
    end
end

docprob=mdl10.DocumentTopicProbabilities; %rowsum is 1 sum(mdl10.DocumentTopicProbabilities,2)=1
%docprob_norm=docprob./repmat(sum(docprob,2),1,size(docprob,2));
[ig,topicid]=max(docprob,[],2);
[igt,tforder]=sort(topicid);
f=figure;
imagesc(docprob(tforder,:));
title('Document-Topic Probabilities');
colorbar;
ylabel('TFs');
set(gcf,'PaperPosition',[ 0 0 5 5], 'PaperPositionMode','manual', 'PaperSize',[5 5]);%  W H
saveas(gcf,sprintf('%s/lda_k%d_TF-topic%s_cells%s%s.pdf',outdir,k,dat,prefix,cf),'pdf');

[TopicWordProb,geneid]=max(mdl10.TopicWordProbabilities,[],2);
[~,gorder]=sort(geneid);
f2=fopen(sprintf('%s/genes_topicid.txt',outdir),'w');   
for ii=1:v
    fprintf(f2,'%s\t%d\t%f\n',gnames{ii},geneid(ii),TopicWordProb(ii));
end
fclose(f2);

f=figure;
imagesc(mdl10.TopicWordProbabilities(gorder,:));
title('Topic-gene Probabilities');
colorbar;
ylabel('Genes');
set(gcf,'PaperPosition',[ 0 0 5 5], 'PaperPositionMode','manual', 'PaperSize',[5 5]);%  W H
saveas(gcf,sprintf('%s/lda_k%d_topic-gene%s_allcells%s%s.pdf',outdir,k,dat,prefix,cf),'pdf');

%end
%clear mdl10

% get tf-topic per cell type:

topicid=[];
upp=0.5 %quantile(mdl10.DocumentTopicProbabilities(:),0.95);
f=figure;
for i=1:length(cells)
    id=st(i):en(i);
    st(i)
    en(i)
    idrm=intersect(idrms,id);
    idkp=setdiff(id,idrm);
    d1prob=mdl10.DocumentTopicProbabilities(idkp,:);  % rowsum is 1, sum(mdl10.DocumentTopicProbabilities,2)=1                            
    %d1prob_norm=d1prob./repmat(sum(d1prob,2),1,size(d1prob,2)); %
    [TopicProb,tfid]=max(d1prob,[],2);
    topicid=[topicid;tfid];
    
    subplot(2,4,i);
    [ig,tforder]=sort(tfid);
    imagesc(d1prob(tforder,:),[0 upp]);
    title(cells{i});
    colorbar;

    %% output TFs-topic id  

    fid=fopen(sprintf('%s/TFs_topicid_%s.txt',outdir,cells{i}),'w');
    for ii=1:size(d1prob,1)
        %if length(idrm)>0 && find(idrm==ii)
        %    continue
        %end
        fprintf(fid,'%s\t%d\t%f\n',regnames{ii},tfid(ii),TopicProb(ii));
    end
    fclose(fid);

end
set(gcf,'PaperPosition',[ 0 0 20 10], 'PaperPositionMode','manual', 'PaperSize',[20 10]);%  W H
saveas(gcf,sprintf('%s/lda_k%d_%s_tf-topicPercell%s%s.pdf',outdir,k,dat,prefix,cf),'pdf');

%[ig,rorder]=sortrows(topicid);
%topicid1=topicid(rorder,:);

% rewiring networks:
%/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/scripts/topicnet/KLDiv.m
addpath('/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/scripts/topicnet/')
regs=regnames;
for i=1:length(cells)
    id=st(i):en(i);    
    %id=((i-1)*n+1):n*i;
    d1prob=mdl10.DocumentTopicProbabilities(id,:);                             
    %d1prob_reg=d1prob_norm(regids,:);
    data1_reg=d(id,:);
    %data1_reg=d(regids,id);
    data1_reg_degree=sum(data1_reg,2); 
    %reg1=regulators(id);
    %nameIDMap1=containers.Map(reg1,1:size(reg1,1));
    %tf1=[]
    %tf2=[]
    for j=(i+1):length(cells)
        id2=st(j):en(j);
        %id2=((i-1)*n+1):n*i;
        d2prob=mdl10.DocumentTopicProbabilities(id2,:);  
        %reg2=regulators(id2);       
        %nameIDMap2=containers.Map(reg2,1:size(reg2,1));                             
        %d2prob_norm=d2prob./repmat(sum(d2prob,2),1,size(d2prob,2));
        %d2prob_reg=d2prob_norm(regids,:);
        data2_reg=d(id2,:);
        %data2_reg=d(regids,id2);
        data2_reg_degree=sum(data2_reg,2); 
        % regs=intersect(reg1,reg2);
        % for r=1:size(regs,1)
        %     tf1=[tf1,nameIDMap1(regs{r})];
        %     tf2=[tf2,nameIDMap2(regs{r})];
        % end
        % d1prob=d1prob(tf1,:);
        % d2prob=d2prob(tf2,:);
        % data1_reg_degree=sum(data1_reg(tf1,:),2); 
        % data2_reg_degree=sum(data2_reg(tf2,:),2); 

        %ddist=pdist2(d1prob,d2prob);
        ddist1=KLDiv(d1prob,d2prob);
        ddist2=KLDiv(d2prob,d1prob);
        %mydist=diag(ddist);
        mydist=(ddist1+ddist2)/2;
        [ig,rorder]=sort(mydist,'descend');
        fid=fopen(sprintf('%s/%s_%s_rewired.txt',outdir,cells{i},cells{j}),'w');
        %%gnames=importdata('/mnt/dv/wid/projects5/Roy-singlecell/js_work/TMF/Shilu_nets/inputs/dense/union_gnames.txt');
        %%gnames=dt1.textdata;
        for ii=1:length(rorder)
           fprintf(fid,'%s\t%f\n',regs{rorder(ii)},mydist(rorder(ii)));
        end
        fclose(fid);
        %mydist_reg=mydist(regids);
        %[ig,rorder_reg]=sort(mydist_reg,'descend');
        upp=max(max(max(d1prob(rorder(1:ntop),:)),max(d2prob(rorder(1:ntop),:))))
        ntop=50
        figure;
        subplot(1,5,1);
        imagesc(d1prob(rorder(1:ntop),:),[0,upp]);
        set(gca,'ytick',[1:ntop]);
        set(gca,'yticklabel',regs(rorder(1:ntop)));
        set(gca,'fontsize',6);
        title(sprintf('%s regulators',cells{i}));
        colorbar;
        subplot(1,5,2);
        imagesc(d2prob(rorder(1:ntop),:),[0,upp]);
        set(gca,'ytick',[1:ntop]);
        set(gca,'yticklabel',regs(rorder(1:ntop)));
        set(gca,'fontsize',6);
        title(sprintf('%s regulators',cells{j}));
        colorbar;
        subplot(1,5,3);
        plot(mydist(rorder(1:ntop)),'o-');
        title(sprintf('regulator rewiring score:%s vs %s',cells{i},cells{j}));
        box(gca,'off');
        %set(gcf,'PaperPosition',[ 0 0 10 ntop/10], 'PaperPositionMode','manual', 'PaperSize',[10 ntop/10]);%  W H
        %saveas(gcf,sprintf('%s/lda_k%d_%s_%s_regulator_rewiring_top%d.pdf',indir,k,cells{i},cells{j},ntop),'pdf');

        maxx=max(max(data1_reg_degree(rorder(1:ntop))),max(data2_reg_degree(rorder(1:ntop))))+1
        
        %figure;
        subplot(1,5,4);
        barh(data1_reg_degree(rorder(1:ntop)));
        xlim([0,maxx])
        set(gca,'ytick',[1:ntop]);
        set(gca,'yticklabel',regs(rorder(1:ntop)));
        set(gca,'fontsize',6);
        title(sprintf('%s regulators',cells{i}));
        subplot(1,5,5);
        barh(data2_reg_degree(rorder(1:ntop)));
        xlim([0,maxx])
        set(gca,'ytick',[1:ntop]);
        set(gca,'yticklabel',regs(rorder(1:ntop)));
        set(gca,'fontsize',6);
        title(sprintf('%s regulators',cells{j}));
        set(gcf,'PaperPosition',[ 0 0 20 6], 'PaperPositionMode','manual', 'PaperSize',[20 6]);%  W H
        saveas(gcf,sprintf('%s/lda_k%d_%s_%s_top%dregulator_heatmap_numberoftargets.pdf',outdir,k,cells{i},cells{j},ntop),'pdf');
        % rewiring genes:
        ntop=50
        figure;
        subplot(1,3,1);
        imagesc(d1prob_norm(rorder(1:ntop),:));
        set(gca,'ytick',[1:ntop]);
        set(gca,'yticklabel',gnames(rorder(1:ntop)));
        set(gca,'fontsize',6);
        title(sprintf('%s regulators',cells{i}));
        colorbar;
        subplot(1,3,2);
        imagesc(d2prob_norm(rorder(1:ntop),:));
        set(gca,'ytick',[1:ntop]);
        set(gca,'yticklabel',gnames(rorder(1:ntop)));
        set(gca,'fontsize',6);
        title(sprintf('%s regulators',cells{j}));
        colorbar;
        subplot(1,3,3);
        plot(mydist(rorder(1:ntop)),'o-');
        title(sprintf('rewiring gene score:%s vs %s',cells{i},cells{j}));
        box(gca,'off');
        set(gcf,'PaperPosition',[ 0 0 10 5], 'PaperPositionMode','manual', 'PaperSize',[10 5]);%  W H
        saveas(gcf,sprintf('%s/lda_k%d_%s_%s_top%dregulator_heatmap_numberoftargets.pdf',indir,k,cells{i},cells{j},ntop),'pdf');

    end
end
find(contains(regnames(rorder),'Nanog'))
find(contains(regnames(rorder),'PROCR'))
find(contains(regnames(rorder),'IRF8'))
find(contains(regnames(rorder),'CSF1R'))
find(contains(regnames(rorder),'CTSG'))
find(contains(regnames(rorder),'MPO'))

find(contains(regnames,'Nanog'))

row=[]
for r=1:n
for i=1:k
    if all(topicid1(r,:)==i)
    row=[row;r]
    end
end
end

regnames1=regnames(rorder)
regnames2=regnames1
topicid2=topicid1
topicid2(row,:)=[]
regnames2(row)=[]
f=figure;
imagesc(topicid2);
set(gca,'xtick',1:length(cells));
set(gca,'xticklabel',cells);
set(gca,'ytick',1:n,'TickLength',[0.001, 0.01]);
set(gca,'yticklabel',regnames1);
set(gca,'fontsize',5);
colorbar;
set(gcf,'PaperPosition',[ 0 0 8 16], 'PaperPositionMode','manual', 'PaperSize',[8 16]);%  W H
saveas(gcf,sprintf('%s/lda_k%d_Buenrostro_cells%s%s_topicID_divergent.pdf',indir,k,prefix,cf),'pdf');

f=figure;
imagesc(topicid1(row,:));
set(gca,'xtick',1:length(cells));
set(gca,'xticklabel',cells);
set(gca,'ytick',1:n,'TickLength',[0.001, 0.01]);
set(gca,'yticklabel',regnames1(row));
set(gca,'fontsize',5);
colorbar;
set(gcf,'PaperPosition',[ 0 0 8 8], 'PaperPositionMode','manual', 'PaperSize',[8 8]);%  W H
saveas(gcf,sprintf('%s/lda_k%d_Buenrostro_cells%s%s_topicID_common.pdf',indir,k,prefix,cf),'pdf');


idtf=[]
for t=1:length(regnames)
        Index = find(ismember(gnames,regnames{t}));
        idtf=[idtf,Index];
end

figure;
for i=1:length(cells)
    subplot(2,4,i);
    id=((i-1)*n+1):n*i;
    d1prob=mdl10.TopicWordProbabilities(id,:);                                      
    d1prob_norm=d1prob./repmat(sum(d1prob,2),1,size(d1prob,2));
    [ig,d1id]=max(d1prob_norm,[],2);
    [ig,order]=sort(d1id);
    d1=d(:,id);
    spy(d1(geneorder,order),'r',2);
    title(cells{i});
end
set(gcf,'PaperPosition',[ 0 0 20 10], 'PaperPositionMode','manual', 'PaperSize',[20 10]);%  W H
saveas(gcf,sprintf('%s/lda_k%d_adjacencymatrix_Buenrostro_cells.pdf',indir,k),'pdf');

regids=importdata(sprintf('%s/regulator_id.txt',indir));
%regnames=regids.textdata;
gnames=dt.textdata;
regnames=gnames(regids);

figure;
for i=1:length(cells)
    subplot(2,4,i);
    id=((i-1)*n+1):n*i;
    d1prob=mdl10.TopicWordProbabilities(id,:);                                      
    d1prob_norm=d1prob./repmat(sum(d1prob,2),1,size(d1prob,2));
    d1prob_reg=d1prob_norm(regids,:);
    [ig,d1d_reg]=max(d1prob_reg,[],2);
    [ig,order_reg]=sort(d1d_reg);
    imagesc(d1prob_reg(order_reg,:),[0 upp]);
    title(sprintf('%s regulators',cells{i}));
    colorbar;
end
set(gcf,'PaperPosition',[ 0 0 20 10], 'PaperPositionMode','manual', 'PaperSize',[20 10]);%  W H
%saveas(gcf,sprintf('%s/lda_k%d_A2S_allreg.pdf',indir,k),'pdf');
saveas(gcf,sprintf('%s/lda_k%d_Buenrostro_allreg.pdf',indir,k),'pdf');






f=figure;
subplot(1,3,1);
imagesc(d1prob_reg(rorder_reg(1:ntop),:));
set(gca,'ytick',[1:ntop]);
set(gca,'yticklabel',regnames(rorder_reg(1:ntop)));
set(gca,'fontsize',5);
title(sprintf('%s regulators',cells{i}));
colorbar;
subplot(1,3,2);
imagesc(d2prob_reg(rorder_reg(1:ntop),:));
set(gca,'ytick',[1:ntop]);
set(gca,'yticklabel',regnames(rorder_reg(1:ntop)));
set(gca,'fontsize',5);
title(sprintf('%s regulators',cells{j}));
colorbar;
subplot(1,3,3);
plot(mydist_reg(rorder_reg(1:ntop)),'o-');
title(sprintf('regulator rewiring score:%s vs %s',cells{i},cells{j}));
box(gca,'off');
set(gcf,'PaperPosition',[ 0 0 10 ntop/10], 'PaperPositionMode','manual', 'PaperSize',[10 ntop/10]);%  W H
saveas(gcf,sprintf('%s/lda_k%d_%s_%s_regulator_rewiring_top%d.pdf',indir,k,cells{i},cells{j},ntop),'pdf');

f=figure;
subplot(1,2,1);
barh(data1_reg_degree(rorder_reg(1:ntop)));
set(gca, 'YDir','reverse')
set(gca,'ytick',[1:ntop]);
set(gca,'yticklabel',regnames(rorder_reg(1:ntop)));
set(gca,'fontsize',5);
title(sprintf('%s regulators:# of targets',cells{i}));
subplot(1,2,2);
barh(data2_reg_degree(rorder_reg(1:ntop)));
set(gca, 'YDir','reverse')
set(gca,'ytick',[1:ntop]);
set(gca,'yticklabel',regnames(rorder_reg(1:ntop)));
set(gca,'fontsize',5);
title(sprintf('%s regulators:# of targets',cells{j}));
set(gcf,'PaperPosition',[ 0 0 10 ntop/10], 'PaperPositionMode','manual', 'PaperSize',[10 ntop/10]);%  W H
saveas(gcf,sprintf('%s/lda_k%d_%s_%s_regulator_rewiring_top%d_ntargets.pdf',indir,k,cells{i},cells{j},ntop),'pdf');



plot(mydist_reg(rorder_reg(1:200)),'o-');
title('regulator rewiring score');
box(gca,'off');
saveas(gcf,sprintf('%s/regulator_rewiring_top200.pdf',indir),'pdf');


dtops=[]
for i=1:length(cells)
    id=((i-1)*n+1):n*i;
    d1prob=mdl10.TopicWordProbabilities(id,:);                                      
    d1prob_norm=d1prob./repmat(sum(d1prob,2),1,size(d1prob,2));
    [ig,d1id]=max(d1prob_norm,[],2);
    dtops=[dtops,d1id];
end

dtopsregs=dtops(regids,:);
f=figure;
imagesc(dtopsregs(1:100,:));
set(gca,'ytick',[1:100]);
set(gca,'TickLength',[0, 0])
set(gca,'yticklabel',regnames(1:100));
set(gca,'xticklabel',cells);
set(gca,'fontsize',6);
title(sprintf('Buenrostro regulators topics with max probablity'));  
colorbar;
set(gcf,'PaperPosition',[ 0 0 8 10], 'PaperPositionMode','manual', 'PaperSize',[8 10]);%  W H
saveas(gcf,sprintf('%s/lda_k%d_maxcluster_Buenrostro_cells.png',indir,k),'png');

