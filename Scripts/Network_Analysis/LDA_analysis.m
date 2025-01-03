function [mdl10]=LDA_analysis(indir,cells,k,prefix,cf,dat,binry)

cells=importdata(cells)

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


n=size(dt.data,1); % tf
v=size(dt.data,2); % gene
regnames=dt.textdata(2:end,1);
gnames=dt.textdata(1,2:end)';


rng('default')
outdir=sprintf('%s/k%d/',indir,k)
mkdir(outdir)
if binry==1
    mdl10=fitlda(d,k);
else
    mdl10=fitlda((ceil(d*100)),k);
end
save(sprintf('%s/lda_model_k%d_%s%s%s.mat',outdir,k,dat,prefix,cf),'mdl10','-v7.3');


f=figure;
subplot(1,2,1)
imagesc(mdl10.DocumentTopicProbabilities);
title('Document-Topic Probabilities');
colorbar;
ylabel('TFs * cells');
subplot(1,2,2)
imagesc(mdl10.TopicWordProbabilities);
title('Topic-gene Probabilities');
colorbar;
ylabel('Genes');
set(gcf,'PaperPosition',[ 0 0 10 5], 'PaperPositionMode','manual', 'PaperSize',[10 5]);%  W H
saveas(gcf,sprintf('%s/lda_k%d_%s_allcells%s%s.pdf',outdir,k,dat,prefix,cf),'pdf');


% plot Document-Topic Probabilities sorted:
docprob=mdl10.DocumentTopicProbabilities; %rowsum is 1 sum(mdl10.DocumentTopicProbabilities,2)=1
[~,topicid]=max(docprob,[],2);
[~,tforder]=sort(topicid);
f=figure;
subplot(1,2,1)
imagesc(docprob(tforder,:));
title('Document-Topic Probabilities');
colorbar;
ylabel('TFs * cells');


% plot Topic-gene Probabilities
[TopicWordProb,geneid]=max(mdl10.TopicWordProbabilities,[],2);
[~,gorder]=sort(geneid);
subplot(1,2,2)
imagesc(mdl10.TopicWordProbabilities(gorder,:));
title('Topic-gene Probabilities');
colorbar;
ylabel('Genes');
set(gcf,'PaperPosition',[ 0 0 10 5], 'PaperPositionMode','manual', 'PaperSize',[10 5]);%  W H
saveas(gcf,sprintf('%s/lda_k%d_topic-gene%s_allcells%s%s.pdf',outdir,k,dat,prefix,cf),'pdf');


% Topic-gene Probabilities
[TopicWordProb,geneid]=max(mdl10.TopicWordProbabilities,[],2);
f2=fopen(sprintf('%s/genes_topicid.txt',outdir),'w');   
for ii=1:length(gnames)
    fprintf(f2,'%s\t%d\t%f\n',gnames{ii},geneid(ii),TopicWordProb(ii));
end
fclose(f2);


% get tf-topic per cell type:
docprob=mdl10.DocumentTopicProbabilities; %rowsum is 1 sum(mdl10.DocumentTopicProbabilities,2)=1
[TopicProb,TFtopicid]=max(docprob,[],2);
upp=0.5 %quantile(mdl10.DocumentTopicProbabilities(:),0.95);
f=figure;
for i=1:length(cells)
    id=st(i):en(i);
    st(i)
    en(i)
    idrm=intersect(idrms,id);
    idkp=setdiff(id,idrm);

    %% output TFs-topic id  
    fid=fopen(sprintf('%s/TFs_topicid_%s.txt',outdir,cells{i}),'w');
    for ii=idkp
        fprintf(fid,'%s\t%d\t%f\n',regulators{ii},TFtopicid(ii),TopicProb(ii));
    end
    fclose(fid);

    d1prob=mdl10.DocumentTopicProbabilities(idkp,:);  % rowsum is 1, sum(mdl10.DocumentTopicProbabilities,2)=1     
    tfs=regulators(idkp);                       
    [~,tfid]=max(d1prob,[],2);
    
    subplot(2,4,i);
    [~,tforder]=sort(tfid);
    imagesc(d1prob(tforder,:),[0 upp]);
    title(cells{i});
    colorbar;

end
set(gcf,'PaperPosition',[ 0 0 20 10], 'PaperPositionMode','manual', 'PaperSize',[20 10]);%  W H
saveas(gcf,sprintf('%s/lda_k%d_%s_tf-topicPercell%s%s.pdf',outdir,k,dat,prefix,cf),'pdf');


regs=regnames;
for i=1:length(cells)
    id=st(i):en(i);    
    d1prob=mdl10.DocumentTopicProbabilities(id,:);                             
    data1_reg=d(id,:);
    data1_reg_degree=sum(data1_reg,2); 
    for j=(i+1):length(cells)
        id2=st(j):en(j);
        d2prob=mdl10.DocumentTopicProbabilities(id2,:);  
        data2_reg=d(id2,:);
        data2_reg_degree=sum(data2_reg,2); 

        ddist1=KLDiv(d1prob,d2prob);
        ddist2=KLDiv(d2prob,d1prob);
        mydist=(ddist1+ddist2)/2;
        [ig,rorder]=sort(mydist,'descend');
        fid=fopen(sprintf('%s/%s_%s_rewired.txt',outdir,cells{i},cells{j}),'w');
        for ii=1:length(rorder)
           fprintf(fid,'%s\t%f\n',regnames{rorder(ii)},mydist(rorder(ii)));
        end
        fclose(fid);
        ntop=50;
        upp=max(max(max(d1prob(rorder(1:ntop),:)),max(d2prob(rorder(1:ntop),:))))
        figure;
        subplot(1,5,1);
        imagesc(d1prob(rorder(1:ntop),:),[0,upp]);
        set(gca,'ytick',[1:ntop]);
        set(gca,'yticklabel',regnames(rorder(1:ntop)));
        set(gca,'fontsize',6);
        title(sprintf('%s regulators',cells{i}));
        colorbar;
        subplot(1,5,2);
        imagesc(d2prob(rorder(1:ntop),:),[0,upp]);
        set(gca,'ytick',[1:ntop]);
        set(gca,'yticklabel',regnames(rorder(1:ntop)));
        set(gca,'fontsize',6);
        title(sprintf('%s regulators',cells{j}));
        colorbar;
        subplot(1,5,3);
        plot(mydist(rorder(1:ntop)),'o-');
        title(sprintf('regulator rewiring score:%s vs %s',cells{i},cells{j}));
        box(gca,'off');

        maxx=max(max(data1_reg_degree(rorder(1:ntop))),max(data2_reg_degree(rorder(1:ntop))))+1
        
        %figure;
        subplot(1,5,4);
        barh(data1_reg_degree(rorder(1:ntop)));
        xlim([0,maxx])
        set(gca,'ytick',[1:ntop]);
        set(gca,'yticklabel',regnames(rorder(1:ntop)));
        set(gca,'fontsize',6);
        set(gca,'Ydir','reverse')
        title(sprintf('%s regulators',cells{i}));
        subplot(1,5,5);
        barh(data2_reg_degree(rorder(1:ntop)));
        xlim([0,maxx])
        set(gca,'ytick',[1:ntop]);
        set(gca,'yticklabel',regnames(rorder(1:ntop)));
        set(gca,'fontsize',6);
        title(sprintf('%s regulators',cells{j}));
        set(gca,'Ydir','reverse')
        set(gcf,'PaperPosition',[ 0 0 20 6], 'PaperPositionMode','manual', 'PaperSize',[20 6]);%  W H
        saveas(gcf,sprintf('%s/lda_k%d_%s_%s_top%dregulator_heatmap_numberoftargets.pdf',outdir,k,cells{i},cells{j},ntop),'pdf');

    end
end

end

