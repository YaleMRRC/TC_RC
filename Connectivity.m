clear all
clc

addpath(genpath('/Users/Mehraveh/Documents/MATLAB/BCT'))
addpath(genpath('/Users/Mehraveh/Documents/MATLAB/Connectivity'))

DataSet = input('what is the data set to use? (''nback'', ''gradCPT'', ''ANT'', ''HCP'', ''Random'') \n','s')

if ~isequal(DataSet,'Random')
    %%% Load data sets and behavior
    load (['/Users/Mehraveh/Documents/MATLAB/Connectivity/Datasets/',DataSet,'_gsr.mat']);
    load (['/Users/Mehraveh/Documents/MATLAB/Connectivity/Datasets/subjname_',DataSet,'.mat']);
    load (['/Users/Mehraveh/Documents/MATLAB/Connectivity/Behaviors/behavior_',DataSet,'.mat']);

elseif isequal(DataSet, 'Random')    
    for subj = 1:20
        for st=1:5
            V = rand(268,1000);
            norm_V = (V-mean(V,2))./std(V,1,2);
            m = real(atanh(corr(norm_V',norm_V')));
            m = m-diag(diag(m));
            m(isnan(m))=0;    
            M{subj,st} = m;
        end
    end
end

l = length(M);

%%
if isequal(DataSet,'nback')
    all_behav = dprime(:,1);
    behav = dprime(:,1);
    behav_label = 'd''';
    subM = 4;
    subN = 7;
elseif isequal(DataSet,'gradCPT')
    all_behav = dprime(:,1);
    behav = dprime(:,1);
    behav_label = 'd''';
    subM = 3;
    subN = 6;
elseif isequal(DataSet,'ANT')
%     behav_input = input('What behavior you want to use for ANT? (''std'', ''cv'') \n','s')
    behav_input = 'cv';
    all_behav = Output(:,2:6);
    if isequal(behav_input,'std')        
        behav = Output(:,3);    
        behav_label = 'Std RT';
    elseif isequal(behav_input,'cv')
        behav = Output(:,3)./Output(:,2);
        behav_label = 'RT CV';
    end    
    subM = 5;
    subN = 8;
%     HOM_size = input('What size HOM analysis to do? (''5'', ''7'') \n')
    HOM_size = 7;
end

metrics = {'Strength', 'Binarized Degree', 'Strength_{positive}', ...
    'Strength_{negative}', 'Participation Coefficient', 'Betweenness Centrality', ...
    'Eigenvector Centrality', 'Within Module Degree Z-score', 'Clustering Coefficient',...
    'Local Assortativity', 'Distance', 'Local Efficiency'}

metrics_sep = {'Strength','', 'Binarized', 'Degree', ...
    'Strength','(positive)', 'Strength','(negative)', ...
    'Participation', 'Coefficient','Betweenness','Centrality',...
    'Eigenvector', 'Centrality','Within Module', 'Degree Z-score',...
    'Clustering','Coefficient','Local','Assortativity', ...
    'Distance','','Local','Efficiency'}

func_met = {'deg','deg_bin','Spos','Sneg','Ppos','BC','V','Z','clust','AssortPos','distsum','locEffic'}


if isequal(DataSet,'ANT')    
    N=M;
    clear M
    if HOM_size == 5                
        %%% 5x5 HOMS
        for subj = 1:l                        
            M{subj,1} = (N{subj,1}+ N{subj,2})/2; %task1
            M{subj,2} = (N{subj,3}+ N{subj,4})/2; %task2
            M{subj,3} = N{subj,5}; %task3
            M{subj,4} = N{subj,7}; %rest1
            M{subj,5} = N{subj,8}; %rest2
        end
    elseif HOM_size == 7
        %%% 7x7 HOM
        for subj = 1:l                        
            M{subj,1} = N{subj,1}; %task1
            M{subj,2} = N{subj,2};
            M{subj,3} = N{subj,3};
            M{subj,4} = N{subj,4}; %task2
            M{subj,5} = N{subj,5}; %task3
            M{subj,6} = N{subj,7}; %rest1
            M{subj,7} = N{subj,8}; %rest2
        end
    end
    clear N
end


runs = length(M(1,:));
sub_ind = 1:l

clear M_diag0
M_diag0 = cell(l,runs);
for subj = 1:l
    for i = 1:runs        
        M_diag0{subj,i} = M{subj,i} - diag(diag(M{subj,i}));               
    end
end

for subj = 1:l
    for i = 1:runs
        M_abs{subj,i} = abs(M_diag0{subj,i});
        M_bin{subj,i} = zeros(268,268);
        M_bin{subj,i}(M_abs{subj,i}>=0.01) = 1;
    end
end

modu_sign = cell(1,runs);
for subj = 1:l
    for i = 1:runs
        [modu_sign{subj,i}, Q_sign{subj,i}] = community_louvain(M_diag0{subj,i},[],[],'negative_sym');
    end
end


% Measures:

%%% Weighted FC
l = length(M)

weightedFC = cell(l,runs);
for subj = 1:l
    for i = 1:runs         
        tmp = triu(M_diag0{subj,i},1);        
        [ii,jj]=find(tmp==0);        
        tmp(ii(find(ii<jj)),jj(find(ii<jj))) = 0.0001;
        tmp = tmp(tmp~=0);    
        weightedFC{subj,i} = tmp;
    end 
end


%%% Degree_weighted
l = length(M)

deg = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        deg{subj,i} = strengths_und(M_diag0{subj,i});
        deg{subj,i} = deg{subj,i}';
    end
end

% clear std_subj
if isequal(DataSet,'nback')    
    for subj=1:l
        tmp = cell2mat(deg(subj,1:3));               
        std_subj{1}(subj,:) = std(tmp,[],2);
    end
    [~,nodes_nback]=sort(mean(std_subj{1}),'ascend')
elseif isequal(DataSet,'gradCPT')
    for subj=1:l
        tmp = cell2mat(deg(subj,1:3));               
        std_subj{2}(subj,:) = std(tmp,[],2);
    end
    [~,nodes_gradCPT]=sort(mean(std_subj{2}),'ascend')
elseif isequal(DataSet,'ANT')
    for subj=1:l
        tmp = cell2mat(deg(subj,1:5));               
        std_subj{3}(subj,:) = std(tmp,[],2);
    end
    [~,nodes_ANT]=sort(mean(std_subj{3}),'ascend')
elseif isequal(DataSet,'HCP')
    for subj=1:l
        tmp = cell2mat(deg(subj,1:3));               
        std_subj(subj,:) = std(tmp,[],2);
    end
    [~,nodes_HCP]=sort(mean(std_subj),'ascend')
end
%%
nodes = [nodes_nback',nodes_gradCPT',nodes_ANT']
figure
clear std_subj_mean
for d=1:3
    std_subj_mean(:,d) = mean(std_subj{d});
end

[c,p] = corr(std_subj_mean)
%%
%%% Degree_binary
deg_bin = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        deg_bin{subj,i} = strengths_und(M_bin{subj,i});
        deg_bin{subj,i} = deg_bin{subj,i}';
    end
end


%%% Strength
Spos = cell(l,runs);
Sneg = cell(l,runs);
Vpos = cell(l,runs);
Vneg = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        [Spos{subj,i}, Sneg{subj,i}, Vpos{subj,i}, Vneg{subj,i}] = strengths_und_sign(M_diag0{subj,i});
        Spos{subj,i} = Spos{subj,i}';
        Sneg{subj,i} = Sneg{subj,i}';
    end
end

%%% Participation coefficient
Ppos = cell(l,runs);
Pneg = cell(1,runs);
for subj = 1:l
    for i = 1:runs
        [Ppos{subj,i}, Pneg{subj,i}] = participation_coef_sign(M_diag0{subj,i},modu_sign{subj,i});
    end
end


%%% Betweennes centrality
%!!!Takes long!!!

% BC = cell(l,runs);
% for subj = 1:l
%     subj
%     for i = 1:1
%         BC{subj,i} = betweenness_wei(M_diag0{subj,i});
%     end    
%     for i=2:5
%         BC{subj,i} = BC{subj,1};
%     end
% end

load(['/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/BC_',DataSet,'_',num2str(runs),'.mat']);


%%% Eigenvector centrality
V = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        V{subj,i} = eigenvector_centrality_und(M_diag0{subj,i});
    end
end

%%% Within module degree Z-score
Z = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        Z{subj,i} = module_degree_zscore((M_diag0{subj,i}),modu_sign{subj,i},0);
    end
end


%%% Clustering coefficient
clust = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        [clust{subj,i}, ~]= clustering_coef_wu_sign(M_diag0{subj,i});
    end
end


%%% Local Assortativity
AssortPos = cell(l,runs);
AssortNeg = cell(1,runs);
for subj = 1:l
    for i = 1:runs
        [AssortPos{subj,i}, AssortNeg{subj,i}] = local_assortativity_wu_sign(M_diag0{subj,i});
    end
end


%%% Binary Distance
dist = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        dist{subj,i} = distance_bin(M_bin{subj,i});
        dist{subj,i}(~isfinite(dist{subj,i})) = 0;
        distsum{subj,i} = sum(dist{subj,i},2);
    end
end



%%% Local efficiency
%%!!!Takes long!!!

% locEffic = cell(l,runs);
% for subj = 1:l
%     subj
%     for i = 1:5
%         locEffic{subj,i} = efficiency_wei(M_diag0{subj,i},2);
%     end
% end

load (['/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/locEffic_',DataSet,'_',num2str(runs),'.mat']);

% %% Are pos and neg correlated?
% 
% for subj = 1:l
%     for st=1:runs
%         [r(subj,st),p(subj,st)] = corr(Ppos{subj,st},Pneg{subj,st});
%     end
% end

%% Exporting to R
clear data_st hom_entries
%%%% Just Degree
mm=1;
net = eval(func_met{mm});    
clear HOM
for subj = 1:l
    degreeVector = [];
    for i = 1:runs
        degreeVector = [degreeVector,net{subj,i}];
    end        
    degreeVector = abs(degreeVector');
        
    norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);    
    hom = real(atanh(corr(norm_degreeVector',norm_degreeVector','type','Pearson')));
    hom = hom-diag(diag(hom));
    hom(isnan(hom))=0;
    HOM(subj,:,:) = hom;    
    
%     HOM(subj,:,:) = corr(degreeVector');  
    
    
    k = (mm-1)*12+1;
    for i = 1:runs-1
        for j = i+1:runs
            hom_entries(subj,k) = HOM(subj,i,j);
            k = k+1;
        end
    end    
end
        
data_st = [subjname,hom_entries,all_behav];
csvwrite(['/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/',DataSet,'_HOMs_',num2str(runs),'_ZTransformed.csv'],data_st);
%%
%%%%% All nets        
clear data_st hom_entries
for mm = 1:12
    net = eval(func_met{mm});    
    clear HOM
    for subj = 1:l
        degreeVector = [];
        for i = 1:runs
            degreeVector = [degreeVector,net{subj,i}];
        end        
        degreeVector = abs(degreeVector');
        
        norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
        hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
        hom = hom-diag(diag(hom));
        hom(isnan(hom))=0;
        HOM(subj,:,:) = hom;
        
%         HOM(subj,:,:) = corr(degreeVector');
        
        
        k = (mm-1)*(sum(1:runs-3)+2)+1
        tc = triu(squeeze(HOM(subj,1:runs-2,1:runs-2)),1)';
        rtc = squeeze(HOM(subj,1:runs-2,runs-1:runs));
        hom_entries(subj,k:k+sum(1:runs-3)-1) = tc(tc~=0);                        
        hom_entries(subj,k+sum(1:runs-3)) = mean(rtc(:));                            
        hom_entries(subj,k+sum(1:runs-3)+1) = squeeze(HOM(subj,runs-1,runs));
        
        

    end                   
end

data_st = [subjname,hom_entries,all_behav];
csvwrite(['/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/',DataSet,'_HOMs_Allnets_',num2str(runs),'_ZTransformed.csv'],data_st);

%%
%%%%% Just weighted FC        
clear data_st hom_entries

mm=1;
net = eval('weightedFC');                    
clear HOM
for subj = 1:l
    degreeVector = [];
    for i = 1:runs
        degreeVector = [degreeVector,net{subj,i}];
    end        
    degreeVector = abs(degreeVector');
        
    norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
    hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
    hom = hom-diag(diag(hom));
    hom(isnan(hom))=0;
    HOM(subj,:,:) = hom;    
    
%     HOM(subj,:,:) = corr(degreeVector');  
    
    k = (mm-1)*12+1;
    for i = 1:runs-1
        for j = i+1:runs
            hom_entries(subj,k) = HOM(subj,i,j);
            k = k+1;
        end
    end    
end
        


data_st = [subjname,hom_entries,all_behav];
csvwrite(['/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/',DataSet,'_HOMs_WeightedFC_',num2str(runs),'_ZTransformed.csv'],data_st);

%% ANOVA tests TC vs. RC & all HOM entries
close all
mkdir /Users/Mehraveh/Documents/MATLAB/Connectivity/Results/
% anova_type = input('What type of ANOVA test to do? (''TCRC'', ''ALL'')\n','s')
anova_type = 'TCRC';
for mm = 1:1
    net = eval(func_met{mm});
    
    clearvars HOM
    for subj = 1:l
        degreeVector = [];
        for i = 1:runs
            degreeVector = [degreeVector,net{subj,i}];
        end        
        degreeVector = abs(degreeVector)';
        
  
        
        norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
        hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
        hom = hom-diag(diag(hom));
        hom(isnan(hom))=0;
        HOM(subj,:,:) = hom;

%         HOM(subj,:,:) = corr(degreeVector');     
                         
        tc = triu(squeeze(HOM(subj,1:runs-2,1:runs-2)),1)';
        TC(subj,:) = tc(tc~=0);
        RC(subj,1) = squeeze(HOM(subj,runs-1,runs));
        rtc = squeeze(HOM(subj,1:runs-2,runs-1:runs));
        RTC(subj,:) = rtc(:);   
    end           
    
    if isequal(anova_type,'TCRC')
        anovaVector = [TC,RC];
        
        format short e
        % axis square
        if runs == 5
            t = {'<T1,T2>', '<T1,T3>', '<T2,T3>', '<R1,R2>'};
        elseif runs == 7
            t = {'<T1,T2>', '<T1,T3>', '<T1,T4>','<T1,T5>',...
                '<T2,T3>','<T2,T4>','<T2,T5>',...
                '<T3,T4>', '<T3,T5>',...
                '<T4,T5>','<R1,R2>'};
        end
    
    elseif isequal(anova_type,'ALL')
        anovaVector = [TC,RC,RTC];
        
        format short e
        % axis square
        if runs == 5
            t = {'<T1,T2>', '<T1,T3>', '<T2,T3>', '<R1,R2>',...
                '<T1,R1>', '<T2,R1>', '<T3,R1>', ...
                '<T1,R2>','<T2,R2>', '<T3,R2>'};
        elseif runs == 7
            t = {'<T1,T2>', '<T1,T3>', '<T1,T4>','<T1,T5>',...
                '<T2,T3>','<T2,T4>','<T2,T5>',...
                '<T3,T4>', '<T3,T5>',...
                '<T4,T5>',...
                '<R1,R2>',...
                '<T1,R1>', '<T2,R1>', '<T3,R1>','<T4,R1>','<T5,R1>', ...
                '<T1,R2>', '<T2,R2>', '<T3,R2>','<T4,R2>','<T5,R2>'};
        end
    end
    figure
    [pval,~,stats] = anova1(anovaVector)
    fontsize = 15;
    rotation = 0;
    if (runs==7 || isequal(anova_type,'ALL'))
        fontsize = 10;
        rotation = 45;
    end
    set(gca,'xtick',1:length(t),'xticklabel',t,'xticklabelrotation',rotation,'FontSize',fontsize,'FontWeight','bold'); 
    h = findobj(gca,'Type','line')
    set(h,'LineWidth',3)
    
    title({metrics{mm}; ['(p<', num2str(pval,'%.2e'),')']}, 'FontSize', 30)
%     saveas(gcf,['/Users/Mehraveh/Documents/MATLAB/Connectivity/Results/ANOVA_with_ZTransform/anova_',func_met{mm},'_',anova_type,'_',num2str(runs),'_',DataSet,'.png'])
% close all
figure
% [results,means] = multcompare(stats,'CType','lsd')
h=dunnett(stats,[],11)
%     h = findobj(gca,'Type','line')
%     set(h,'LineWidth',3)
    
end

%% Canonical Correlaiton for TASK vs REST, i.e. RTC
close all
figure
for mm = 1:1
    net = eval(func_met{mm});
      
    for subj = 1:l
        degreeVector = [];
        for i = 1:runs
            degreeVector = [degreeVector,net{subj,i}];
        end        
        degreeVector = abs(degreeVector');                        
        
        norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
        hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
        hom = hom-diag(diag(hom));
        hom(isnan(hom))=0;
        HOM(subj,:,:) = hom;

%         HOM(subj,:,:) = corr(degreeVector');     

        
        %%% Claculating the lambda dissimilarity
        Task = degreeVector(1:runs-2,:)';
        Rest = degreeVector(runs-1:runs,:)'; 
        
        tc = triu(squeeze(HOM(subj,1:runs-2,1:runs-2)),1)';
        TC(subj,:) = tc(tc~=0);
        RC(subj,1) = squeeze(HOM(subj,runs-1,runs));
        rtc = squeeze(HOM(subj,1:runs-2,runs-1:runs));
        RTC(subj,:) = rtc(:);   
        
        [A,B,r(subj,:),U,V,stats] = canoncorr(Task,Rest);        
        subplot(subM,subN,subj)
        plot(U(:,1),V(:,1),'.')
        Wilks(subj,:) = stats.Wilks;
        Fstat(subj,:) = stats.F;
        Pvals(subj,:) = stats.p*2;       
        title(['wilk = ', num2str(Wilks(subj,1),'%.2f'),', p<',num2str(Pvals(subj,1),'%.2f')],...
            'fontsize',10)
        set(gca,'fontsize',10)
    end               
end

save(['/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/Wilks_',DataSet,'_',num2str(runs),'.mat'],'Wilks');
%% Wilks visualziation
DataSets = {'nback_5','gradCPT_5','ANT_5','ANT_7','HCP_4'}
clear A
for i=1:length(DataSets)
    d = DataSets{i};
    load (['/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/Wilks_',d,'.mat'])
    A{i} = Wilks(:,1);
end

figure
c = colormap(lines)
DataSets = {'n-back_{5\times 5}','gradCPT_{5\times5}','ANT_{5\times5}','ANT_{7\times7}','HCP_{4\times4}'}
for i=1:length(DataSets)
    H=bar(i,mean(A{i})), hold on
    set(H,'FaceColor',c(i,:),'Edgecolor',c(i,:))
    errorbar(i,mean(A{i}),std(A{i}),'.','Color',c(i,:),'linewidth',3,'CapSize',10)    
end
    set(gca,'FontSize',10)
    yy= mean(A{i}) + std(A{i});
    ylim([0,1.1])
    set(gca,'Xtick',1:length(DataSets),'XtickLabels',DataSets,...
        'fontsize',25,'XtickLabelRotation',45)
    axis('square','xy');
    box('off')
    
    
%     [H,starY]=sigstar(groups,pva

%% How are eigen values and HOM entries related?
% close all
figure
for mm = 1:1
    net = eval(func_met{mm});    
    clearvars HOM
    for subj = 1:l
        degreeVector = [];
        for i = 1:runs
            degreeVector = [degreeVector,net{subj,i}];
        end        
        degreeVector = abs(degreeVector');
        
        norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
        hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
        hom = hom-diag(diag(hom));
        hom(isnan(hom))=0;
        HOM(subj,:,:) = hom;


        HOM(subj,:,:) = corr(degreeVector');  
                           
        %%% Calculating the eigendecomposition

        EigenHOM(subj,:) = sort(eig(squeeze(HOM(subj,:,:))),'descend');
                                      
    end           
    
    [C,P] = corr(EigenHOM,[TC,RTC,RC]);
%     [C,P] = corr(EigenHOM,[mean(TC,2),mean(RTC,2),RC]);
  
mat = C;  
imagesc(mat);
colormap ([hot])
% colorbar
textStrings = num2str(mat(:),'%0.1f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:size(mat,2),1:size(mat,1));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center','fontsize',25);
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
% textColors = repmat(mat(:) < midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color

textColors = repmat(mat(:) < 0.01,1,3).*repmat([1 1 1],size(mat,1)*size(mat,2),1);

set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors

if runs == 5
    X_label = {'TC','TC','TC','RTC','RTC','RTC','RTC','RTC','RTC','RC'};
    Y_label = {'\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5'}  ;  
elseif runs == 7
    X_label = {'TC','TC','TC','TC','TC','TC','TC','TC','TC','TC',...
        'RTC','RTC','RTC','RTC','RTC','RTC','RTC','RTC','RTC','RTC',...
        'RC'};
    Y_label = {'\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5','\lambda_6','\lambda_7'}  ;  
end
% X_label = {'TC_{avg}','RTC_{avg}','RC'};


set(gca,'XTick',1:size(mat,2),...                         %# Change the axes tick marks
        'XTickLabel',[X_label],'XtickLabelRotation',45,...  %#   and tick labels
        'YTick',1:size(mat,1),...
        'YTickLabel',[Y_label],...
        'TickLength',[0 0],'Fontsize',30,'fontweight','bold');

M=mat;
caxis([-1,1])


end

    
%% PCA analysis
clear scores coeff latent explained
clc
close all
for mm = 1:1
    net = eval(func_met{mm});    
    clearvars HOM
    for subj = 1:l
        degreeVector = [];
        for i = 1:runs
            degreeVector = [degreeVector,net{subj,i}];
        end        
        degreeVector = abs(degreeVector');
        
        norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
        hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
        hom = hom-diag(diag(hom));
        hom(isnan(hom))=0;
        HOM(subj,:,:) = hom;
                           
        
        %%% Calculating the eigendecomposition
        [coeff(subj,:,:),scores(subj,:,:),latent(subj,:),tsquared,explained(subj,:)] = pca(degreeVector(1:3,:)');
        
%         [coeff,scores,latent,tsquared,explained] = eig(hom);
% 
%         [Eig,V] = eig(hom);
%         EigenHOM(subj,:) = sort(eig(squeeze(HOM(subj,:,:))),'descend');
                                      
    end           
      
end



%% Marvin's idea for Figure 2
close all
load(['/Users/Mehraveh/Documents/MATLAB/Connectivity/R_results/TC_RC_pvals_',DataSet,'_',num2str(runs),'.mat'])

clear Srho_state_all
%
figure('Color',[1 1 1]);
% set(gca,'Position',[0 0 1 1])
sub_ind = 1:l

for mm = 1:12
    net = eval(func_met{mm});
    
    clearvars HOM
    for subj = 1:l    
        degreeVector = [];
        for i = 1:runs
            degreeVector = [degreeVector,net{subj,i}];
        end        
        degreeVector = abs(degreeVector');
        norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
%     
        hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
        hom = hom-diag(diag(hom));
        hom(isnan(hom))=0;
        HOM(subj,:,:) = hom;
% 
%         HOM(subj,:,:) = corr(degreeVector');    
        
    end
   
    contr_mtrc = zeros(l,1);
    for subj = 1:l
%         HOM0(subj,:,:) = squeeze(HOM(subj,:,:))
        HOM0(subj,:,:) = squeeze(HOM(subj,:,:));
        HOM(subj,:,:) = squeeze(HOM(subj,:,:))+ eye(runs);
    end
    

    
    TC = sum(squeeze(sum(HOM0(:,1:runs-2,1:runs-2),2)),2)/((runs-2)*(runs-3));
    RC = sum(squeeze(sum(HOM0(:,runs-1:runs,runs-1:runs),2)),2)/2;
    RTC = sum(squeeze(sum(HOM0(:,1:runs-2,runs-1:runs),2)),2)/(2*(runs-2));
    contr_mtrc =  TC - RC;    

    HOM_all{mm} = squeeze(sum(HOM));
    HOM_all{mm} = HOM_all{mm} ./length(subjname);
    
    newmap = hot;
    subtightplot (12,10,(mm-1)*10+2,[0.009,0.009],0.005,0.01)
    hh=imshow(HOM_all{mm}, 'InitialMagnification', 800,'Border','tight','Colormap',newmap), hold on
    rectangle('Position',[0.5,0.5,runs,runs],'LineWidth',4)
    ylabel({metrics_sep{mm*2-1};metrics_sep{mm*2}},'FontSize', 15,'fontweight','bold','rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')

    subtightplot (12,10,(mm-1)*10+3,[0.009,0.009],0.005,0.01)
    barData=[TC,RC,RTC];
    c = [[0,176,244];[0,137,39];[121,218,76]]./255;
    for i=1:3
        H=bar(i,mean(barData(:,i))), hold on
        set(H,'FaceColor',c(i,:),'Edgecolor',c(i,:))
        errorbar(i,mean(barData(:,i)),std(barData(:,i)),'.','Color',c(i,:),'linewidth',3,'CapSize',10)    
    end
    set(gca,'FontSize',10)
    yy= mean(barData(:,1)) + std(barData(:,1));
    ylim([0,yy])
    set(gca,'Xtick',1:3)
    axis('square','xy');
    groups = {[1,2]}
%     [H,starY]=sigstar(groups,pval(mm,1),yy);
    ylim([0.01,1.07])
%     ylim([-1,yy])
    
    A(mm,:) = TC - RC;
end

[i,j]=max(mean(A,2))



%% TC vs RC visualizaiton after review
close all
load(['/Users/Mehraveh/Documents/MATLAB/Connectivity/R_results/TC_RC_pvals_',DataSet,'_',num2str(runs),'.mat'])

clear Srho_state_all HOM HOM0 TC RC 
%
figure('Color',[1 1 1]);
% set(gca,'Position',[0 0 1 1])
sub_ind = 1:l

for mm = 1:1
    net = eval(func_met{mm});
    
    clearvars HOM
    for subj = 1:l    
        degreeVector = [];
        for i = 1:runs
            degreeVector = [degreeVector,net{subj,i}];
        end        
        degreeVector = abs(degreeVector)';
        
  
        
        norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
        hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
        hom = hom-diag(diag(hom));
        hom(isnan(hom))=0;
        HOM(subj,:,:) = hom;

                         
        tc = triu(squeeze(HOM(subj,1:runs-2,1:runs-2)),1)';
        TC(subj,:) = tc(tc~=0);
        RC(subj,1) = squeeze(HOM(subj,runs-1,runs));      
        
    end
   
    contr_mtrc = zeros(l,1);
    for subj = 1:l
        HOM0(subj,:,:) = squeeze(HOM(subj,:,:));
        HOM(subj,:,:) = squeeze(HOM(subj,:,:))+ eye(runs);        
    end
    
      
    HOM_all = squeeze(sum(HOM));
    HOM_all = HOM_all./length(subjname);
    
    newmap = hot;
%     figure
    subtightplot (1,1,1,[0.009,0.009],0.005,0.01)
    hh=imshow(HOM_all, 'InitialMagnification', 800,'Border','tight','Colormap',newmap), hold on    
    caxis([min(min(HOM_all))+0.2,max(max(HOM_all(1:runs-2,1:runs-2)-eye(runs-2)))+0.05])
    col=colorbar('fontsize',20)
    
    
   
  

%     rectangle('Position',[0.5,0.5,2,2],`'LineWidth',4)
    
    figure
    subtightplot (1,1,1,[0.4,0.4],0.1,0.1)
    barData=[TC,RC];
    
    c = [repmat([0,176,244],[sum(1:runs-3),1]);[0,137,39]]./255;
    for i=1:size(barData,2)
        H=bar(i,mean(barData(:,i))), hold on
        set(H,'FaceColor',c(i,:),'Edgecolor',c(i,:))
        errorbar(i,mean(barData(:,i)),std(barData(:,i)),'.','Color',c(i,:),'linewidth',3,'CapSize',10)    
    end
    set(gca,'FontSize',40)
    yy= mean(barData(:,1)) + std(barData(:,1));
    ylim([0,1.5])
    set(gca,'Xtick',[])
    set(gca,'Ytick',[0,0.5,1,1.5])
    axis('square','xy');    
    
    box('off')
    
end



%% Generating a stand-alone colorbar
figure
colormap hot
h=colorbar;
set(h,'fontsize',30,'fontweight','bold','Ticks',[0 2.0000e-01 4.0000e-01 6.0000e-01 8.0000e-01 1])
set(h,'location','southoutside')
x1=get(gca,'position');
x=get(h,'Position');
x(3)=0.8;
set(h,'Position',x)
set(gca,'position',x1)

%% Normality
figure
subplot(1,2,1)
hist(TC)
subplot(1,2,2)
hist(RC)
%% Subject 24 with labels of each row and column 5x5 (nback)
label={'Task1','Task2','Task3','Task4','Task5'};
figure('Color',[1 1 1]);
subj=24
title(subj);
newmap = hot;
imshow(squeeze(HOM(subj,:,:)), 'InitialMagnification', 5000,'Colormap',newmap), hold on
rectangle('Position',[0.5,0.5,runs,runs],'LineWidth',4)
h = xlabel (['Subj ' num2str(subj)], 'FontSize', 50,'Fontweight','normal');
% set(h, 'Units', 'Normalized');
% pos = get(h, 'Position');
% set(h, 'Position', pos + [0, -0.2, 0])

%% Are these measures correlated subject wise?

clearvars netCor_subj netPval_subj C P
mm_indice = 1:12;
C=zeros(12,12);
P=zeros(12,12);
for mm1=mm_indice
    mm1
    for mm2=mm_indice
        net1 = eval(func_met{mm1});            
        net2 = eval(func_met{mm2});
        netCor_subj=[];
        netPval_subj=[];
        if (mm1 == 6 ||mm1 == 12 ||mm2 == 6 ||mm2 == 12)
            runs=1
        end
        for subj = 1:l
            for st = 1:runs        
                [netCor_subj(subj,st)] = real(corr(net1{subj,st},net2{subj,st}));
                try
                    [netCor_subj(subj,st), netPval_subj(subj,st)] = corr(net1{subj,st},net2{subj,st});
                end
            end
        end
        C(mm1,mm2) = mean(netCor_subj(:));
        P(mm1,mm2) = mean(netPval_subj(:));
    end
end

labels = {'Strength','Bin Deg','Strength_{+}','Strength_{-}','Part','BC','Eign','Within-Deg Z','Clust','Assort','Dist','Effic'}


C=tril(C,-1);
P=tril(P,-1);
C(C==0)=nan;
P(P==0)=nan;
mat=C;


figure
close all
CC = (othercolor('Spectral4',25))
imagesc(mat);

colormap([CC])
caxis([-1,1])
textStrings = num2str(mat(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:length(mat));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center');
midValue = 0.3 %mean(get(gca,'CLim'));  %# Get the middle value of the color range
% textColors = repmat(mat(:) < midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                         %#   the background color
% tmp = abs(C)>[min(abs(C(P<0.01)))];
% textColors(tmp(:),:)=repmat([0,0,1],size(textColors(tmp(:),:),1),1);

% set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
textWeight = repmat({'normal'},144,1);
textSize = repmat(15,144,1);
tmp = abs(C)>[min(abs(C(P<0.01)))];
textWeight(tmp(:),:)=repmat({'bold'},size(textWeight(tmp(:),:),1),1);
textSize(tmp(:),:)=repmat(18,size(textSize(tmp(:),:),1),1);
set(hStrings,{'FontWeight'},textWeight)

set(hStrings,{'FontSize'},num2cell(textSize,2))
set(gca,...
        'XTick',1:length(mat),...
        'Xticklabel',(metrics),...
        'XTickLabelRotation',45,...
        'YTick',1:length(mat),...
        'YTickLabel',(metrics),...
        'YTickLabelRotation',0,...
        'TickLength',[0 0],'Fontsize',15,'fontweight','bold');
axis('square')

%% Tensoring

for subj = 1:l
    for p=1:268
        for st=1:5
            DegTensor(subj,p,st) = Deg{subj,st}(p);
        end
    end
end

[Ue,Se,sve] = mlsvd(DegTensor);
%
for n = 1:3
    subplot(1,3,n);
    semilogy(sve{n},'.-');
    xlim([1 length(sve{n})])
    grid on
end
%
clc
h = waitbar(0,'tensor')
for P=195:size(Ue{2},2)
waitbar(P/size(Ue{2},2),h)
Uentrunc{1} = Ue{1}(:,1:l);
Uentrunc{2} = Ue{2}(:,1:P);
Uentrunc{3} = Ue{3}(:,1:5);
Sentrunc = Se(1:l, 1:P, 1:5);

Deg_approx=lmlragen(Uentrunc,Sentrunc);
% frob(Deg_approx-DegTensor)/frob(DegTensor)

reduced_Deg = tmprod(DegTensor,Uentrunc{2}',2);

clear deg
for subj = 1:l
    for p=1:268
        for st=1:5
            deg{subj,st}(p,1) = Deg_approx(subj,p,st) ;
        end
    end
end



clc
clearvars task rest tr contr_mtrc Srho_state0

% sub_ind = [2:4,6:24,26:28] % After removing common subjects with ANT

% sub_ind = [2,3,6:8,10:24,26:28] % After removing common subjects with ANT and gradCPT
% deg = eval(func_met{1});

Srho_state = cell (l,1);
for subj = sub_ind
    temp = [deg{subj,1}, deg{subj,2}, deg{subj,3}, deg{subj,4}, deg{subj,5}];
    Srho_state{subj} = corr(abs(temp));
end

%
contr_mtrc = zeros(l,1);
for subj = sub_ind
    Srho_state0{subj} = Srho_state{subj} - eye(5);
    task(subj,1) = sum(sum(Srho_state0{subj}(1:3,1:3)));
    rest(subj,1) = sum(sum(Srho_state0{subj}(4:5,4:5)));
    tr (subj,1)= sum(sum(Srho_state0{subj}(1:3,4:5)));
    contr_mtrc(subj) = task(subj) - rest(subj);
end


% clearvars cor_score pval
for i = 1:6
    [cor_score(i,P),pval(i,P)] = corr(Output(1:l,i),contr_mtrc(1:l), 'type','Pearson');
end
end
%
scores = {'Accuracy','Mean RT', 'Std RT', 'Alerting effect', 'Orienting effect', ' Conflict effect'}
figure
for i=1:6
    subplot(3,2,i)
    plot(cor_score(i,:))
    xlabel ('metric: Task - Rest', 'FontSize',30);
    ylabel([scores{i}], 'FontSize', 30);
end
    figure
plot(pval)

%%
S = randn(2,2,3);
U = {rand(4,2), rand(5,2), rand(6,3)};
T = lmlragen(U,S);
[Ue,Se,sve] = mlsvd(T)
%%
Uentrunc{1} = Ue{1}(:,1:2);
Uentrunc{2} = Ue{2}(:,1:2);
Uentrunc{3} = Ue{3}(:,1:3);
Sentrunc = Se(1:2, 1:2, 1:3);
T=lmlragen(Uentrunc,Sentrunc)

%% Related to behavior??
clear TC RC RTC
for mm = 1:1
    net = eval(func_met{mm});
    
    clearvars HOM
    for subj = 1:l
        degreeVector = [];
        for i = 1:runs
            degreeVector = [degreeVector,net{subj,i}];
        end        
        degreeVector = abs(degreeVector');
        
        norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
        hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
        hom = hom-diag(diag(hom));
        hom(isnan(hom))=0;
        HOM(subj,:,:) = hom;
       
                           
        tc = triu(squeeze(HOM(subj,1:runs-2,1:runs-2)),1)';
        TC(subj,:) = tc(tc~=0);
        RC(subj,1) = squeeze(HOM(subj,runs-1,runs));
        rtc = squeeze(HOM(subj,1:runs-2,runs-1:runs));
        RTC(subj,:) = rtc(:);   
        
        %%% Calculating the eigendecomposition

%         EigenHOM(subj,:) = sort(eig(squeeze(HOM(subj,:,:))),'descend');
                
        
        

    end           
    
    z = behav;
    x = mean(TC,2) 
    y = RC;
    
    
    [c,p] = corr(x,y)
    b = glmfit([x,y],z);
    
    figure('Color',[1 1 1]);
    plot3(x,y,z, '.k', 'markersize',40); hold on    
    yfit = glmval(b,[x,y],'identity');
    plot3(x,y,yfit,'r','LineWidth',6);
    set(gca,'fontsize',30)
    xlabel ('TC', 'FontSize',30);    
    ylabel ('RC', 'FontSize',30);    
    ylabel(behav_label, 'FontSize', 30);
    title(['r = ', num2str(c,'%.5f'), ' p < ' ...
        num2str(p,'%.5f')], 'FontSize', 30);

%     y = behav;
%     x = mean(TC,2) - RC;
%     
%     
%     [c,p] = corr(x,y)
%     [bls,bint,r,rint,stat_lin] = regress(y,[ones(length(x),1) x]);
% 
%     reg = regstats(y,x,'linear');    
%     figure('Color',[1 1 1]);
%     plot(x,y, '.k', 'markersize',40); hold on
%     plot(x,bls(1)+bls(2)*x,'r','LineWidth',6);
%     set(gca,'fontsize',30)
%     xlabel ('metric: TC - RC', 'FontSize',30);    
%     ylabel(behav_label, 'FontSize', 30);
%     title(['r = ', num2str(c,'%.5f'), ' p < ' ...
%         num2str(p,'%.5f')], 'FontSize', 30);

  
end

%% LOCALIZATION
atlas = input('what network atlas to use? (''Shen'', ''Yeo'') \n','s')

if isequal(atlas,'Shen')
    load label.mat
elseif isequal(atlas,'Yeo')
    load Emilyseg268_mapto_Yeo7network_reordered.mat

    label(:,1)=1:268;
    for p=1:268
        tmp = find(net_label(p,:)==1);
        if ~isempty(tmp)
            label(p,2) = tmp
        else
            label(p,2) = 0;
        end
    end
end

clear l_B M_sub
M_sub = cell(1,max(label(:,2)));
for i = 1:max(label(:,2))
    M_sub {i} = label(label(:,2) == i,1)';
end

for i = 1:max(label(:,2))
    l_B(i) = length(find(label(:,2)==i));
end

comb = []
c_selected = []
p_selected = []
h1 = waitbar(0,'exhaustive localization {inner}...')
h2 = waitbar(0,'exhaustive localization {outer}...')

for num = 1:max(label(:,2))
    waitbar(num/max(label(:,2)),h2);
    clear deg_B 
    nets=nchoosek(1:max(label(:,2)),num);
    if num ~= 1
        N_all = M_sub(nets);
    else
        N_all = M_sub(nets)';
    end
    for run = 1:size(nets,1)
        waitbar(run/size(nets,1),h1)
        N=cell2mat(squeeze(N_all(run,:)));
        for subj = 1:l
            degreeVector = [];
            for i = 1:runs                
                net{subj,i} = strengths_und(M_diag0{subj,i}(N,N))';
                degreeVector = [degreeVector,net{subj,i}];
            end        
            degreeVector = abs(degreeVector');

            norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
            hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
            hom = hom-diag(diag(hom));
            hom(isnan(hom))=0;
            HOM(subj,:,:) = hom;


            tc = triu(squeeze(HOM(subj,1:runs-2,1:runs-2)),1)';
            TC(subj,:) = tc(tc~=0);
            RC(subj,1) = squeeze(HOM(subj,runs-1,runs));
            rtc = squeeze(HOM(subj,1:runs-2,runs-1:runs));
            RTC(subj,:) = rtc(:);   


        end           
    
        y = behav;
        x = mean(TC,2) - RC;
        
        [c,p] = corr(x,y, 'type','Pearson');
        if (isequal(DataSet,'nback') || isequal(DataSet,'gradCPT'))
            if c >= 0 && p < 0.05
                comb = [comb;mat2cell(nets(run,:),1,length(nets(run,:)))];
                c_selected = [c_selected;num2cell(c)];
                p_selected = [p_selected;num2cell(p)];
            end
        elseif (isequal(DataSet,'ANT') || isequal(DataSet,'HCP'))
            if c <= 0 && p < 0.05    
                comb = [comb;mat2cell(nets(run,:),1,length(nets(run,:)))];
                c_selected = [c_selected;num2cell(c)];
                p_selected = [p_selected;num2cell(p)];
            end
        end
    end
end


[i,j] = max(abs(cell2mat(c_selected)));
most_significant_combination = comb{j}
most_significant_combination_C = c_selected{j}
most_significant_combination_P = p_selected{j}


j=1;
min_sufficient_combination = comb{j}
min_sufficient_combination_C = c_selected{j}
min_sufficient_combination_P = p_selected{j}
%%

Shen = [];
for i = [4,5,6,7,8,10]
    Shen = union(Shen,M_sub{i});
end

Yeo = [];
for i = [1,2,3,4]
    Yeo = union(Yeo,M_sub{i});
end

%%
length(intersect(Shen,Yeo)*2)/(length(Shen)+length(Yeo))
%% GLM coefficients
clear all
clc
close all 
file = dir(['/Users/Mehraveh/Documents/MATLAB/Connectivity/R_results/GLM*'])
file(3)=[];
file_order=[4,16,14,15,13,3,11,12,7,6,1,5,10,9,8,2]
figure
for data=1:length(file)
    subtightplot (4,4,data,[0.04,0.009],0.03,0.2)
    load(file(file_order(data)).name);
    if ismember(data,[1,6,11,16])
%         RC = coeff(:,size(cioeff,2));
%         TC = mean(coeff(:,1:size(coeff,2)-1),2);
        TC = coeff(:,1);        
        RC = coeff(:,2);
        
    else
%         RC = coeff(size(coeff,1));
%         TC = mean(coeff(1:size(coeff,1)-1),1);
        TC = coeff(1);
        RC = coeff(2);        

    end
        
    barData=[TC,RC];
    c = [[0,176,244];[0,137,39]]./255;
    for i=1:2
        H=bar(i,mean(barData(:,i))), hold on
        set(H,'FaceColor',c(i,:),'Edgecolor',c(i,:));
        errorbar(i,mean(barData(:,i)),std(barData(:,i)),'.','Color',c(i,:),'linewidth',3,'CapSize',10);
    end
    set(gca,'FontSize',10)        
    set(gca,'Xtick',1:2,'XtickLabels',[],...        
        'fontsize',15);
    
    axis('square','xy');
    box('off')
%     grid on
   
end
%% Siyuan Subjects
load('HCP_subj_names_515.csv')
load('HCP_subj_common_names.mat')
load('HCP_subj_common_names_with_behavior_withinScanner_717.mat')
common = intersect(common_717_withinScannerBehav,HCP_subj_names_515);
find(HCP_subj_names_515==setdiff(HCP_subj_names_515,common))
%%
clc
for mm = 1:12
    net = eval(func_met{mm});
    for subj=1:l
        if any(net{subj,4}<0)
            mm
            break
        end
    end
end

%% Significance test
clear h p
close all
for st = 1:runs
%     figure
    for subj = 1:l
%         subplot(6,7,subj)        
%         hist(abs(deg{subj,st}));
        degreeVector = (deg{subj,st}');                
        [h(subj,st),p(subj,st)] = jbtest(degreeVector);        
    end
end
%%
[h,p] = jbtest(behav)
[h,p] = jbtest(mean(TC,2))
[h,p] = jbtest(RC)
%%
[c,p] = corr(mean(TC,2),RC, 'type', 'Pearson')
%% node visualization

% close all
% figure
% c = [[0,176,244];[0,137,39];[121,218,76]]./255;
% H=bar(mean(std_subj)), hold on
% errorbar(mean(std_subj),std(std_subj),'.','linewidth',1,'CapSize',3)    
%% Motion analysis
clear tb
if isequal(DataSet, 'gradCPT')
    load('/Users/Mehraveh/Documents/MATLAB/Connectivity/gradCPT_motion_n25.mat');
    for subj = 1:25   
        tb(subj,1) = str2num(gradCPT_subs{subj})
    end
    subjInd = arrayfun(@(x)find(tb==x,1),subjname);
    
    task_mot = mean_f2f_task(subjInd);
    rest_mot = mean_f2f_rest(subjInd);    

elseif isequal(DataSet, 'ANT')
    load('/Users/Mehraveh/Documents/MATLAB/Connectivity/ANT_motion_n44.mat');
    for subj = 1:44
        tbnums{subj}(strfind(tbnums{subj}, 't')) = []
        tbnums{subj}(strfind(tbnums{subj}, 'b')) = []
        tb(subj,1) = str2num(tbnums{subj})
    end
    
    subjInd = arrayfun(@(x)find(tb==x,1),subjname);
    task_mot = ANT_motion(subjInd,3);
    rest_mot = ANT_motion(subjInd,6);

elseif isequal(DataSet, 'HCP')
    for i = 1:length(subjname)
        i
        tmp1 = load(['/Users/Mehraveh/Documents/MATLAB/Connectivity/Motion/tfMRI_WM_LR/',num2str(common(i)),'_Movement_RelativeRMS_mean.txt']);
        tmp2 = load(['/Users/Mehraveh/Documents/MATLAB/Connectivity/Motion/tfMRI_WM_RL/',num2str(common(i)),'_Movement_RelativeRMS_mean.txt']);
        task_mot(i,1) = (tmp1+tmp2)/2;
        tmp1 = load(['/Users/Mehraveh/Documents/MATLAB/Connectivity/Motion/rfMRI_REST1_LR/',num2str(common(i)),'_movement_rms_mean.txt']);
        tmp2 = load(['/Users/Mehraveh/Documents/MATLAB/Connectivity/Motion/rfMRI_REST1_RL/',num2str(common(i)),'_movement_rms_mean.txt']);
        rest_mot(i,1) = (tmp1+tmp2)/2;
    end
end


net = eval(func_met{1});
    
clearvars HOM
for subj = 1:l    
    degreeVector = [];
    for i = 1:runs
        degreeVector = [degreeVector,net{subj,i}];
    end        
    degreeVector = abs(degreeVector)';



    norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
    hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
    hom = hom-diag(diag(hom));
    hom(isnan(hom))=0;
    HOM(subj,:,:) = hom;


    tc = triu(squeeze(HOM(subj,1:runs-2,1:runs-2)),1)';
    TC(subj,:) = tc(tc~=0);
    RC(subj,1) = squeeze(HOM(subj,runs-1,runs));      

end
    

%%

[c,p] = corr(task_mot,mean(TC,2))
[c,p] = corr(task_mot,RC)

[c,p] = corr(rest_mot,mean(TC,2))
[c,p] = corr(rest_mot,RC)
