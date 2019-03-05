clear all
clc

addpath(genpath('/home/mehraveh/documents/MATLAB/BCT'))
addpath(genpath('/home/mehraveh/documents/MATLAB/Connectivity'))

load '/home/mehraveh/documents/MATLAB/Connectivity/HCP_subj_common.mat'
% load /home/mehraveh/documents/MATLAB/Parcellation/HCP_data.mat  % Read the files
load /home/mehraveh/documents/MATLAB/Connectivity/HCP_4runs_WM_REST1.mat  % Read the files
load /home/mehraveh/documents/MATLAB/Connectivity/pmat_data.txt
load('/home/mehraveh/documents/MATLAB/Connectivity/HCP_subj_common_names.mat')

DataSet = 'HCP'
% T=M;
% clear M
% 
% for task = 1:length(T)
%     task
%     for lr=1:2
%         for subj = 1:length(T{task,lr})
%             timecourse = T{task,lr}{subj};
%             timecourse = timecourse';
%             norm_timecourse = bsxfun(@rdivide,(timecourse-repmat(mean(timecourse,2),[1,size(timecourse,2)])),std(timecourse,1,2));            
%             temp = real(atanh(corr(norm_timecourse',norm_timecourse')));       
%             temp = temp-diag(diag(temp));
%             temp(isnan(temp))=0;
%             M{subj,task,lr} = temp;
%         end
%     end
% end
% 
% % save /home/mehraveh/documents/MATLAB/Connectivity/HCP_matrices.mat M
% 
% 
% for task = 1:9
%     task
%     LR = HCP_subj_common{task,1};
%     RL = HCP_subj_common{task,2};
%     N(:,task,1) = M(LR,task,1);
%     N(:,task,2) = M(RL,task,2);
%     
% end

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
subjname = common;

% clear M
% M(:,1:2) = squeeze(N(:,3,:));
% M(:,3:4) = squeeze(N(:,9,:));


l = length(M);
runs = length(M(1,:));
sub_ind = 1:l;


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
        M_bin{subj,i}(M_abs{subj,i}>=0.4) = 1;
    end
end

modu_sign = cell(1,runs);
for subj = 1:l
    for i = 1:runs
        [modu_sign{subj,i}, Q_sign{subj,i}] = modularity_louvain_und_sign(M_diag0{subj,i});
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
fprintf('Strength\n')
deg = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        deg{subj,i} = strengths_und(M_diag0{subj,i});
        deg{subj,i} = deg{subj,i}';
    end
end
%%
%%% Degree_binary
fprintf('Degree_binary\n')
deg_bin = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        deg_bin{subj,i} = strengths_und(M_bin{subj,i});
        deg_bin{subj,i} = deg_bin{subj,i}';
    end
end


%%% Strength
fprintf('Pos/Neg Strength\n')
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
fprintf('Participation\n')
Ppos = cell(l,runs);
Pneg = cell(1,runs);
for subj = 1:l
    for i = 1:runs
        [Ppos{subj,i}, Pneg{subj,i}] = participation_coef_sign(M_diag0{subj,i},modu_sign{subj,i});
    end
end



%%% Betweennes centrality
%!!!Takes long!!!
fprintf('BC\n')
% BC = cell(l,runs);
% for subj = 1:l
%     subj
%     for i = 1:runs
%         BC{subj,i} = betweenness_wei(M_diag0{subj,i});
%     end
% end
load (['/home/mehraveh/documents/MATLAB/Connectivity/NetworkMeasures/BC_',DataSet,'_',num2str(runs),'.mat']);


%%% Eigenvector centrality
fprintf('Eigen\n')
V = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        V{subj,i} = eigenvector_centrality_und(M_diag0{subj,i});
    end
end

%%% Within module degree Z-score
fprintf('Z\n')
Z = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        Z{subj,i} = module_degree_zscore((M_diag0{subj,i}),modu_sign{subj,i},0);
    end
end


%%% Clustering coefficient
fprintf('Clust\n')
clust = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        clust {subj,i} = real(clustering_coef_wu(M_diag0{subj,i}));
    end
end


%%% Local Assortativity
fprintf('Assort\n')
AssortPos = cell(l,runs);
AssortNeg = cell(1,runs);
for subj = 1:l
    for i = 1:runs
        [AssortPos{subj,i}, AssortNeg{subj,i}] = local_assortativity_wu_sign(M_diag0{subj,i});
    end
end


%%% Binary Distance
fprintf('Dist\n')
dist = cell(l,runs);
for subj = 1:l
    for i = 1:runs
        dist{subj,i} = distance_bin(M_bin{subj,i});
        dist{subj,i}(~isfinite(dist{subj,i})) = 0;
        distsum{subj,i} = sum(dist{subj,i},2);
    end
end



% !!!Takes long!!!
% fprintf('Effic\n')
% effic = cell(l,runs);
% for subj = 401:600
%     subj
%     for i = 1:runs
%         i
%         effic{subj,i} = efficiency_wei(M_diag0{subj,i},1);
%     end
% end
% 
% save (['/home/mehraveh/documents/MATLAB/Connectivity/NetworkMeasures/locEffic__401_600_',DataSet,'_',num2str(runs),'.mat'],'effic');
% for subj = 1:l
%     for st=1:runs
%         [r(subj,st),p(subj,st)] = corr(Ppos{subj,st},Pneg{subj,st});
%     end
% end

%% Exporting to R

%%%% importing behavior
load('/Users/Mehraveh/Documents/MATLAB/Connectivity/Behaviors/717Subj_Behavior_within_Scanner.mat')
load('/Users/Mehraveh/Documents/MATLAB/Connectivity/Datasets/HCP_subj_common_names_with_behavior_withinScanner_717.mat')
removals = find (ismember(common,setdiff(common,common_717_withinScannerBehav)));
%%
l=length(common)-length(removals)
body = 14 % 0-back
medRT0 = avg_WM_behav(body,B,l);
body = 38 % 2-back
medRT2 = avg_WM_behav(body,B,l);

body = 15 % 0-back
medRT0_target = avg_WM_behav(body,B,l);
body = 39 % 2-back
medRT2_target = avg_WM_behav(body,B,l);

body = 16 % 0-back
medRT0_nontarget = avg_WM_behav(body,B,l);
body = 40 % 2-back
medRT2_nontarget = avg_WM_behav(body,B,l);


all_behav = [medRT0,medRT0_target,medRT0_nontarget,...
    medRT2,medRT2_target,medRT2_nontarget];
all_behav(all_behav==0) = nan;

close all

t = {'median RT_0bk (body)','median RT_0bk (face)','median RT_0bk (place)','median RT_0bk (tool)',...
    'median RT target_0bk (body)','median RT target_0bk (face)', 'median RT target_0bk (place)','median RT target_0bk (tool)',...    
    'median RT nontarget_0bk (body)','median RT nontarget_0bk (face)','median RT nontarget_0bk (place)','median RT nontarget_0bk (tool)',...    
    'median RT_2bk (body)','median RT_2bk (face)','median RT_2bk (place)','median RT_2bk (tool)',...
    'median RT target_2bk (body)','median RT target_2bk (face)', 'median RT target_2bk (place)','median RT target_2bk (tool)',...    
    'median RT nontarget_2bk (body)','median RT nontarget_2bk (face)','median RT nontarget_2bk (place)','median RT nontarget_2bk (tool)'}
            

%% Just Degree
clear data_st hom_entries
mm=1;
net = eval(func_met{mm});    
net(removals,:)=[];
l=length(net);
subjname = common;
subjname(removals)=[];

clear HOM
for subj = 1:l
    degreeVector = [];
    for i = 1:runs
        degreeVector = [degreeVector,net{subj,i}];
    end        
    degreeVector = (degreeVector');
        
    % without z-transform
%     HOM(subj,:,:) = corr(degreeVector');    
 
    % with z-transform
    norm_degreeVector = bsxfun(@rdivide,(degreeVector-repmat(mean(degreeVector,2),[1,size(degreeVector,2)])),std(degreeVector,1,2));            
    hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
    hom = hom-diag(diag(hom));
    hom(isnan(hom))=0;
    HOM(subj,:,:) = hom;    
    
    k = (mm-1)*12+1;
    for i = 1:runs-1
        for j = i+1:runs
            hom_entries(subj,k) = HOM(subj,i,j);
            k = k+1;
        end
    end    
end
        
data_st = [subjname,hom_entries,all_behav];
csvwrite(['/home/mehraveh/documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/',DataSet,'_HOMs_',num2str(runs),'_ZTransformed_noABS.csv'],data_st);%%

%% All nets        
clear data_st hom_entries

subjname = common;
subjname(removals)=[];

for mm = 1:1
    net = eval(func_met{mm});    
    net(removals,:)=[];
    l=length(net);
    clear HOM
    for subj = 1:l
        if subj==242
            continue
        end
        degreeVector = [];
        for i = 1:runs       
            degreeVector = [degreeVector,net{subj,i}];
        end        
        
        degreeVector = abs(degreeVector');
        
        % without z-transform
%         HOM(subj,:,:) = corr(degreeVector');    
        
        % with z-transform
        norm_degreeVector = bsxfun(@rdivide,(degreeVector-repmat(mean(degreeVector,2),[1,size(degreeVector,2)])),std(degreeVector,1,2));            
        hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
        hom = hom-diag(diag(hom));
        hom(isnan(hom))=0;
        HOM(subj,:,:) = hom;
        
        
        k = (mm-1)*(sum(1:runs-3)+2)+1;
        tc = triu(squeeze(HOM(subj,1:runs-2,1:runs-2)),1)';
        rtc = squeeze(HOM(subj,1:runs-2,runs-1:runs));
        hom_entries(subj,k:k+sum(1:runs-3)-1) = tc(tc~=0);                        
        hom_entries(subj,k+sum(1:runs-3)) = mean(rtc(:));                            
        hom_entries(subj,k+sum(1:runs-3)+1) = squeeze(HOM(subj,runs-1,runs));
        
        

    end                   
end

data_st = [subjname,hom_entries,all_behav];
% csvwrite(['/home/mehraveh/documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/',DataSet,'_HOMs_Allnets_',num2str(runs),'_ZTransformed.csv'],data_st);


%% Just weighted FC        
clear data_st hom_entries

mm=1;
net = eval('weightedFC');                    
clear HOM
net(removals,:)=[];
l=length(net);
subjname = common;
subjname(removals)=[];

for subj = 1:l
    degreeVector = [];
    for i = 1:runs
        degreeVector = [degreeVector,net{subj,i}];
    end        
    degreeVector = abs(degreeVector');
        
    norm_degreeVector = bsxfun(@rdivide,(degreeVector-repmat(mean(degreeVector,2),[1,size(degreeVector,2)])),std(degreeVector,1,2)); 
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
csvwrite(['/home/mehraveh/documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/',DataSet,'_HOMs_WeightedFC_',num2str(runs),'_ZTransformed.csv'],data_st);

%% ANOVA tests TC vs. RC & all HOM entries
close all
mkdir /home/mehraveh/documents/MATLAB/Connectivity/Results/
mkdir /home/mehraveh/documents/MATLAB/Connectivity/Results/ANOVA_with_ZTransform
% anova_type = input('What type of ANOVA test to do? (''TCRC'', ''ALL'')\n','s')
anova_type = 'ALL'
% anova_type = 'TCRC'

for mm =1:1
    net = eval(func_met{mm});
    
    clearvars HOM
    for subj = 1:l            
        if sum(isnan(net{subj,2}))>0
            continue
        end
        degreeVector = [];
        for i = 1:runs
            degreeVector = [degreeVector,net{subj,i}];
        end        
        degreeVector = abs(degreeVector');
        
%         HOM(subj,:,:) = corr(degreeVector');       

        norm_degreeVector = bsxfun(@rdivide,(degreeVector-repmat(mean(degreeVector,2),[1,size(degreeVector,2)])),std(degreeVector,1,2));            

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
    
    if isequal(anova_type,'TCRC')
        anovaVector = [TC,RC];
        
        format short e
        % axis square
        if runs == 4
            t = {'<T1,T2>', '<R1,R2>'};
        end
    
    elseif isequal(anova_type,'ALL')
        anovaVector = [TC,RC,RTC];
        format short e
        % axis square
        if runs == 4
            t = {'<T1,T2>', '<R1,R2>', '<T1,R1>', '<T1,R2>',...
                '<T2,R1>', '<T2,R2>';}     
        end
    end
    
    [corrstat,F{mm}] = anova1(anovaVector)
    fontsize = 15;
    rotation = 0;   
   if (isequal(anova_type,'ALL'))
        fontsize = 15;
        rotation = 45;
    end
    set(gca,'xtick',1:length(t),'xticklabel',t,'xticklabelrotation',rotation,'FontSize',fontsize,'FontWeight','bold'); 
    h = findobj(gca,'Type','line')
    set(h,'LineWidth',3)
    
    title({metrics{mm}; ['(p<', num2str(corrstat,'%.2e'),')']}, 'FontSize', 30)
    saveas(gcf,['/home/mehraveh/documents/MATLAB/Connectivity/Results/ANOVA_with_ZTransform/anova_',func_met{mm},'_',anova_type,'_',num2str(runs),'_',DataSet,'.png'])

    
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
        
        norm_degreeVector = bsxfun(@rdivide,(degreeVector-repmat(mean(degreeVector,2),[1,size(degreeVector,2)])),std(degreeVector,1,2));            

        hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
        hom = hom-diag(diag(hom));
        hom(isnan(hom))=0;
        HOM(subj,:,:) = hom;

        
        %%% Claculating the lambda dissimilarity
        Task = degreeVector(1:runs-2,:)';
        Rest = degreeVector(runs-1:runs,:)'; 
        
        tc = triu(squeeze(HOM(subj,1:runs-2,1:runs-2)),1)';
        TC(subj,:) = tc(tc~=0);
        RC(subj,1) = squeeze(HOM(subj,runs-1,runs));
        rtc = squeeze(HOM(subj,1:runs-2,runs-1:runs));
        RTC(subj,:) = rtc(:);   
        
        [A,B,r,U,V,stats] = canoncorr(Task,Rest);        
%         subplot(10,10,subj)
%         plot(U(:,1),V(:,1),'.')
        Wilks(subj,:) = stats.Wilks;
        Fstat(subj,:) = stats.F;
        Pvals(subj,:) = stats.p*2;       
%         title(['wilk = ', num2str(Wilks(subj,1),'%.2f'),', p<',num2str(Pvals(subj,1),'%.2f')],...
%             'fontsize',10)
%         set(gca,'fontsize',10)
    end               
end

save(['/home/mehraveh/documents/MATLAB/Connectivity/NetworkMeasures/Wilks_',DataSet,'_',num2str(runs),'.mat'],'Wilks');
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
        
%         norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
%         hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
%         hom = hom-diag(diag(hom));
%         hom(isnan(hom))=0;
%         HOM(subj,:,:) = hom;


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

    
%%

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
        
%         norm_degreeVector = (degreeVector-mean(degreeVector,2))./std(degreeVector,1,2);
%         hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
%         hom = hom-diag(diag(hom));
%         hom(isnan(hom))=0;
%         HOM(subj,:,:) = hom;


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
        [coeff,scores,latent,tsquared,explained] = pca(HOM(subj,:,:));



%% Marvin's idea for Figure 2
close all
load(['/home/mehraveh/documents/MATLAB/Connectivity/Results/TC_RC_pvals_',DataSet,'_',num2str(runs),'.mat'])

clear Srho_state_all
%
figure('Color',[1 1 1]);
% set(gca,'Position',[0 0 1 1])
sub_ind = 1:l



for mm = 1:11

    net = eval(func_met{mm});
    
    clearvars HOM
    for subj = 1:l    
        degreeVector = [];
        for i = 1:runs
            degreeVector = [degreeVector,net{subj,i}];
        end        
        degreeVector = abs(degreeVector');
        norm_degreeVector = bsxfun(@rdivide,(degreeVector-repmat(mean(degreeVector,2),[1,size(degreeVector,2)])),std(degreeVector,1,2));            
    
        hom = real(atanh(corr(norm_degreeVector',norm_degreeVector')));
        hom = hom-diag(diag(hom));
        hom(isnan(hom))=0;
        HOM(subj,:,:) = hom;

%         HOM(subj,:,:) = corr(degreeVector');    
        
    end
   
    contr_mtrc = zeros(l,1);
    for subj = 1:l
        HOM0(subj,:,:) = squeeze(HOM(subj,:,:));
        HOM(subj,:,:) = squeeze(HOM(subj,:,:)) + eye(runs);
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
        h=errorbar(i,mean(barData(:,i)),std(barData(:,i)),'.','Color',c(i,:),'linewidth',3)  ;
%         h.CapSize=10;
    end
    set(gca,'FontSize',10)
    yy= mean(barData(:,1)) + std(barData(:,1));
%     ylim([0.1,yy])
    set(gca,'Xtick',1:3)
    axis('square','xy');
    groups = {[1,2]}
%     [H,starY]=sigstar(groups,pval(mm,1),yy);
%     ylim([0.01,1.07])
    
    A(mm,:) = TC - RC;
end

[i,j]=max(mean(A,2))



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

%% Related to behavior??

%% IQ - outside scanner
load('/home/mehraveh/documents/MATLAB/Connectivity/HCP_subj_common_names_with_PMAT_715.mat')
removals = find (ismember(common,setdiff(common,common_715_PMAT)));
behav_subj = arrayfun(@(x)find(pmat_data(:,1)==x,1),common_715_PMAT);
behav = pmat_data(behav_subj,2)

clear TC RC RTC
for mm = 1:1
    net = eval(func_met{mm});
    net(removals,:)=[];
    l=length(net);
    clearvars HOM
    for subj = 1:l
        degreeVector = [];
        for i = 1:runs
            degreeVector = [degreeVector,net{subj,i}];
        end        
        degreeVector = abs(degreeVector');
        
        norm_degreeVector = bsxfun(@rdivide,(degreeVector-repmat(mean(degreeVector,2),[1,size(degreeVector,2)])),std(degreeVector,1,2));            
        
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
        
        %%% Calculating the eigendecomposition

%         EigenHOM(subj,:) = sort(eig(squeeze(HOM(subj,:,:))),'descend');
                
        
        

    end           
    
    y = behav;
    x = mean(TC,2)% - RC;
    [c,p] = corr(x,y)
    [bls,bint,r,rint,stat_lin] = regress(y,[ones(length(x),1) x]);

    reg = regstats(y,x,'linear');    
    figure('Color',[1 1 1]);
    plot(x,y, '.k', 'markersize',40); hold on
    plot(x,bls(1)+bls(2)*x,'r','LineWidth',6);
    set(gca,'fontsize',30)
    xlabel ('metric: TC - RC', 'FontSize',30);    
    ylabel('IQ', 'FontSize', 30);
    title(['r = ', num2str(c,'%.5f'), ' p < ' ...
        num2str(p,'%.5f')], 'FontSize', 30);

    
  
end

%% Withinscanner behavior - outside scanner

load('/home/mehraveh/documents/MATLAB/Parcellation/2017/Behavior/717Subj_Behavior_within_Scanner.mat')
load('/home/mehraveh/documents/MATLAB/Connectivity/HCP_subj_common_names_with_behavior_withinScanner_717.mat')
removals = find (ismember(common,setdiff(common,common_717_withinScannerBehav)));

l=length(common)-length(removals)
body = 14 % 0-back
medRT0 = avg_WM_behav(body,B,l);
body = 38 % 2-back
medRT2 = avg_WM_behav(body,B,l);

body = 15 % 0-back
medRT0_target = avg_WM_behav(body,B,l);
body = 39 % 2-back
medRT2_target = avg_WM_behav(body,B,l);

body = 16 % 0-back
medRT0_nontarget = avg_WM_behav(body,B,l);
body = 40 % 2-back
medRT2_nontarget = avg_WM_behav(body,B,l);


all_behav = [medRT0,medRT0_target,medRT0_nontarget,...
    medRT2,medRT2_target,medRT2_nontarget];
all_behav(all_behav==0) = nan;

close all

t = {'median RT_0bk (body)','median RT_0bk (face)','median RT_0bk (place)','median RT_0bk (tool)',...
    'median RT target_0bk (body)','median RT target_0bk (face)', 'median RT target_0bk (place)','median RT target_0bk (tool)',...    
    'median RT nontarget_0bk (body)','median RT nontarget_0bk (face)','median RT nontarget_0bk (place)','median RT nontarget_0bk (tool)',...    
    'median RT_2bk (body)','median RT_2bk (face)','median RT_2bk (place)','median RT_2bk (tool)',...
    'median RT target_2bk (body)','median RT target_2bk (face)', 'median RT target_2bk (place)','median RT target_2bk (tool)',...    
    'median RT nontarget_2bk (body)','median RT nontarget_2bk (face)','median RT nontarget_2bk (place)','median RT nontarget_2bk (tool)'}
            
[corrstat,F] = anova1(all_behav)
fontsize = 15;
rotation = 45; 

set(gca,'xtick',1:length(t),'xticklabel',t,'xticklabelrotation',rotation,'FontSize',fontsize,'FontWeight','bold'); 
h = findobj(gca,'Type','line')
set(h,'LineWidth',3)

% title({metrics{mm}; ['(p<', num2str(corrstat,'%.2e'),')']}, 'FontSize', 30)
% saveas(gcf,['/home/mehraveh/documents/MATLAB/Connectivity/Results/ANOVA_with_ZTransform/anova_',func_met{mm},'_',anova_type,'_',num2str(runs),'_',DataSet,'.png'])
%
% l=length(common)-length(removals)
% body = 2 % 0-back
% medRT0 = avg_WM_behav(body,B,l);
% body = 26 % 2-back
% medRT2 = avg_WM_behav(body,B,l);
% 
% body = 3 % 0-back
% medRT0_target = avg_WM_behav(body,B,l);
% body = 27 % 2-back
% medRT2_target = avg_WM_behav(body,B,l);
% 
% body = 4 % 0-back
% medRT0_nontarget = avg_WM_behav(body,B,l);
% body = 28 % 2-back
% medRT2_nontarget = avg_WM_behav(body,B,l);
% 
% 
% all_behav_ACC = [medRT0,medRT0_target,medRT0_nontarget,...
%     medRT2,medRT2_target,medRT2_nontarget];
% all_behav_ACC(all_behav_ACC==0) = nan;
% 
% close all
% 
% t = {'ACC_0bk (body)','ACC_0bk (face)','ACC_0bk (place)','ACC_0bk (tool)',...
%     'ACC target_0bk (body)','ACC target_0bk (face)', 'ACC target_0bk (place)','ACC target_0bk (tool)',...    
%     'ACC nontarget_0bk (body)','ACC nontarget_0bk (face)','ACC nontarget_0bk (place)','ACC nontarget_0bk (tool)',...    
%     'ACC_2bk (body)','ACC_2bk (face)','ACC_2bk (place)','ACC_2bk (tool)',...
%     'ACC target_2bk (body)','ACC target_2bk (face)', 'ACC target_2bk (place)','ACC target_2bk (tool)',...    
%     'ACC nontarget_2bk (body)','ACC nontarget_2bk (face)','ACC nontarget_2bk (place)','ACC nontarget_2bk (tool)'}
%             
% [corrstat,F] = anova1(all_behav_ACC)
% fontsize = 15;
% rotation = 45; 
% 
% set(gca,'xtick',1:length(t),'xticklabel',t,'xticklabelrotation',rotation,'FontSize',fontsize,'FontWeight','bold'); 
% h = findobj(gca,'Type','line')
% set(h,'LineWidth',3)

% title({metrics{mm}; ['(p<', num2str(corrstat,'%.2e'),')']}, 'FontSize', 30)
% saveas(gcf,['/home/mehraveh/documents/MATLAB/Connectivity/Results/ANOVA_with_ZTransform/anova_',func_met{mm},'_',anova_type,'_',num2str(runs),'_',DataSet,'.png'])

%%

behav=nanmean(all_behav(:,:),2)
% behav = mean([medRT0medRT2],2)
clear TC RC RTC
for mm = 1:1
    net = eval(func_met{mm});
    net(removals,:)=[];
    l=length(net);
    clearvars HOM
    for subj = 1:l
        degreeVector = [];
        for i = 1:runs
            degreeVector = [degreeVector,net{subj,i}];
        end        
        degreeVector = abs(degreeVector');
        
        norm_degreeVector = bsxfun(@rdivide,(degreeVector-repmat(mean(degreeVector,2),[1,size(degreeVector,2)])),std(degreeVector,1,2));            
        
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
    behav(395)=nan;
    y = behav(~isnan(behav));
    x = mean(TC,2) - RC;
    x = x(~isnan(behav));
    [c,p] = corr(x,y)
    [bls,bint,r,rint,stat_lin] = regress(y,[ones(length(x),1) x]);

    reg = regstats(y,x,'linear');    
    figure('Color',[1 1 1]);
    plot(x,y, '.k', 'markersize',40); hold on
    plot(x,bls(1)+bls(2)*x,'r','LineWidth',6);
    set(gca,'fontsize',30)
    xlabel ('metric: RC ', 'FontSize',30);    
    ylabel('Median RT (mean_{0-bk,2-bk})', 'FontSize', 30);
    title(['r = ', num2str(c,'%.5f'), ' p < ' ...
        num2str(p,'%.5f')], 'FontSize', 30);

    
  

end

%% Are these measures correlated subject wise?

clearvars netCor_subj netPval_subj C P
C=zeros(12,12)
P=zeros(12,12)
for mm1=1:11
    mm1
    for mm2=1:11
        net1 = eval(func_met{mm1});            
        net2 = eval(func_met{mm2});
        netCor_subj=[];
        netPval_subj=[];
        for subj = 1:l
            for st = 1:runs              
                [netCor_subj(subj,st)] = real(corr(net1{subj,st},net2{subj,st}));
                try
                    [netCor_subj(subj,st), netPval_subj(subj,st)] = corr(net1{subj,st},net2{subj,st});
                end
            end
        end
        C(mm1,mm2) = nanmean(netCor_subj(:));
        P(mm1,mm2) = nanmean(netPval_subj(:));
    end
end


labels = {'Strength','Bin Deg','Strength_{+}','Strength_{-}','Part','BC','Eign','Within-Deg Z','Clust','Assort','Dist','Effic'}


C=tril(C,-1);
P=tril(P,-1);
C(C==0)=nan;
P(P==0)=nan;
mat=C;

%%

close all
figure
CC = (othercolor('Spectral4',25))
imagesc(mat);

colormap([CC])

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
textSize(tmp(:),:)=repmat(16,size(textSize(tmp(:),:),1),1);
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
caxis([-1,1])
axis('square')
%% LOCALIZATION
atlas = input('what network atlas to use? (''Shen'', ''Yeo'') \n','s')

if isequal(atlas,'Shen')
    load /home/mehraveh/documents/MATLAB/label.mat
elseif isequal(atlas,'Yeo')
    load /home/mehraveh/documents/MATLAB/Emilyseg268_mapto_Yeo7network_reordered.mat

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
l=length(M)
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

            norm_degreeVector = bsxfun(@rdivide,(degreeVector-repmat(mean(degreeVector,2),[1,size(degreeVector,2)])),std(degreeVector,1,2));            
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

        
        
        TC(removals,:)=[];
        RC(removals,:)=[];

        behav(395)=nan;
        y = behav(~isnan(behav));
        x = mean(TC,2) - RC;
        x = x(~isnan(behav));
        
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
%% Motion analysis
clear task_mot rest_mot
addpath(genpath('/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/movement_RMS_means/'))
load('/home/mehraveh/documents/MATLAB/Connectivity/HCP_subj_common_low_motion_Siyuan_names.mat')
common = intersect(subjname,subjname);

for i = 1:length(common)
    i
    tmp1 = load(['/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/movement_RMS_means/tfMRI_WM_LR/',num2str(common(i)),'_Movement_RelativeRMS_mean.txt']);
    tmp2 = load(['/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/movement_RMS_means/tfMRI_WM_RL/',num2str(common(i)),'_Movement_RelativeRMS_mean.txt']);
    task_mot(i,1) = tmp1;%(tmp1+tmp2)/2;
    tmp1 = load(['/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/movement_RMS_means/rfMRI_REST2_LR/',num2str(common(i)),'_Movement_RelativeRMS_mean.txt']);
    tmp2 = load(['/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/movement_RMS_means/rfMRI_REST2_RL/',num2str(common(i)),'_Movement_RelativeRMS_mean.txt']);
    rest_mot(i,1) = tmp2;%(tmp1+tmp2)/2;
end

%%
subjInd = arrayfun(@(x)find(subjname==x,1),common);

[c,p] = corr(task_mot,TC(subjInd))
[c,p] = corr(task_mot,RC(subjInd))

[c,p] = corr(rest_mot,TC(subjInd))
[c,p] = scatter(rest_mot,RC(subjInd))
