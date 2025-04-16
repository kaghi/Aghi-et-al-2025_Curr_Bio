%Phil Mendonca - Jan 2024. Analyze spatial distribution of AZs that display
%distinct plasticity mechanisms. Script written for Krisha's data analysis

clear all
close all
clc

data_all=cell(1);

%Load data
for i=1:10 %number of sheets
data_all{i}=readmatrix('AZ_Locations_Ib_ONLY.xls','Sheet',i);
end

Dunn=[];
Silho=[];
Nclust=[];
Sclust=[];
Aclust=[];
Dclust=[];
Oclust_CH_Kmeans_all=[];
Oclust_DB_Kmeans_all=[];
Oclust_S_Kmeans_all=[];
Oclust_CH_Link_all=[];
Oclust_DB_Link_all=[];
Oclust_S_Link_all=[];

%Cycle through threshold distance to find out which is the optimal cluster size
th_range=200:10:3000; %threshold range for each cluster size. This value is used for the hiearchical cluster analysis

for t=1:length(th_range)
th_dist=th_range(t);

for i=1:length(data_all)  %for each NMJ
data=data_all{i};
    x=data(:,1);
    y=data(:,2);
    coord=[x y];
    plast_abs=data(:,4);
    plast_bin=data(:,5);
    plast02=data(:,6);
    plast5=data(:,7);

    plast=plast_abs;

Dist=pdist(coord); % Distance between all datapoints
Sq=squareform(Dist);

% Hierchical cluster analysis
Z = linkage(coord,'ward');
T = cluster(Z,'cutoff',th_dist,'Criterion','distance'); %Index of cluster for each event
n_clusters=length(unique(T)); %Number of cluster

% Dunn's index
Dunn(i,t)=dunns(n_clusters, Sq, T);

%Silhouete index
Silho(i,t)=mean(silhouette(coord,T));

% Number of clusters for each th_dist
Nclust(i,t)=n_clusters;

% Size of cluster for each th_dist
size_cluster_tmp=[];
for n=1:n_clusters
size_cluster_tmp(n)=length(find(T==n)); % array with the number of points per cluster
end
Sclust(i,t)=mean(size_cluster_tmp);

% Centroid of all clusters
centroid=[];
    for j=1:n_clusters
    idx=[];
    idx=find(T==j);

    cent_x=mean(x(idx)); cent_y=mean(y(idx));
    centroid(j,1)=cent_x; centroid(j,2)=cent_y;
    end

% Area and boundary of each cluster
k_c=[];area_c=[];x_events=[];y_events=[];
for j=1:n_clusters
    idx=[]; c=[];
    idx=find(T==j);
 
    x_tmp=x(idx,:);y_tmp=y(idx,:);
    x_events{j}=x_tmp; y_events{j}=y_tmp;
    
    if size(x_tmp,1)<=2%if there are only 2 events or less, use boundary instead of hull
    [k_c{j},area_c{j}] = boundary(x_tmp,y_tmp);
    else
    [k_c{j},area_c{j}] = convhull(x_tmp,y_tmp);  
    end
    area_c{j}=area_c{j}/10^6; %converting to um2 so the values don't look absurd
end

% Average cluster area
area_c_array=cell2mat(area_c);
Aclust(i,t)=mean(area_c_array);

% Average cluster density (events per area)
area_density=size_cluster_tmp./area_c_array;
Dclust(i,t)=mean(area_density);
end
end %------------------------------------------------------------------

%Plots -------------------------------------------------------------
%Plot data - Dunn's index on Y axis, Distance on X axis, use Shade Plot to show all cells
figure(1)
for i=1:length(data_all)  %for each NMJ
    plot(th_range,Dunn(i,:));
hold on
end
title('Dunns Index for each cell')

figure(2)
for i=1:length(data_all)  %for each NMJ
    plot(th_range,Silho(i,:));
hold on
end
title('Dunns Index for each cell')

figure(3)
stdshade(Dunn,0.3,'b',th_range,1);hold on
ylabel('Dunn index')
xlabel('Cluster size (nm)')
title('Average Dunns Index')

figure(4)
stdshade(Silho,0.3,'r',th_range,1);hold on
ylabel('Silhouette Index')
xlabel('Cluster size (nm)')
title('Average Silhouette Index')

figure(5)
stdshade(Nclust,0.3,'k',th_range,1);hold on
title('Number of clusters')

figure(6)
stdshade(Sclust,0.3,'g',th_range,1);hold on
title('Average events per cluster')

figure(7)
stdshade(Aclust,0.3,'b',th_range,1);hold on
title('Average cluster area (um2)')

figure(8)
stdshade(Dclust(:,3:end),0.3,'r',th_range(:,3:end),1);hold on %some initial columns have inf (area=0, clusters made of just two events, so skip it)
title('Average cluster density ')
ylabel('AZs/um2')
xlabel('Cluster size (nm)')

cluster_min=700;
cluster_max=900;

figure(9)% summary figure with most important parts
subplot(3,1,1)
stdshade(Dunn,0.3,'b',th_range,1);hold on
ylabel('Dunn index')
xlabel('Cluster size (nm)')

    % Draw range where all 3 metrics indicate stable clusters
    x = [cluster_min cluster_max cluster_max cluster_min];
    y = [0 0 0.2 0.2];
    fill(x,y,'y','LineStyle','none','FaceAlpha',0.5)

xlim([th_range(1)-10 th_range(end)])
ylim([0.05 0.2])

subplot(3,1,2)
stdshade(Silho,0.3,'r',th_range,1);hold on
ylabel('Silhouette Index')
xlabel('Cluster size (nm)')

    % Draw range where all 3 metrics indicate stable clusters
    x = [cluster_min cluster_max cluster_max cluster_min];
    y = [0 0 1 1];
    fill(x,y,'y','LineStyle','none','FaceAlpha',0.5)

xlim([th_range(1)-10 th_range(end)])
ylim([0.5 0.8])

subplot(3,1,3)
stdshade(Dclust(:,3:end),0.3,'g',th_range(:,3:end),1);hold on
ylabel('Cluster densitiy (AZs/um2)')
xlabel('Cluster size (nm)')

    % Draw range where all 3 metrics indicate stable clusters
    x = [cluster_min cluster_max cluster_max cluster_min];
    y = [0 0 2000 2000];
    fill(x,y,'y','LineStyle','none','FaceAlpha',0.5)

xlim([th_range(1)-10 th_range(end)])
tmp=Dclust(:,3:end);
ylim([min(Dclust(:)) 1800])



%% 2 Plot example bouton with clustered AZs

clear all
close all
clc

cluster_size=850; %nm

data_all=cell(1);

%Load data
for i=1:10 %number of shtees
data_all{i}=readmatrix('AZ_Locations_Ib_ONLY.xls','Sheet',i);
end

th_dist=cluster_size;
for i=1:length(data_all) %for each NMJ
data=data_all{i};
    x=data(:,1);
    y=data(:,2);
    coord=[x y];
    plast_abs=data(:,4);
    plast_bin=data(:,5);
    plast02=data(:,6);
    plast5=data(:,7);

    plast=plast5;

Dist=pdist(coord);
Sq=squareform(Dist);

Z = linkage(coord,'ward');
T = cluster(Z,'cutoff',th_dist,'Criterion','distance'); %Index of cluster for each position
n_clusters=length(unique(T)); %Number of cluster

% Force a pattern (debugging, check if the script works)
% T=sort(T);
% plast02=sort(plast02);
% plast5=sort(plast5);
% plast=plast5;

% Get centroid of all clusters
centroid=[];
    for j=1:n_clusters
    idx=[];
    idx=find(T==j);

    cent_x=mean(x(idx)); cent_y=mean(y(idx));
    centroid(j,1)=cent_x; centroid(j,2)=cent_y;
    end

%Area and boundary of each cluster
k_c=[];area_c=[];x_events=[];y_events=[];
for j=1:n_clusters
    idx=[]; c=[];
    idx=find(T==j);
 
    x_tmp=x(idx,:);y_tmp=y(idx,:);
    x_events{j}=x_tmp; y_events{j}=y_tmp;
    
    if size(x_tmp,1)<=2%if there are only 2 events or less, use boundary instead of hull
    [k_c{j},area_c{j}] = boundary(x_tmp,y_tmp);
    else
    [k_c{j},area_c{j}] = convhull(x_tmp,y_tmp);  
    end
end

figure(i)
%Plot data - use Pr heatmaped to plot location of release

tp=0.4;
scatter(x,y,[],plast,'filled','MarkerFaceAlpha',tp);hold on
colormap(jet) %

plot(centroid(:,1),centroid(:,2),'k.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10);hold on

%Plot clusters
for j=1:n_clusters
x_temp=x_events{j}; y_temp=y_events{j};
k_temp=k_c{j};
plot(x_temp(k_temp),y_temp(k_temp),'k');hold on
end

xlabel('X position'); ylabel('Y position')
title('Release Map')
colorbar;
axis square
end


%% 3a Cluster data analysis (spatial pattern)

clear all
close all
clc

MC_interations=100;

data_all=cell(1);
S=struct;

%Load data
for i=1:10 %number of shtees
data_all{i}=readmatrix('AZ_Locations_Ib_ONLY.xls','Sheet',i);
S(i).table=data_all{i};
end

th_dist=850;
for i=1:length(data_all) %for each NMJ
data=data_all{i};
    x=data(:,1);
    y=data(:,2);
    coord=[x y];
    plast_abs=data(:,4);
    plast_bin=data(:,5);
    plast02=data(:,6);
    plast5=data(:,7);

    plast=plast5;

Dist=pdist(coord);
Sq=squareform(Dist);

Z = linkage(coord,'ward');
T = cluster(Z,'cutoff',th_dist,'Criterion','distance'); %Index of cluster for each position
n_clusters=length(unique(T)); %Number of cluster

% Force a pattern (debugging, check the script works)
% T=sort(T);
% plast02=sort(plast02);
% plast5=sort(plast5);

% Get centroid of all clusters
centroid=[];
    for j=1:n_clusters
    idx=[];
    idx=find(T==j);

    cent_x=mean(x(idx)); cent_y=mean(y(idx));
    centroid(j,1)=cent_x; centroid(j,2)=cent_y;
    end

%Area and boundary of each cluster
k_c=[];area_c=[];x_events=[];y_events=[];
for j=1:n_clusters
    idx=[]; c=[];
    idx=find(T==j);
 
    x_tmp=x(idx,:);y_tmp=y(idx,:);
    x_events{j}=x_tmp; y_events{j}=y_tmp;
    
    if size(x_tmp,1)<=2%if there are only 2 events or less, use boundary instead of hull
    [k_c{j},area_c{j}] = boundary(x_tmp,y_tmp);
    else
    [k_c{j},area_c{j}] = convhull(x_tmp,y_tmp);  
    end
end

% Shuffled position
coord_rand=[];
rand_table=[];
for k=1:MC_interations %define number of iterations
r_idx=randperm(size(coord,1));
r_idx=r_idx';
%tmp=coord;
%coord_rand{k}=tmp([r_idx],:);
rand_table(:,k)=r_idx; %This is the index of the shuffled position, use it to extract Pr
end

% Extract values of each cluster, like mean Pr, Std Pr, Difference of 
% cluster Pr to global pr, difference between Pr 0.2 and 5hz

mean_Pr02=[];
std_Pr02=[];
mean_Pr5=[];
std_Pr5=[];
mean_Pr_abs=[];
std_Pr_abs=[];
fac_frac=[];

for n=1:n_clusters %for the number of cluster for each cell
Pr02_mean=plast02(find(T==n));
mean_Pr02(n)=mean(Pr02_mean);
std_Pr02(n)=std(Pr02_mean);

Pr5=plast5(find(T==n));
mean_Pr5(n)=mean(Pr5);
std_Pr5(n)=std(Pr5);

Pr_abs=plast_abs(find(T==n));
mean_Pr_abs(n)=mean(Pr_abs);
std_Pr_abs(n)=std(Pr_abs);
fac_frac(n)=length(find(Pr_abs>0))/length(Pr_abs);% facilitating fraction;
end


mean_Pr02_rand_all=[];
std_Pr02_rand_all=[];
mean_Pr5_rand_all=[];
std_Pr5_rand_all=[];
mean_Pr_abs_rand_all=[];
std_Pr_abs_rand_all=[];
fac_frac_rand_all=[];
%Extract values for suffled versions
for n=1:n_clusters 
    for k=1:size(rand_table,2)% number of MC iterations
    idx=rand_table(:,k);
        
    plast02_rand=plast02(idx);% this is where the index gets shuffled
    Pr02_mean_rand=plast02_rand(find(T==n));
    mean_Pr02_rand_all(n,k)=mean(Pr02_mean_rand);
    std_Pr02_rand_all(n,k)=std(Pr02_mean_rand);

    plast5_rand=plast5(idx);% this is where the index gets shuffled
    Pr5_rand=plast5_rand(find(T==n));
    mean_Pr5_rand_all(n,k)=mean(Pr5_rand);
    std_Pr5_rand_all(n,k)=std(Pr5_rand);

    plast_abs_rand=plast_abs(idx);% this is where the index gets shuffled
    Pr_abs_rand=plast_abs_rand(find(T==n));
    mean_Pr_abs_rand_all(n,k)=mean(Pr_abs_rand);
    std_Pr_abs_rand_all(n,k)=std(Pr_abs_rand);
    fac_frac_rand_all(n,k)=length(find(Pr_abs_rand>0))/length(Pr_abs_rand);% facilitating fraction;
    end
end
mean_Pr02_rand=[];
std_Pr02_rand=[];
mean_Pr5_rand=[];
std_Pr5_rand=[];
mean_Pr_abs_rand=[];
std_Pr_abs_rand=[];

mean_Pr02_rand=mean(mean_Pr02_rand_all,2);
std_Pr02_rand=mean(std_Pr02_rand_all,2);

mean_Pr5_rand=mean(mean_Pr5_rand_all,2);
std_Pr5_rand=mean(std_Pr5_rand_all,2);

mean_Pr_abs_rand=mean(mean_Pr_abs_rand_all,2);
std_Pr_abs_rand=mean(std_Pr_abs_rand_all,2);
fac_frac_rand=mean(fac_frac_rand_all,2);

%Pass values to struct
S(i).n_clusters=n_clusters;
S(i).i_clusters=T;
S(i).coord_rand=rand_table;

S(i).mean_Pr02=mean_Pr02;
S(i).std_Pr02=std_Pr02;
S(i).mean_Pr5=mean_Pr5;
S(i).std_Pr5=std_Pr5;
S(i).mean_Pr_abs=mean_Pr_abs;
S(i).std_Pr_abs=std_Pr_abs;
S(i).fac_frac=fac_frac;

S(i).mean_Pr02_rand_all=mean_Pr02_rand_all;
S(i).std_Pr02_rand_all=std_Pr02_rand_all;
S(i).mean_Pr5_rand_all=mean_Pr5_rand_all;
S(i).std_Pr5_rand_all=std_Pr5_rand_all;
S(i).mean_Pr_abs_rand_all=mean_Pr_abs_rand_all;
S(i).std_Pr_abs_rand_all=std_Pr_abs_rand_all;
S(i).fac_frac_rand_all=fac_frac_rand_all;

S(i).mean_Pr02_rand=mean_Pr02_rand;
S(i).std_Pr02_rand=std_Pr02_rand;
S(i).mean_Pr5_rand=mean_Pr5_rand;
S(i).std_Pr5_rand=std_Pr5_rand;
S(i).mean_Pr_abs_rand=mean_Pr_abs_rand;
S(i).std_Pr_abs_rand=std_Pr_abs_rand;
S(i).fac_frac_rand=fac_frac_rand;
end

%% 3b Plot data

close all
clearvars -except S

for c=3%1:10 %choose cell to plot
%c=3; %cell

Pr02_mean=S(c).mean_Pr02;
Pr02_mean_rand=S(c).mean_Pr02_rand;

Pr5_mean=S(c).mean_Pr5;
Pr5_mean_rand=S(c).mean_Pr5_rand;

Pr_abs_mean=S(c).mean_Pr_abs;
Pr_abs_mean_rand=S(c).mean_Pr_abs_rand;

Pr02_std=S(c).std_Pr02;
Pr02_std_rand=S(c).std_Pr02_rand;

Pr5_std=S(c).std_Pr5;
Pr5_std_rand=S(c).std_Pr5_rand;

Pr_abs_std=S(c).std_Pr_abs;
Pr_abs_std_rand=S(c).std_Pr_abs_rand;

fac_frac=S(c).fac_frac;
fac_frac_rand=S(c).fac_frac_rand;

f1=figure(1);%------------------------------------------------------

subplot(3,3,1)
for i=1:size(Pr02_mean,2) %for each cluster
plot([1 2], [Pr02_mean(i) Pr02_mean_rand(i)],'b--o'); hold on
end
title('Mean Pr at 0.2 Hz')
xticklabels({'Original data','MC version'})

subplot(3,3,2)
for i=1:size(Pr5_mean,2) %for each cluster
plot([1 2], [Pr5_mean(i) Pr5_mean_rand(i)],'b--o'); hold on
end
title('Mean Pr at 5 Hz')
xticklabels({'Original data','MC version'})

subplot(3,3,3)
for i=1:size(Pr_abs_mean,2) %for each cluster
plot([1 2], [Pr_abs_mean(i) Pr_abs_mean_rand(i)],'b--o'); hold on
end
title('Mean Pr difference')
xticklabels({'Original data','MC version'})

subplot(3,3,4)
for i=1:size(Pr02_std,2) %for each cluster
plot([1 2], [Pr02_std(i) Pr02_std_rand(i)],'r--o'); hold on
end
title('Std Pr at 0.2 Hz')
xticklabels({'Original data','MC version'})

subplot(3,3,5)
for i=1:size(Pr5_std,2) %for each cluster
plot([1 2], [Pr5_std(i) Pr5_std_rand(i)],'r--o'); hold on
end
title('Std Pr at 5 Hz')

subplot(3,3,6)
for i=1:size(Pr_abs_std,2) %for each cluster
plot([1 2], [Pr_abs_std(i) Pr_abs_std_rand(i)],'r--o'); hold on
end
title('Std Pr difference')
xticklabels({'Original data','MC version'})

subplot(3,3,7)
for i=1:size(fac_frac,2) %for each cluster
plot([1 2], [fac_frac(i) fac_frac_rand(i)],'g--o'); hold on
end
title('Fraction of Facilitating Synapses')
xticklabels({'Original data','MC version'})


figure(2)
for i=1:size(Pr02_mean,2) %for each cluster
plot([1 2], [Pr02_mean(i) Pr5_mean(i)],'b--o'); hold on
end
xticklabels({'Mean Pr at 0.2 Hz','Mean Pr at 5 Hz'})
end































