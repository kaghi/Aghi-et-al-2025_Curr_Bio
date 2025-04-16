
x_coords = table2array(AZLocationsIbONLY(:,1));
y_coords = table2array(AZLocationsIbONLY(:,2));
plasticity_label = table2array(AZLocationsIbONLY(:,5));
clr = [0.6350 0.0780 0.1840;0 0.4470 0.7410;0.4660 0.6740 0.1880];
gscatter(x_coords,y_coords,plasticity_label,clr);
legend('Depressing','No change','Facilitating')
title('NMJ 1: 0.2 Hz to 5 Hz 1')

x_coords = table2array(AZLocationsIbONLYS1(:,1));
y_coords = table2array(AZLocationsIbONLYS1(:,2));
plasticity_label = table2array(AZLocationsIbONLYS1(:,5));
clr = [0.6350 0.0780 0.1840;0 0.4470 0.7410;0.4660 0.6740 0.1880];
gscatter(x_coords,y_coords,plasticity_label,clr);
legend('Depressing','No change','Facilitating')
title('NMJ 2: 0.2 Hz to 5 Hz 1')

x_coords = table2array(AZLocationsIbONLYS2(:,1));
y_coords = table2array(AZLocationsIbONLYS2(:,2));
plasticity_label = table2array(AZLocationsIbONLYS2(:,5));
clr = [0.6350 0.0780 0.1840;0 0.4470 0.7410;0.4660 0.6740 0.1880];
gscatter(x_coords,y_coords,plasticity_label,clr);
legend('Depressing','No change','Facilitating')
title('NMJ 3: 0.2 Hz to 5 Hz 1')

x_coords = table2array(AZLocationsIbONLYS3(:,1));
y_coords = table2array(AZLocationsIbONLYS3(:,2));
plasticity_label = table2array(AZLocationsIbONLYS3(:,5));
clr = [0.6350 0.0780 0.1840;0 0.4470 0.7410;0.4660 0.6740 0.1880];
gscatter(x_coords,y_coords,plasticity_label,clr);
legend('Depressing','No change','Facilitating')
title('NMJ 4: 0.2 Hz to 5 Hz 1')

x_coords = table2array(AZLocationsIbONLYS4(:,1));
y_coords = table2array(AZLocationsIbONLYS4(:,2));
plasticity_label = table2array(AZLocationsIbONLYS4(:,5));
clr = [0.6350 0.0780 0.1840;0 0.4470 0.7410;0.4660 0.6740 0.1880];
gscatter(x_coords,y_coords,plasticity_label,clr);
legend('Depressing','No change','Facilitating')
title('NMJ 5: 0.2 Hz to 5 Hz 1')

x_coords = table2array(AZLocationsIbONLYS5(:,1));
y_coords = table2array(AZLocationsIbONLYS5(:,2));
plasticity_label = table2array(AZLocationsIbONLYS5(:,5));
clr = [0.6350 0.0780 0.1840;0 0.4470 0.7410;0.4660 0.6740 0.1880];
gscatter(x_coords,y_coords,plasticity_label,clr);
legend('Depressing','No change','Facilitating')
title('NMJ 6: 0.2 Hz to 5 Hz 1')

x_coords = table2array(AZLocationsIbONLYS6(:,1));
y_coords = table2array(AZLocationsIbONLYS6(:,2));
plasticity_label = table2array(AZLocationsIbONLYS6(:,5));
clr = [0.6350 0.0780 0.1840;0 0.4470 0.7410;0.4660 0.6740 0.1880];
gscatter(x_coords,y_coords,plasticity_label,clr);
legend('Depressing','No change','Facilitating')
title('NMJ 7: 0.2 Hz to 5 Hz 1')

x_coords = table2array(AZLocationsIbONLYS7(:,1));
y_coords = table2array(AZLocationsIbONLYS7(:,2));
plasticity_label = table2array(AZLocationsIbONLYS7(:,5));
clr = [0.6350 0.0780 0.1840;0 0.4470 0.7410;0.4660 0.6740 0.1880];
gscatter(x_coords,y_coords,plasticity_label,clr);
legend('Depressing','No change','Facilitating')
title('NMJ 8: 0.2 Hz to 5 Hz 1')

x_coords = table2array(AZLocationsIbONLYS8(:,1));
y_coords = table2array(AZLocationsIbONLYS8(:,2));
plasticity_label = table2array(AZLocationsIbONLYS8(:,5));
clr = [0.6350 0.0780 0.1840;0 0.4470 0.7410;0.4660 0.6740 0.1880];
gscatter(x_coords,y_coords,plasticity_label,clr);
legend('Depressing','No change','Facilitating')
title('NMJ 9: 0.2 Hz to 5 Hz 1')

x_coords = table2array(AZLocationsIbONLYS9(:,1));
y_coords = table2array(AZLocationsIbONLYS9(:,2));
plasticity_label = table2array(AZLocationsIbONLYS9(:,5));
clr = [0.6350 0.0780 0.1840;0 0.4470 0.7410;0.4660 0.6740 0.1880];
gscatter(x_coords,y_coords,plasticity_label,clr);
legend('Depressing','No change','Facilitating')
title('NMJ 10: 0.2 Hz to 5 Hz 1')


%%%Lets try and see where silent sites are in the 0.2 Hz and 5 Hz
%%%comparison

%NMJ1
x_coords = table2array(AZLocationsIbONLY(:,1));
y_coords = table2array(AZLocationsIbONLY(:,2));

lowFvals = table2array(AZLocationsIbONLY(:,6));
rows = find(lowFvals==0);

silentlowx = x_coords(rows);
silentlowy = y_coords(rows);

highFvals = table2array(AZLocationsIbONLY(:,7));
rows1 = find(highFvals==0);

silenthighx = x_coords(rows1);
silenthighy = y_coords(rows1);

scatter(silentlowx,silentlowy,30,[0.8500 0.3250 0.0980],'filled');
hold on 
scatter(silenthighx,silenthighy,30,[0.4940 0.1840 0.5560],'filled');
legend('0.2 Hz silent sites','5 Hz 1 silent sites')
title('NMJ 1: Silent sites at 0.2 Hz and 5 Hz 1')
%NMJ2
x_coords = table2array(AZLocationsIbONLYS1(:,1));
y_coords = table2array(AZLocationsIbONLYS1(:,2));

lowFvals = table2array(AZLocationsIbONLYS1(:,6));
rows = find(lowFvals==0);

silentlowx = x_coords(rows);
silentlowy = y_coords(rows);

highFvals = table2array(AZLocationsIbONLYS1(:,7));
rows1 = find(highFvals==0);

silenthighx = x_coords(rows1);
silenthighy = y_coords(rows1);

scatter(silentlowx,silentlowy,30,[0.8500 0.3250 0.0980],'filled');
hold on 
scatter(silenthighx,silenthighy,30,[0.4940 0.1840 0.5560],'filled');
legend('0.2 Hz silent sites','5 Hz 1 silent sites')
title('NMJ 2: Silent sites at 0.2 Hz and 5 Hz 1')

%NMJ3
x_coords = table2array(AZLocationsIbONLYS2(:,1));
y_coords = table2array(AZLocationsIbONLYS2(:,2));

lowFvals = table2array(AZLocationsIbONLYS2(:,6));
rows = find(lowFvals==0);

silentlowx = x_coords(rows);
silentlowy = y_coords(rows);

highFvals = table2array(AZLocationsIbONLYS2(:,7));
rows1 = find(highFvals==0);

silenthighx = x_coords(rows1);
silenthighy = y_coords(rows1);

scatter(silentlowx,silentlowy,30,[0.8500 0.3250 0.0980],'filled');
hold on 
scatter(silenthighx,silenthighy,30,[0.4940 0.1840 0.5560],'filled');
legend('0.2 Hz silent sites','5 Hz 1 silent sites')
title('NMJ 3: Silent sites at 0.2 Hz and 5 Hz 1')

%NMJ4
x_coords = table2array(AZLocationsIbONLYS3(:,1));
y_coords = table2array(AZLocationsIbONLYS3(:,2));

lowFvals = table2array(AZLocationsIbONLYS3(:,6));
rows = find(lowFvals==0);

silentlowx = x_coords(rows);
silentlowy = y_coords(rows);

highFvals = table2array(AZLocationsIbONLYS3(:,7));
rows1 = find(highFvals==0);

silenthighx = x_coords(rows1);
silenthighy = y_coords(rows1);

scatter(silentlowx,silentlowy,30,[0.8500 0.3250 0.0980],'filled');
hold on 
scatter(silenthighx,silenthighy,30,[0.4940 0.1840 0.5560],'filled');
legend('0.2 Hz silent sites','5 Hz 1 silent sites')
title('NMJ 4: Silent sites at 0.2 Hz and 5 Hz 1')

%NMJ5
x_coords = table2array(AZLocationsIbONLYS4(:,1));
y_coords = table2array(AZLocationsIbONLYS4(:,2));

lowFvals = table2array(AZLocationsIbONLYS4(:,6));
rows = find(lowFvals==0);

silentlowx = x_coords(rows);
silentlowy = y_coords(rows);

highFvals = table2array(AZLocationsIbONLYS4(:,7));
rows1 = find(highFvals==0);

silenthighx = x_coords(rows1);
silenthighy = y_coords(rows1);

scatter(silentlowx,silentlowy,30,[0.8500 0.3250 0.0980],'filled');
hold on 
scatter(silenthighx,silenthighy,30,[0.4940 0.1840 0.5560],'filled');
legend('0.2 Hz silent sites','5 Hz 1 silent sites')
title('NMJ 5: Silent sites at 0.2 Hz and 5 Hz 1')

%NMJ6
x_coords = table2array(AZLocationsIbONLYS5(:,1));
y_coords = table2array(AZLocationsIbONLYS5(:,2));

lowFvals = table2array(AZLocationsIbONLYS5(:,6));
rows = find(lowFvals==0);

silentlowx = x_coords(rows);
silentlowy = y_coords(rows);

highFvals = table2array(AZLocationsIbONLYS5(:,7));
rows1 = find(highFvals==0);

silenthighx = x_coords(rows1);
silenthighy = y_coords(rows1);

scatter(silentlowx,silentlowy,30,[0.8500 0.3250 0.0980],'filled');
hold on 
scatter(silenthighx,silenthighy,30,[0.4940 0.1840 0.5560],'filled');
legend('0.2 Hz silent sites','5 Hz 1 silent sites')
title('NMJ 6: Silent sites at 0.2 Hz and 5 Hz 1')

%NMJ7
x_coords = table2array(AZLocationsIbONLYS6(:,1));
y_coords = table2array(AZLocationsIbONLYS6(:,2));

lowFvals = table2array(AZLocationsIbONLYS6(:,6));
rows = find(lowFvals==0);

silentlowx = x_coords(rows);
silentlowy = y_coords(rows);

highFvals = table2array(AZLocationsIbONLYS6(:,7));
rows1 = find(highFvals==0);

silenthighx = x_coords(rows1);
silenthighy = y_coords(rows1);

scatter(silentlowx,silentlowy,30,[0.8500 0.3250 0.0980],'filled');
hold on 
scatter(silenthighx,silenthighy,30,[0.4940 0.1840 0.5560],'filled');
legend('0.2 Hz silent sites','5 Hz 1 silent sites')
title('NMJ 7: Silent sites at 0.2 Hz and 5 Hz 1')

%NMJ8
x_coords = table2array(AZLocationsIbONLYS7(:,1));
y_coords = table2array(AZLocationsIbONLYS7(:,2));

lowFvals = table2array(AZLocationsIbONLYS7(:,6));
rows = find(lowFvals==0);

silentlowx = x_coords(rows);
silentlowy = y_coords(rows);

highFvals = table2array(AZLocationsIbONLYS7(:,7));
rows1 = find(highFvals==0);

silenthighx = x_coords(rows1);
silenthighy = y_coords(rows1);

scatter(silentlowx,silentlowy,30,[0.8500 0.3250 0.0980],'filled');
hold on 
scatter(silenthighx,silenthighy,30,[0.4940 0.1840 0.5560],'filled');
legend('0.2 Hz silent sites','5 Hz 1 silent sites')
title('NMJ 8: Silent sites at 0.2 Hz and 5 Hz 1')

%NMJ9
x_coords = table2array(AZLocationsIbONLYS8(:,1));
y_coords = table2array(AZLocationsIbONLYS8(:,2));

lowFvals = table2array(AZLocationsIbONLYS8(:,6));
rows = find(lowFvals==0);

silentlowx = x_coords(rows);
silentlowy = y_coords(rows);

highFvals = table2array(AZLocationsIbONLYS8(:,7));
rows1 = find(highFvals==0);

silenthighx = x_coords(rows1);
silenthighy = y_coords(rows1);

scatter(silentlowx,silentlowy,30,[0.8500 0.3250 0.0980],'filled');
hold on 
scatter(silenthighx,silenthighy,30,[0.4940 0.1840 0.5560],'filled');
legend('0.2 Hz silent sites','5 Hz 1 silent sites')
title('NMJ 9: Silent sites at 0.2 Hz and 5 Hz 1')

%NMJ10
x_coords = table2array(AZLocationsIbONLYS9(:,1));
y_coords = table2array(AZLocationsIbONLYS9(:,2));

lowFvals = table2array(AZLocationsIbONLYS9(:,6));
rows = find(lowFvals==0);

silentlowx = x_coords(rows);
silentlowy = y_coords(rows);

highFvals = table2array(AZLocationsIbONLYS9(:,7));
rows1 = find(highFvals==0);

silenthighx = x_coords(rows1);
silenthighy = y_coords(rows1);

scatter(silentlowx,silentlowy,30,[0.8500 0.3250 0.0980],'filled');
hold on 
scatter(silenthighx,silenthighy,30,[0.4940 0.1840 0.5560],'filled');
legend('0.2 Hz silent sites','5 Hz 1 silent sites')
title('NMJ 10: Silent sites at 0.2 Hz and 5 Hz 1')



