

Attp2Ib = table2array(MergedAttp2Ib);
Attp2Ib = str2double(Attp2Ib);
Attp2IbPre = Attp2Ib(:,1);
Attp2IbPost = Attp2Ib(:,2);

histogram(Attp2IbPre,20);
histfit(Attp2IbPre,20,'exponential');
ll= fitdist(Attp2IbPre,'exponential');
hold on
histfit(Attp2IbPost,20,'exponential');

histogram(Attp2IbPost,20);
title('Attp2 Ib Pre vs Post LY354740')
txt = {'n = 1285 AZs'};
text(0.355,95,txt)
legend('Pre-LY354740','Post-LY354740')
xlabel('Probability of release (Pr)');
ylabel('Count')
Attp2IbSort = sortrows(Attp2Ib);

%Attp2 Is
Attp2Is = table2array(MergedAttp2Is);
Attp2Is = str2double(Attp2Is);
Attp2IsPre = Attp2Is(:,1);
Attp2IsPost = Attp2Is(:,2);
histogram(Attp2Is(:,1),20);
hold on
histogram(Attp2Is(:,2),20);
title('Attp2 Is Pre vs Post LY354740')
txt = {'n = 308 AZs'};
text(0.355,95,txt)
legend('Pre-LY354740','Post-LY354740')
xlabel('Probability of release (Pr)');
ylabel('Count')
Attp2IsSort = sortrows(Attp2Is);



