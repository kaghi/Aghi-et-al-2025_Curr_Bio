
%Compare 0.2 Hz evoked Pr distribution
%Ib Pr
Ib_Pr = cell(10,1);
for nn = 1:10
Ib_Pr{nn,1} = ExperimentSet_Reduced(3).Grouped_Data_Reduced(nn).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(1).Recording.Evoked_Pr 
end

Conc_Ib_Pr = Ib_Pr{1,1};
for uu = 2:10
AddPr = Ib_Pr{uu,1}
Conc_Ib_Pr = [Conc_Ib_Pr AddPr]
end

Ib_Pr_HF_1 = cell(10,1)
for nn = 1:10
Ib_Pr_HF_1{nn,1} = ExperimentSet_Reduced(3).Grouped_Data_Reduced(nn).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(4).Recording(1).Evoked_Pr 
end

Conc_Ib_Pr_HF_1 = Ib_Pr_HF_1{1,1}
for uu = 2:10
AddPrHf = Ib_Pr_HF_1{uu,1}
Conc_Ib_Pr_HF_1 = [Conc_Ib_Pr_HF_1 AddPrHf]
end
%Compare 0.2Hz and 5Hz first train in Ib

scatter(Conc_Ib_Pr,Conc_Ib_Pr_HF_1,'.');
xlabel('0.2 Hz Pr')
ylabel('5 Hz Pr Train 1')
title('Comparison between conditions in Ib input')
xlim([0 1])
ylim([0 1])
x = [0;0.05;1];
y = x
hold on
plot(x,y);


%Is Pr

Is_Pr = cell(10,1);
for nn = 1:5
Is_Pr{nn,1} = ExperimentSet_Reduced(3).Grouped_Data_Reduced(nn).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(1).Recording.Evoked_Pr 
end
for nn = 7:10
Is_Pr{nn,1} = ExperimentSet_Reduced(3).Grouped_Data_Reduced(nn).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(1).Recording.Evoked_Pr 
end

Conc_Is_Pr = Is_Pr{1,1}
for uu = 2:10
AddPr = Is_Pr{uu,1}
Conc_Is_Pr = [Conc_Is_Pr AddPr]
end

Is_Pr_HF_1 = cell(10,1)
for nn = 1:5
Is_Pr_HF_1{nn,1} = ExperimentSet_Reduced(3).Grouped_Data_Reduced(nn).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(4).Recording(1).Evoked_Pr 
end
for nn = 7:10
Is_Pr_HF_1{nn,1} = ExperimentSet_Reduced(3).Grouped_Data_Reduced(nn).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(4).Recording(1).Evoked_Pr 
end

Conc_Is_Pr_HF_1 = Is_Pr_HF_1{1,1}
for uu = 2:10
AddPrHf = Is_Pr_HF_1{uu,1}
Conc_Is_Pr_HF_1 = [Conc_Is_Pr_HF_1 AddPrHf]
end

%Compare 0.2Hz and 5Hz first train in Is
ComparePr = [Conc_Is_Pr; Conc_Is_Pr_HF_1];
ComparePr = ComparePr.';
scatter(Conc_Is_Pr,Conc_Is_Pr_HF_1,'.');
xlabel('0.2 Hz Pr')
ylabel('5 Hz Pr Train 1')
title('Comparison between conditions in Is input')
xlim([0 1])
ylim([0 1])
x = [0;0.05;1];
y = x
hold on
plot(x,y);

% ticksPr = [1;2];
% col1 = ones(length(Conc_Is_Pr),1);
% Conc_Is_Pr_Col = [col1 Conc_Is_Pr.'];
% col2 = ones(length(Conc_Is_Pr_HF_1),1);
% col2 = col2 + 1;
% Conc_Is_Pr_Col_HF_1 = [col2 Conc_Is_Pr_HF_1.'];
% 
% conckey = [Conc_Is_Pr_Col;Conc_Is_Pr_Col_HF_1];
% scatter(conckey(:,1),conckey(:,2),'.');
% xlim([0 3])
% xticks([1 2])
% xticklabels({'0.2 Hz','5 Hz Train 1'}) 

h1 = histogram(Conc_Ib_Pr);
h1.FaceColor = 'red';
h1.BinWidth = 0.05;
hold on
h2 = histogram(Conc_Is_Pr);
h2.FaceColor = 'blue';
h2.BinWidth = 0.05;

%Compare Spont Fs distribution
%Ib Fs
Ib_Fs = cell(10,1);
for nn = 1:10
Ib_Fs{nn,1} = ExperimentSet_Reduced(3).Grouped_Data_Reduced(nn).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(1).Recording.Spont_Fs 
end

Conc_Ib_Fs = Ib_Fs{1,1};
for uu = 2:10
AddPr = Ib_Fs{uu,1}
Conc_Ib_Fs = [Conc_Ib_Fs AddPr]
end


%Is Fs

Is_Fs = cell(10,1);
for nn = 1:5
Is_Fs{nn,1} = ExperimentSet_Reduced(3).Grouped_Data_Reduced(nn).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(1).Recording.Spont_Fs
end
for nn = 7:10
Is_Fs{nn,1} = ExperimentSet_Reduced(3).Grouped_Data_Reduced(nn).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(1).Recording.Spont_Fs 
end

Conc_Is_Fs = Is_Fs{1,1}
for uu = 2:10
AddPr = Is_Fs{uu,1}
Conc_Is_Fs = [Conc_Is_Fs AddPr]
end



h3 = histogram(Conc_Ib_Fs);
h3.FaceColor = 'red';
h3.BinWidth = 0.05;
hold on
h4 = histogram(Conc_Is_Fs);
h4.FaceColor = 'blue';
h4.BinWidth = 0.05;

%look at 5 hz data 


Is_HF_Data = cell(10,5);
for ll = 1:5
    for kk = 1:5
        Is_HF_Data{ll,kk} =  ExperimentSet_Reduced(3).Grouped_Data_Reduced(ll).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(4).Recording(kk).Evoked_Pr 
    end
end

for ll = 7:10
    for kk = 1:5
        Is_HF_Data{ll,kk} =  ExperimentSet_Reduced(3).Grouped_Data_Reduced(ll).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(4).Recording(kk).Evoked_Pr 
    end
end
newmat = ones(1,length(Is_HF_Data{1,1}));

%Animal 1
newmat1 = ones(1,length(Is_HF_Data{1,1}));
for ii = 1:5
    newmat1= [newmat1;Is_HF_Data{1,ii}];
end
newmat1(1,:) = []
traincol = [1;2;3;4;5];
newmat1 = [traincol newmat1];
%Animal 2
newmat2 = ones(1,length(Is_HF_Data{2,1}));
for ii = 1:5
    newmat2= [newmat2;Is_HF_Data{2,ii}];
end
newmat2(1,:) = []
% traincol = [1;2;3;4;5];
% newmat2 = [traincol newmat2];
%Animal 3
newmat3 = ones(1,length(Is_HF_Data{3,1}));
for ii = 1:5
    newmat3= [newmat3;Is_HF_Data{3,ii}];
end
newmat3(1,:) = []
% traincol = [1;2;3;4;5];
% newmat3 = [traincol newmat3];
%Animal 4
newmat4 = ones(1,length(Is_HF_Data{4,1}));
for ii = 1:5
    newmat4= [newmat4;Is_HF_Data{4,ii}];
end
newmat4(1,:) = []
% traincol = [1;2;3;4;5];
% newmat4 = [traincol newmat4];
%Animal 5
newmat5 = ones(1,length(Is_HF_Data{5,1}));
for ii = 1:5
    newmat5= [newmat5;Is_HF_Data{5,ii}];
end
newmat5(1,:) = []
% traincol = [1;2;3;4;5];
% newmat5 = [traincol newmat5];

plot(newmat5(:,1), newmat5(:,2:end)) 
%Animal 7
newmat7 = ones(1,length(Is_HF_Data{7,1}));
for ii = 1:5
    newmat7= [newmat7;Is_HF_Data{7,ii}];
end
newmat7(1,:) = []
% traincol = [1;2;3;4;5];
% newmat7 = [traincol newmat7];

plot(newmat7(:,1), newmat7(:,2:end))
%Animal 8
newmat8 = ones(1,length(Is_HF_Data{8,1}));
for ii = 1:5
    newmat8= [newmat8;Is_HF_Data{8,ii}];
end
newmat8(1,:) = []
% traincol = [1;2;3;4;5];
% newmat8 = [traincol newmat8];

plot(newmat8(:,1), newmat8(:,2:end))
%Animal 9
newmat9 = ones(1,length(Is_HF_Data{9,1}));
for ii = 1:5
    newmat9= [newmat9;Is_HF_Data{9,ii}];
end
newmat9(1,:) = []
% traincol = [1;2;3;4;5];
% newmat9 = [traincol newmat9];

plot(newmat9(:,1), newmat9(:,2:end))
%Animal 10
newmat10 = ones(1,length(Is_HF_Data{10,1}));
for ii = 1:5
    newmat10= [newmat10;Is_HF_Data{10,ii}];
end
newmat10(1,:) = []
% traincol = [1;2;3;4;5];
% newmat10 = [traincol newmat10];

plot(newmat10(:,1), newmat10(:,2:end))

newmat1 = [newmat1 newmat2 newmat3 newmat4 newmat5 newmat7 newmat8 newmat9 newmat10];
Is_consolidated_HF_Pr = newmat1;

x = [0;0.05;1];
y = x
%Animal1
subplot(1,5,1)
scatter(newmat1(1,2:end), newmat1(2,2:end),'.')
xlabel('Train 1: Pr')
ylabel('Train 2: Pr')
xlim([0 0.5])
ylim([0 0.5])
hold on
plot(x,y);
subplot(1,5,2)
scatter(newmat1(2,2:end), newmat1(3,2:end),'.')
xlabel('Train 2: Pr')
ylabel('Train 3: Pr')
xlim([0 0.5])
ylim([0 0.5])
hold on 
plot(x,y);
subplot(1,5,3)
scatter(newmat1(3,2:end), newmat1(4,2:end),'.')
xlabel('Train 3: Pr')
ylabel('Train 4: Pr')
xlim([0 0.5])
ylim([0 0.5])
hold on
plot(x,y)
subplot(1,5,4)
scatter(newmat1(4,2:end), newmat1(5,2:end),'.')
xlabel('Train 4: Pr')
ylabel('Train 5: Pr')
xlim([0 0.5])
ylim([0 0.5])
hold on
plot(x,y)
subplot(1,5,5)
scatter(newmat1(1,2:end), newmat1(5,2:end),'.')
xlabel('Train 1: Pr')
ylabel('Train 5: Pr')
xlim([0 0.5])
ylim([0 0.5])
hold on
plot(x,y)
hold on

%Ib analysis
Ib_HF_Data = cell(10,5);
for ll = 1:10
    for kk = 1:5
        Ib_HF_Data{ll,kk} =  ExperimentSet_Reduced(3).Grouped_Data_Reduced(ll).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(4).Recording(kk).Evoked_Pr 
    end
end
%Animal 1
newmat1 = ones(1,length(Ib_HF_Data{1,1}));
for ii = 1:5
    newmat1= [newmat1;Ib_HF_Data{1,ii}];
end
newmat1(1,:) = []
traincol = [1;2;3;4;5];
newmat1 = [traincol newmat1];
%Animal 2
newmat2 = ones(1,length(Ib_HF_Data{2,1}));
for ii = 1:5
    newmat2= [newmat2;Ib_HF_Data{2,ii}];
end
newmat2(1,:) = []
% traincol = [1;2;3;4;5];
% newmat2 = [traincol newmat2];
%Animal 3
newmat3 = ones(1,length(Ib_HF_Data{3,1}));
for ii = 1:5
    newmat3= [newmat3;Ib_HF_Data{3,ii}];
end
newmat3(1,:) = []
% traincol = [1;2;3;4;5];
% newmat3 = [traincol newmat3];
%Animal 4
newmat4 = ones(1,length(Ib_HF_Data{4,1}));
for ii = 1:5
    newmat4= [newmat4;Ib_HF_Data{4,ii}];
end
newmat4(1,:) = []
% traincol = [1;2;3;4;5];
% newmat4 = [traincol newmat4];
%Animal 5
newmat5 = ones(1,length(Ib_HF_Data{5,1}));
for ii = 1:5
    newmat5= [newmat5;Ib_HF_Data{5,ii}];
end
newmat5(1,:) = []
% traincol = [1;2;3;4;5];
% newmat5 = [traincol newmat5];

plot(newmat5(:,1), newmat5(:,2:end)) 
%Animal 7
newmat7 = ones(1,length(Ib_HF_Data{7,1}));
for ii = 1:5
    newmat7= [newmat7;Ib_HF_Data{7,ii}];
end
newmat7(1,:) = []
% traincol = [1;2;3;4;5];
% newmat7 = [traincol newmat7];

plot(newmat7(:,1), newmat7(:,2:end))
%Animal 8
newmat8 = ones(1,length(Ib_HF_Data{8,1}));
for ii = 1:5
    newmat8= [newmat8;Ib_HF_Data{8,ii}];
end
newmat8(1,:) = []
% traincol = [1;2;3;4;5];
% newmat8 = [traincol newmat8];

plot(newmat8(:,1), newmat8(:,2:end))
%Animal 9
newmat9 = ones(1,length(Ib_HF_Data{9,1}));
for ii = 1:5
    newmat9= [newmat9;Ib_HF_Data{9,ii}];
end
newmat9(1,:) = []
% traincol = [1;2;3;4;5];
% newmat9 = [traincol newmat9];

plot(newmat9(:,1), newmat9(:,2:end))
%Animal 10
newmat10 = ones(1,length(Ib_HF_Data{10,1}));
for ii = 1:5
    newmat10= [newmat10;Ib_HF_Data{10,ii}];
end
newmat10(1,:) = []
% traincol = [1;2;3;4;5];
% newmat10 = [traincol newmat10];

newmat1 = [newmat1 newmat2 newmat3 newmat4 newmat5 newmat7 newmat8 newmat9 newmat10];
Ib_consolidated_HF_Pr = newmat1;

x = [0;0.05;1];
y = x
%Animal1
subplot(1,5,1)
scatter(newmat1(1,2:end), newmat1(2,2:end),'.')
xlabel('Train 1: Pr')
ylabel('Train 2: Pr')
xlim([0 0.5])
ylim([0 0.5])
hold on
plot(x,y);
subplot(1,5,2)
scatter(newmat1(2,2:end), newmat1(3,2:end),'.')
xlabel('Train 2: Pr')
ylabel('Train 3: Pr')
xlim([0 0.5])
ylim([0 0.5])
hold on 
plot(x,y);
subplot(1,5,3)
scatter(newmat1(3,2:end), newmat1(4,2:end),'.')
xlabel('Train 3: Pr')
ylabel('Train 4: Pr')
xlim([0 0.5])
ylim([0 0.5])
hold on
plot(x,y)
subplot(1,5,4)
scatter(newmat1(4,2:end), newmat1(5,2:end),'.')
xlabel('Train 4: Pr')
ylabel('Train 5: Pr')
xlim([0 0.5])
ylim([0 0.5])
hold on
plot(x,y)
subplot(1,5,5)
scatter(newmat1(1,2:end), newmat1(5,2:end),'.')
xlabel('Train 1: Pr')
ylabel('Train 5: Pr')
xlim([0 0.5])
ylim([0 0.5])
hold on
plot(x,y)
hold on