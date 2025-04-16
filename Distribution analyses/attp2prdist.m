%% Combine all NMJ data into one (Attp2)
    AllLFPrIb = [];
for k=1:numel(fn)
    TempLFPrIBD = LFPrIb.(fn{k});
    AllLFPrIb = vertcat(AllLFPrIb,TempLFPrIBD);
end

AllHF1PrIb = [];
for k=1:numel(fn)
    TempHFPrIBD = HFPr1Ib.(fn{k});
    AllHF1PrIb = vertcat(AllHF1PrIb,TempHFPrIBD);
end
AllHF1PrIb(AllHF1PrIb>=1)=NaN;

AllHF2PrIb = [];
for k=1:numel(fn)
    TempHFPrIBD = HFPr2Ib.(fn{k});
    AllHF2PrIb = vertcat(AllHF2PrIb,TempHFPrIBD);
end
AllHF2PrIb(AllHF2PrIb>=1)=NaN;

AllHF3PrIb = [];
for k=1:numel(fn)
    TempHFPrIBD = HFPr3Ib.(fn{k});
    AllHF3PrIb = vertcat(AllHF3PrIb,TempHFPrIBD);
end
AllHF3PrIb(AllHF3PrIb>=1)=NaN;

AllHF4PrIb = [];
for k=1:numel(fn)
    TempHFPrIBD = HFPr4Ib.(fn{k});
    AllHF4PrIb = vertcat(AllHF4PrIb,TempHFPrIBD);
end
AllHF4PrIb(AllHF4PrIb>=1)=NaN;

AllHF5PrIb = [];

for k=1:numel(fn)
    TempHFPrIBD = HFPr5Ib.(fn{k});
    AllHF5PrIb = vertcat(AllHF5PrIb,TempHFPrIBD);
end
AllHF5PrIb(AllHF5PrIb>=1)=NaN;

    ax1  = subplot(1,6,1)
    hist(AllLFPrIb);
    xlim([0 1])
    ylim([0 1000])
    ax2 = subplot(1,6,2)
    hist(AllHF1PrIb);
    xlim([0 1])
    ylim([0 1000])
    ax3 = subplot(1,6,3)
    hist(AllHF2PrIb);
    xlim([0 1])
    ylim([0 1000])
    ax4 = subplot(1,6,4)
    hist(AllHF3PrIb);
    xlim([0 1])
    ylim([0 1000])
    ax5 = subplot(1,6,5)
    hist(AllHF4PrIb);
    ylim([0 1000])
    xlim([0 1])
    ax6 = subplot(1,6,6)
    hist(AllHF5PrIb);
    xlim([0 1])
    ylim([0 1000])
    

%     pd1 = fitdist(AllLFPrIbD,'Exponential');
%     pd2 = fitdist(AllHF1PrIbD,'Exponential');
%     pd3 = fitdist(AllHF2PrIbD,'Exponential');
%     pd4 = fitdist(AllHF3PrIbD,'Exponential');
%     pd5 = fitdist(AllHF4PrIbD,'Exponential');
%     pd6 = fitdist(AllHF5PrIbD,'Exponential');

    
    c1 = uisetcolor
    c2 = uisetcolor
    c3 = uisetcolor
    
    ax1  = subplot(2,6,1)
    histogram(AllLFPrIb,'BinWidth',0.05);
    xlim([0 1])
    ylim([0 800])
    xlabel('Pr')
    ylabel('Count')
    ax2 = subplot(2,6,2)
    histogram(AllHF1PrIb,'BinWidth',0.05);
    xlim([0 1])
    ylim([0 800])
    ax3 = subplot(2,6,3)
    histogram(AllHF2PrIb,'BinWidth',0.05);
    xlim([0 1])
    ylim([0 800])
    ax4 = subplot(2,6,4)
    histogram(AllHF3PrIb,'BinWidth',0.05);
    xlim([0 1])
    ylim([0 800])
    ax5 = subplot(2,6,5)
    histogram(AllHF4PrIb,'BinWidth',0.05);
    ylim([0 800])
    xlim([0 1])
    ax6 = subplot(2,6,6)
    histogram(AllHF5PrIb,'BinWidth',0.05);
    xlim([0 1])
    ylim([0 800])
    ax7 = subplot(2,6,[7 8 9]);
    [f1,x1]=ecdf(AllLFPrIb);
    plot(x1,f1,'Color',c1,'LineWidth',2)
    hold on 
    [f2,x2]=ecdf(AllHF1PrIb);
    plot(x2,f2,'Color',c2,'LineWidth',2)
    xlabel('Pr')
    ylabel('Probability')
    legend('0.2 Hz Evoked Pr','5 Hz 1 Evoked Pr')

    ax8 = subplot(2,6,[10 11 12]);
    [f2,x2]=ecdf(AllHF1PrIb);
    plot(x2,f2,'Color',c2,'LineWidth',2)
    hold on 
    [f3,x3]=ecdf(AllHF5PrIb);
    plot(x3,f3,'Color',c3,'LineWidth',2)
    xlabel('Pr')
    ylabel('Probability')
    legend('5 Hz 1 Evoked Pr','5 Hz 5 Evoked Pr')
    [h,p] = kstest2(AllHF1PrIb,AllLFPrIb);
    [h2,p2] = kstest2(AllHF1PrIb,AllHF5PrIb)

% subplot(6,5,1)
% plot(rand(100,1)); hold on
% plot(rand(100,1),'r')
% subplot(6,5,2)
% plot(rand(100,1)); hold on
% plot(rand(100,1),'r')
% subplot(6,5,3)
% plot(rand(100,1)); hold on
% plot(rand(100,1),'r')
% subplot(6,5,4)
% plot(rand(100,1)); hold on
% plot(rand(100,1),'r')
% subplot(5,6,4)
% plot(rand(100,1)); hold on
% plot(rand(100,1),'r')
