plot(EventDetectionStruct.EventDetectionSettings_Output.TemplateGeneration.All_Episodes_DeltaFF0_Mean)
plot(EventDetectionStruct.EventDetectionSettings_Output.TemplateGeneration.All_Episodes_DeltaFF0_Mean_Norm)  

%Load in EventDetectionStruct from _EventDetectionData_Ep_[Number].mat
%Calculate area under dF/F whole NMJ mean
mean_trace = EventDetectionStruct.EventDetectionSettings_Output.TemplateGeneration.All_Episodes_DeltaFF0_Mean;
mean_trace_norm = EventDetectionStruct.EventDetectionSettings_Output.TemplateGeneration.All_Episodes_DeltaFF0_Mean_Norm;

max_x = length(mean_trace);
x_vec = 1:1:max_x;
Int = trapz(x_vec,mean_trace);
Int_norm = trapz(x_vec,mean_trace_norm);

area(x_vec,mean_trace_norm)
xlabel('Frame Number')
ylabel('mean NMJ dF/F')