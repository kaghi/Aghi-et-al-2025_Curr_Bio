
All_Locations = cell(1,10);
for ii = 1:10
S= ExperimentSet_Reduced(3).All_Pair_Structures(ii).QuaSOR_STORM_Pair_Structure_Sorted_Verified;
SizeLoc = length(S);
temp_mat = zeros(SizeLoc,3);
    for kk = 1:SizeLoc
    Temp= getfield(S(kk),'PixelMatched_QuaSOR_Coord')
    temp_mat(kk,1) = Temp(1);
    temp_mat(kk,2) = Temp(2);
    temp_mat(kk,3) = getfield(S(kk),'BoutonType')
    end
All_Locations{ii} = temp_mat;
end  

T = cell2table(All_Locations);
NewTab = table
for nn = 1:10
    NewTable =
end




