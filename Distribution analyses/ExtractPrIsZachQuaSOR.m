
All_Is_Data = [];
for ll = 1:10
NMJNo = ll;
if isempty(ExperimentSet_Reduced(3).Grouped_Data_Reduced(NMJNo).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(1).Evoked_Pr);
All_Is_Data{ll} = [];
else
    Low_Pr_Is = ExperimentSet_Reduced(3).Grouped_Data_Reduced(NMJNo).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(1).Recording.Evoked_Pr;
    High_Pr_Is = [];
        for nn = 1:5
        High_Pr_Is{nn} = ExperimentSet_Reduced(3).Grouped_Data_Reduced(NMJNo).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(4).Recording(nn).Evoked_Pr;
        end
    Coords = [];
    BoutonType = [];
        for ii = 1:length(ExperimentSet_Reduced(3).All_Pair_Structures(NMJNo).QuaSOR_STORM_Pair_Structure_Sorted_Verified)
        Coords{ii} = ExperimentSet_Reduced(3).All_Pair_Structures(NMJNo).QuaSOR_STORM_Pair_Structure_Sorted_Verified(ii).PixelMatched_STORM_Coord;
        BoutonType{ii} = ExperimentSet_Reduced(3).All_Pair_Structures(NMJNo).QuaSOR_STORM_Pair_Structure_Sorted_Verified(ii).BoutonType;
        end
    BArray = cell2mat(BoutonType);
    CoordX = [];
    CoordY = [];
        for kk = 1:length(Coords)
        CoordX(kk) = Coords{1,kk}(1);
        CoordY(kk) = Coords{1,kk}(2);
        end
    IsInd = find(BArray(1,:) == 2);
    CArray = [];
    CArray = [CoordX(IsInd); CoordY(IsInd)];
    IsArray = [CArray; Low_Pr_Is];
        for nnn = 1:5
        IsArray = [IsArray; High_Pr_Is{nnn}];
        end
    IsArray = IsArray.';
    All_Is_Data{ll} = IsArray;
end
end
