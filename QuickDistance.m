function L_S_Array2=QuickDistance(features)
[MyRows,NumberOfColumns]=size(features);
NumberOfSubjects=MyRows/2;
D=pdist2(features(1:NumberOfSubjects,:),features(NumberOfSubjects+1:MyRows,:),'cosine');
MinD=min(min(D));
MaxD=max(max(D));
ScaledD=(D-MinD)/(MaxD-MinD);
S=1-ScaledD;
L_S_Array2=zeros(power(NumberOfSubjects,2),2);
count=1;
for i = 1:NumberOfSubjects
    for j = 1:NumberOfSubjects
        L_S_Array2(count,1) = S(i,j);
        if i == j;L_S_Array2(count,2)=1.0;else,L_S_Array2(count,2)=0.0;end
        count=count+1;
    end
end
return