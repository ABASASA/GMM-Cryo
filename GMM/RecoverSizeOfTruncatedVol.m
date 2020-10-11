function [Arecover] = RecoverSizeOfTruncatedVol(ATruncated,AFull)

A_padded     = AFull;
for j=1:length(AFull{1})
    if j<=length(ATruncated{1})
        A_padded{1}{j}     = ATruncated{1}{j};
    else
        A_padded{1}{j}     = zeros(size(AFull{1}{j}));
    end
end
Arecover = A_padded;

end