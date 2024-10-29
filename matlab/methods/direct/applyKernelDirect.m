function xCompanion = applyKernelDirect(xSlice, xCompanion, mtable)
    nSize = length(mtable);
    for ii=1:nSize
        xCompanion(ii) = prod(xSlice.' .^ mtable(ii,:));
    end
end