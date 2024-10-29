function xCompanionNew = applyTableKernelPM(xCompanion, xSlice, M, mtable)
    xCompanionNew = zeros(size(xCompanion));
    xCompanionNew(1:M) = xSlice;
    jj = M;
    for ii=1:size(mtable,1)
        jj = jj + 1;
        xCompanionNew(jj) = xCompanionNew(mtable(ii, 1)) * xSlice(mtable(ii, 2));
    end
end