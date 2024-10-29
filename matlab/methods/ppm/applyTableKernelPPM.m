function xCompanionNew = applyTableKernelPPM(xCompanion, xNew, K, mreassign, mpermute, mproduct)
    arguments
        xCompanion (:,:) 
        xNew (1,1)
        K (1,1) double {mustBeInteger, mustBePositive}
        mreassign
        mpermute
        mproduct
    end

    [N, ~] = size(xCompanion);
    xCompanionNew = zeros(N, 1);
    if isa(xNew,'sym')
        xCompanionNew = sym(xCompanionNew);
    end
    xPowers = xNew.^(1:K);

    % Do Re-assign
    for ii=1:size(mreassign, 1)
        xCompanionNew(mreassign(ii,1)) = xPowers(mreassign(ii,2));
    end

    % Do Permute
    for ii=1:size(mpermute, 1)
        xCompanionNew(mpermute(ii,1)) = xCompanion(mpermute(ii,2));
    end

    % Do Products
    for ii=1:size(mproduct, 1)
        xCompanionNew(mproduct(ii,1)) = xCompanion(mproduct(ii,2)) * xPowers(mproduct(ii,3));
    end
    
end