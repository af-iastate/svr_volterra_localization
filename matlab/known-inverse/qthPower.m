function [hOut, hOutStruct] = qthPower(hInput, MInput, KInput, Q)


    NInput = length(hInput);
    hInputMultiIndices = makeIndices(MInput, 1:KInput);

    hCurrentMultiIndices = hInputMultiIndices;
    hCurrent = hInput;
    NCurrent = NInput;
    
    hNext = makeKernelStruct(MInput, KInput, hCurrentMultiIndices);
    [NCurrentMultiIndices, ~] = size(hCurrentMultiIndices);
    for ii=1:NCurrentMultiIndices
        multiIndex = hCurrentMultiIndices(ii,:);
        hNext.set(multiIndex, hCurrent(ii));
    end

    for q=2:Q

        KCurrent = q*KInput;
        [NCurrentMultiIndices, ~] = size(hCurrentMultiIndices);
        hNextMultiIndices = makeIndices(MInput, 1:KCurrent);
        hNext = makeKernelStruct(MInput, KCurrent, hNextMultiIndices);

        % Do Convolution
        NConv = NInput + NCurrent - 1;
        for n=0:NConv-1
            for m=0:n
                mCurrent = m;
                mInput = n - m;   % flipped and shifted
                if mInput >= 0  && mInput < NInput && mCurrent < NCurrentMultiIndices
                    multiIndexNext = hCurrentMultiIndices(mCurrent + 1, :) ...
                        + hInputMultiIndices(mInput + 1, :);
                    hValueNext = hCurrent(mCurrent + 1) * hInput(mInput + 1);
                    hNext.set( ...
                        multiIndexNext, ...
                        hNext.get(multiIndexNext) + hValueNext ...
                    );
                end
            end
        end
        hCurrent = hNext.flatten();
        hCurrentMultiIndices = hNext.multiIndices;
        NCurrent = length(hCurrent);
    end
    
    hOut = hCurrent;
    hOutStruct = hNext;
end


function multiIndices = makeIndices(M, p)
    multiIndices = cell([length(p), 1]);
    for ii=1:length(p)
        tauArr = generateTauUpperTri(M, p(ii));
        [N, ~] = size(tauArr);
        multiIndices{ii} = arrayfun(@(ii) ...
            arrayfun(@(x) sum(tauArr(ii,:)==x), 0:M-1), ...
            1:N, ...
            'UniformOutput', false ...
        );
        multiIndices{ii} = cell2mat(multiIndices{ii}.');
    end
    multiIndices = vertcat(multiIndices{:});
    
end


function h = makeKernelStruct(M, K, multiIndices)
    h = struct;
    h.M = M;
    h.K = K;
    h.f_hash = @num2str;
    h.multiIndices = multiIndices;
    [N, ~] = size(multiIndices);
    keys = arrayfun(@(ii) ...
        h.f_hash(multiIndices(ii, :)), ...
        1:N, ...
        'UniformOutput', false ...
    );
    values = num2cell(zeros([N, 1]));
    h.map = containers.Map(keys, values);
    function y = f_get(x)
        y = h.map(h.f_hash(x));
    end
    function f_set(x, y)
        h.map(h.f_hash(x)) = y;
    end
    function hFlat = f_flatten()
        [NMulti, ~] = size(multiIndices);
        hFlat = arrayfun(@(ii) ...
            f_get(multiIndices(ii, :)), ...
            1:NMulti ...
        );
    end
    h.get = @f_get;
    h.set = @f_set;
    h.flatten = @f_flatten;
end

