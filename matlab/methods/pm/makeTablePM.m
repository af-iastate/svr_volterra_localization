function prodTab = makeTablePM(M, K)
% From "Algorithm 2" given by 
%     H. Enzinger, K. Freiberger, G. Kubin and C. Vogel, 
%     "Fast time-domain volterra filtering," 
%     2016 50th Asilomar Conference on Signals, Systems and Computers, 
%     Pacific Grove, CA, USA, 2016, pp. 225-228, 
%     doi: 10.1109/ACSSC.2016.7869029. 
    MM = M - 1;
    nSize = sum(arrayfun(@(k) nchoosek(MM+k, k), 2:K));
    prodTab = zeros(nSize, 2);
    ii = 1;
    jj = 1;
    for k = 2:K 
        m = zeros(1, k+1);
        goto_start = true;
        while (goto_start)  % start:
            prodTab(ii,:) = [jj, m(k+1) + 1];
            ii = ii + 1;
            if (m(k+1) < MM)
                m(k+1) = m(k+1) + 1;
                continue;
            end
            jj = jj + 1;
            goto_start = false;
            for kk=k-1:-1:1
                if (m(kk+1) < MM)
                    m(kk+1) = m(kk+1) + 1;
                    for pkk=kk+1:k
                        m(pkk+1) = m(pkk);
                    end
                    goto_start = true;
                    break;
                end
            end
            if goto_start
                continue;
            end
            
            break;
        end
    end
end