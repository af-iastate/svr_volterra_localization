clc; clearvars; close all;

directory_name = 'simulations/svr_swap';

M_scl = 2;
fcs = [19e6, 27e6, 35e6];
nFFTs = [64, 128, 256];
upsample_amnts = [1, 2, 4, 8, 10, 15, 20, 30, 40, 60, 120];
M_svrs = 1:8;
nFFT = 128;

for fc_test=fcs
    for fc_train=fcs
        for upsample_amnt=upsample_amnts
            for M_svr=M_svrs
                
    
                mat_filename = sprintf('svr_res_fc_test_%g_fc_train_%g_nFFT_%g_up_%d_Msvr_%g.mat', fc_test, fc_train, nFFT, upsample_amnt, M_svr);
                mat_filename = strrep(mat_filename, '+', '');
                
                fullpath = [directory_name, '/', mat_filename];
                
                fprintf('--- NEW SIM:  fc_test=%g  fc_train=%g  upsamp=%g  nFFT=%g  Msvr=%g ---\n\n', fc_test, fc_train, upsample_amnt, nFFT, M_svr)
                if ~exist(fullpath, 'file')
                    simRound_swap_ofdm_thru_PA_with_svr_predist(fc_test, fc_train, upsample_amnt, nFFT, M_svr, fullpath)
                end
                drawnow
                pause(.05)
            end
        end
    end
end

