clc; clearvars; close all;

directory_name = 'simulations/rls_lcl_ofdm_swap';

M_predists = 1:8;
K_predist = 3;
fcs = [19e6, 27e6, 35e6];
nFFTs = [64, 128, 256];
upsample_amnts = [1, 2, 4, 8, 10, 15, 20, 30, 40, 60, 120];

for fc_test=fcs
    for fc_train=fcs
        for upsample_amnt=upsample_amnts
            for nFFT=[nFFTs(2)]
                for M_predist=M_predists
                    mat_filename = sprintf('rls_lcl_ofdm_swap_fc_test_%g_fc_train_%g_nFFT_%g_up_%d_M_%g_K_%g.mat', fc_test, fc_train, nFFT, upsample_amnt, M_predist, K_predist);
                    mat_filename = strrep(mat_filename, '+', '');
                    
                    fullpath = [directory_name, '/', mat_filename];
                    
                    fprintf('--- NEW SIM:  fc_test=%g  fc_train=%g  upsamp=%g  nFFT=%g  M=%g  K=%g ---\n\n', fc_test, fc_train, upsample_amnt, nFFT, M_predist, K_predist)
                    if ~exist(fullpath, 'file')
                        close all;
                        simRound_swap_ofdm_thru_PA_with_rls_lcl_predist(fc_test, fc_train, upsample_amnt, nFFT, M_predist, K_predist, fullpath)
                        drawnow
                        pause(0.25)
                    end
                end
            end
        end
    end
end

% fc = 19e6;
% nFFT = 128;
% upsample_amnt = 2;
% M_svrs = 2* upsample_amnts;
% for M_svr=M_svrs
%     mat_filename = sprintf('svr_res_fc_%g_nFFT_%g_up_%d_Msvr_%g.mat', fc, nFFT, upsample_amnt, M_svr);
%     mat_filename = strrep(mat_filename, '+', '');
% 
%     fullpath = [directory_name, '/', mat_filename];
% 
%     fprintf('--- NEW SIM:  fc=%g  upsamp=%g  nFFT=%g  Msvr=%g ---\n\n', fc, upsample_amnt, nFFT, M_svr)
%     if ~exist(fullpath, 'file')
%         simRound_ofdm_thru_PA_with_svr_predist(fc, upsample_amnt, nFFT, M_svr, fullpath)
%     end
% end