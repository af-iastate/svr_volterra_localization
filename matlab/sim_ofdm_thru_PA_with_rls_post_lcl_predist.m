clc; clearvars; close all;

directory_name = 'simulations/rls_lcl_post_ofdm2';

M_predists = 1:17;
K_predist = 3;
fcs = [19e6, 27e6, 35e6];
nFFTs = [64, 128, 256];
upsample_amnts = [1, 2, 4, 8, 10, 15, 20, 30, 40, 60, 120];

for fc=fcs
    for nFFT=[nFFTs(2)]
        for M_predist=M_predists
            mat_filename = sprintf('rls_lcl_post_ofdm_fc_%g_nFFT_%g_M_%g_K_%g.mat', fc, nFFT, M_predist, K_predist);
            mat_filename = strrep(mat_filename, '+', '');
            
            fullpath = [directory_name, '/', mat_filename];
            
            fprintf('--- NEW SIM:  fc=%g  nFFT=%g  M=%g  K=%g ---\n\n', fc, nFFT, M_predist, K_predist)
            if ~exist(fullpath, 'file')
                close all;
                simRound_ofdm_thru_PA_with_rls_post_lcl_predist(fc, nFFT, M_predist, K_predist, fullpath)
                pause(3)
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