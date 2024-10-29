clc; clearvars; %close all;

directory_name = 'simulations/svr_resolution';

M_scl = 2;
fcs = [19e6, 27e6, 35e6];
nFFTs = [64, 128, 256];
upsample_amnts = [1, 2, 4, 8, 10, 15, 20, 30, 40, 60, 120];

% for fc=fcs
%     for upsample_amnt=upsample_amnts
%         for nFFT=nFFTs
% 
%             upfactor = 120;
%             svr_factor = upfactor / upsample_amnt;
%             M_svr = M_scl * (upfactor / svr_factor);
% 
%             mat_filename = sprintf('svr_res_fc_%g_nFFT_%g_up_%d_Msvr_%g.mat', fc, nFFT, upsample_amnt, M_svr);
%             mat_filename = strrep(mat_filename, '+', '');
% 
%             fullpath = [directory_name, '/', mat_filename];
% 
%             fprintf('--- NEW SIM:  fc=%g  upsamp=%g  nFFT=%g  Msvr=%g ---\n\n', fc, upsample_amnt, nFFT, M_svr)
%             if ~exist(fullpath, 'file')
%                 simRound_ofdm_thru_PA_with_svr_predist(fc, upsample_amnt, nFFT, M_svr, fullpath)
%             end
%         end
%     end
% end

fc = 27e6;
nFFT = 128;
M_svrs = 1:8;
for upsample_amnt=upsample_amnts
    for jj=1:length(fcs)
        fc = fcs(jj);
    
        for ii=1:length(M_svrs)
            M_svr = M_svrs(ii);
            mat_filename = sprintf('svr_res_fc_%g_nFFT_%g_up_%d_Msvr_%g.mat', fc, nFFT, upsample_amnt, M_svr);
            mat_filename = strrep(mat_filename, '+', '');
            
            fullpath = [directory_name, '/', mat_filename];
            
            fprintf('--- NEW SIM:  fc=%g  upsamp=%g  nFFT=%g  Msvr=%g ---\n\n', fc, upsample_amnt, nFFT, M_svr)
            if ~exist(fullpath, 'file')
                simRound_ofdm_thru_PA_with_svr_predist(fc, upsample_amnt, nFFT, M_svr, fullpath)
            end
            drawnow;
            pause(5)
        end
    end
end