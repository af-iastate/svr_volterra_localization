clc; clearvars; close all;

directory_name = 'simulations/gbl_multitone';

M_predists = 1:17;
K_predist = 3;
nFFTs = [64, 128, 256];

for nFFT=nFFTs
    for M_predist=M_predists
            mat_filename = sprintf('gbl_multitone_M_%g_nFFT_%g.mat', M_predist, nFFT);
            mat_filename = strrep(mat_filename, '+', '');
            
            fullpath = [directory_name, '/', mat_filename];
            
            fprintf('--- NEW SIM:  Mpred=%g  nFFT=%g ---\n\n', M_predist, nFFT)
            if ~exist(fullpath, 'file')
                close all;
                simRound_multitone_thru_PA_train(M_predist, K_predist, nFFT, fullpath)
                pause(3)
            end
    end
end
