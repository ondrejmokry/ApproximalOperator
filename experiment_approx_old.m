% audio_inpainting tests Dörfler Datasets
% this script reproduces the results presented in the paper. 
% Please note, that the LTFAT toolbox (http://ltfat.sourceforge.net/) has 
% to be included and added to the path.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('NSGToolbox_V010');
nsgt_startup;

global approximal;
approximal = true;

global new_settings;
new_settings = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('methods');
if (exist('dgtreal') ~= 2)
    error('Please include LTFAT toolbox');
else
    warning off; % some methods in the ltfat toolbox might produce warnings
end

snum = [1,2,3,4];
p = 0.8; % percentage of coefficients to delete

tauvec = [7.5e-5 1e-4 2.5e-4 5e-4 7.5e-4 1e-3 2.5e-3 5e-3 1e-2 2.5e-2 5e-2 7.5e-2 1e-1 2.5e-1 5e-1 7.5e-1 1 5 10];

transvec = [{'gab'},{'erb'},{'wav'}];

algovec = [{'forward-backward'},{'douglas-rachford'}];
synthesisvec = [{'synthesis'}, {'analysis'}];

% Gabor
settings.a = 160;
settings.M = 3125;
settings.width = 23.22;

% ERBLET
settings.ebins = 18;
settings.qvar = 0.08;

% Wavelet
settings.fmin = 100;
settings.wbins = 120;
settings.bw = 3;

settings.thres = 'lasso';
settings.plotting = 0;

res = zeros(length(snum),length(algovec),length(synthesisvec),length(transvec),length(tauvec));
res_relnorm = cell(length(snum),length(algovec),length(synthesisvec),length(transvec),length(tauvec));

zzz = 1;
zzztotal = numel(res(:));
fix(clock)
for iii = 1:length(snum)
    % get signal
    filename = ['./wav/sig_' num2str(snum(iii)) '.wav'];
    [s,fs] = audioread(filename);
    s = s(1:100000);
    L = length(s);
    % get mask
    rng(10);
    Mask = rand(size(s)) > p;
    depl_s = Mask.*s;
    M = logical(1-Mask);
    snr_m = @(sol) 20 *log10(std(s(M))/std(sol(M)-s(M)));
    
    for jjj = [2 1]
        settings.algo = algovec{jjj};
        for kkk = 1:2
            settings.synthesis = mod(kkk,2);
            for lll = 1:length(transvec)
                settings.trans = transvec{lll};
                tmpmax = 0;
                for mmm = 1:length(tauvec)
                    tau = tauvec(mmm);
                    [sol,relnormvec] = audio_inpainting(s,depl_s,Mask,fs,tau,settings);
                    tmp = snr_m(sol);
                    res (iii,jjj,kkk,lll,mmm) = tmp;
                    res_relnorm{iii,jjj,kkk,lll,mmm} = relnormvec;
                    fprintf('\n -- %d / % d -- \n',zzz,zzztotal);
                    fprintf('%d: %s (%5s) - SNR: %2.2f',snum(iii),settings.trans,settings.thres, snr_m(sol));
                    zzz = zzz + 1;
                    
                    if tmp > tmpmax
                        tmpmax = tmp;
                    elseif abs(tmp) < 0.66*tmpmax && abs(tmp) > eps
                        fprintf('\ndecreasing value, cont.\n');
                        zzz = zzz + (length(tauvec)-mmm);
                        break;
                    end
                end
            end
            
            % save results
            if approximal
                if new_settings
                    save('results/02 approximal/Experiments_DRvsFistavsSynvsAna.mat','res');
                else
                    save('results/01 approximal/Experiments_DRvsFistavsSynvsAna.mat','res');
                end
            else
                if new_settings
                    save('results/02 proximal/Experiments_DRvsFistavsSynvsAna.mat','res');
                else
                    save('results/01 proximal/Experiments_DRvsFistavsSynvsAna.mat','res');
                end
            end
        end
    end
end