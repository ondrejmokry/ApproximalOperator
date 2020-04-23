addpath('NSGToolbox_V010');
nsgt_startup;

%% load the signal
snum = 2;
% 1: Strings
% 2: Piano
% 3: Percussion
% 4: Jazz

filename = ['./wav/sig_' num2str(snum) '.wav'];
[ orig_sig, fs ] = audioread(filename);
if snum==4
    orig_sig = orig_sig(1:200000);
else
    orig_sig = orig_sig(1:100000);
end

%% mask of the reliable samples
s = fs;
f = fs + 800;
mask = ones(size(orig_sig));
mask(s:f) = 0;
mask = logical(mask);
% p = 0.8; % percentage of coefficients to delete
% rng(10);
% mask = rand(size(orig_sig)) > p;

%% degrade the signal
depl_sig = orig_sig.*mask;

%% settings
tau = 0.1;

transvec = [{'gab'},{'erb'},{'wav'}];
transnum = 1;
threshvec = [{'lasso'},{'ew'},{'wgl'},{'pew'}];
threshnum = 1;
algovec = [{'forward-backward'},{'douglas-rachford'}];
algonum = 2;
synthesisvec = [{'synthesis'}, {'analysis'}];
synthesisnum = 2;

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

switch snum
    case 1
        settings.perslength = 72;
        settings.percussionpersistence = 0;
    case 2
        settings.perslength = 48;
        settings.percussionpersistence = 0;
    case 3
        settings.perslength = 24;
        settings.percussionpersistence = 1;
    case 4
        settings.perslength = 48;
        settings.percussionpersistence = 0;
end

% trans & thresh
settings.trans = transvec{transnum};
settings.thres = threshvec{threshnum};

% plot
settings.plotting = 0;

% analysis / synthesis
settings.synthesis = strcmp('synthesis',synthesisvec{synthesisnum});

% algorithm
settings.algo = algovec{algonum};

%% SNR
M = ~mask;
snr_m = @(sol) 20 *log10(std(orig_sig(M))/std(sol(M)-orig_sig(M)));

%% reconstruction
addpath('Methods');
global approximal;
approximal = true; %#ok<NASGU>
[ sol1, rel_norm_vec1, obj1 ] = audio_inpainting( orig_sig, depl_sig, mask, fs, tau, settings );
fprintf('\n%d: %s (%5s) approximal solution, SNR: %2.2f\n\n',snum,settings.trans,settings.thres, snr_m(sol1));
sol1 = real(sol1);

approximal = false;
[ sol2, rel_norm_vec2, obj2 ] = audio_inpainting( orig_sig, depl_sig, mask, fs, tau, settings );
fprintf('\n%d: %s (%5s) proximal solution,   SNR: %2.2f\n\n',snum,settings.trans,settings.thres, snr_m(sol2));
sol2 = real(sol2);

%% plot
figure;
% plot(orig_sig(s-50:f+50)); hold on;
% plot(sol(s-50:f+50));
plot(orig_sig); hold on;
plot(sol1); hold on;
plot(sol2);
legend('original signal','approximal solution','proximal solution')