function [ sol, rel_norm_vec, obj ] = audio_inpainting( orig_sig, depl_sig, mask, fs, tau, settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global approximal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trans = settings.trans; % 'erb';
thres = settings.thres; % 'lasso';


if ~isfield(settings, 'algo')
    algo = 'douglas-rachford';
else
    algo = settings.algo;
end

if ~isfield(settings,'synthesis')
    synthesis = 1;
else
    synthesis = settings.synthesis;
end

plotting = settings.plotting;
fignum = 10;
verbose = 0;

[r,~] = size(depl_sig);
if (r==1)
    depl_sig = depl_sig';
end
L = length(depl_sig);

if size(tau,2)~=1
    warmstart = 1;
    numintervals = 3;
    wfactor = logspace(tau(1),tau(2),numintervals);
    tau = 1;
else
    warmstart = 0;
end

switch trans
    case 'sqrt'
        %square root warping
        warpfun_sqrt = @(x) sign(x).*((1+abs(x)).^(1/2)-1);
        invfun_sqrt = @(x) sign(x).*((1+abs(x)).^2-1);
        bwmul = 1.5/3;
        fmax = fs/2;
        bins = 2;
        fmin = 0;
        [g,a,fc]=warpedfilters(warpfun_sqrt,invfun_sqrt,fs,fmin,fmax,bins,L,'bwmul',bwmul,'complex','fractional');
        if (L ~= filterbanklength(L,a))
            error('filterbanklength not equal L');
        end
        gd = filterbankdual(g,a,L);
        Psi = @(x) c2mat(filterbank(x,g,a),a(:,2));
        Psit= @(x) real(ifilterbank(mat2c(x,a(:,2)),gd,a));
        fprintf('bins: %d, fmin: %d, bwmul: %2.3f, red: %3.2f \n',bins,fmin,bwmul,sum(a(:,2))/L);
    case 'gab'
        a = settings.a;
        M = settings.M;
        width = settings.width;
        g = {'hann', round(width*1e-3*fs)};
        gd= gabdual(g,a,M,L);
        Psi = @(x) dgtreal(x,g,a,M,L);
        Psit= @(x) idgtreal(x,gd,a,M,L);
    case 'erb'
        bins = settings.ebins;
        Qvar = settings.qvar;
        [g,shift,M] = nsgerbwin(bins,fs,L,'Qvar',Qvar);
        gd = nsdual(g,shift,M);
        Psi = @(x) myERBlettransform (x,g,shift,M);
        Psit= @(x) myERBinvtransform (x,gd,shift,M,L);
        if plotting
            fprintf('bins: %d, Qvar: %f, red: %3.2f\n',bins,Qvar,sum(M)/L);
        end
    case 'mywav'
        wname = 'warpedGaussian';
        fmin = settings.fmin2;
        bw = settings.bw2;
        overlap = settings.overlap;
        a = 2;
        epsilon = 1e-12;
        [g,am,gd,~,fbl,~] = waveletframe5(wname,fmin,bw,overlap,a,L,fs,epsilon);
        disp(fbl);
        gd_transpose = gd';
        g_transpose = g';
        Psi  = @(x) c2mat(wvlttf(x,g_transpose,am),L./am);
        Psit = @(x) wvltitf(mat2c(x,L./am), gd_transpose, am);
        red = sum(am)/L;
        fprintf('fmin: %d, bins: %d, bw:%f, red: %3.2f\n',fmin,overlap,bw,red);
    case 'wav'
        fmin = settings.fmin;
        bins = settings.wbins;
        bw = settings.bw;
        win = 'hann';
        
        if strcmp(win,'wp2inp')
            [g,shift,M,fb] = nswvltwp2win(fmin,bw,bins,fs,L);
            gd = nsdual(g,shift,M);
            Psi = @(x) c2mat(nsgtf(x,g,shift,M),M);
            Psit= @(x) nsigtf(mat2c(x,M),gd,shift,L);
            fprintf('Framebounds: %f %f',fb(1),fb(2));
        else
            [~,g,shift,M,Ls,fb,tgtfl] = mywvltanalysis(orig_sig,fmin,fs,bins,bw,win);
            red = sum(cellfun(@numel,g))/L;
            Psi =  @(x) mywvltanalysis (x,fmin,fs,bins,bw,win);
            Psit = @(x) mywvltsynthesis(x,g,shift,M,Ls,fb,tgtfl);
            if(fb(1) < eps)
                error('lower fb = 0');
            end
            if plotting
                fprintf('fmin: %d, bins: %d, bw:%f, lfb:%f, red: %3.2f\n',fmin,bins,bw,fb(1),red);
            end
        end
    otherwise
        error('Only Gabor transform supported');
end

switch thres
    case {'lasso','ew'}
        % setting the function f1 (l1 norm of tf coefficients)
        paraml1.verbose = 0;
        paraml1.method = thres;
        if synthesis
            % synthethis model:
            f1.prox = @(x,T) myproxl1(x,tau*T,paraml1);
        else
            % analysis model:
            paraml1.A = Psi;
            paraml1.At = Psit;
            f1.prox = @(x,T) myproxl1_analysis(x,tau*T,paraml1);
        end
        objnorm = @(x) 0.5*norm((mask.*Psit(x)) - depl_sig)^2 + tau*norm(x,1);
    case {'grouplasso','gl','wgl','pew'}
        if strcmp(trans,'erb') || strcmp(trans,'wav') || strcmp(trans,'mywav') 
            % with erblets and different time shifts, only persistence in
            % time direction is possible; no frequency persistence due to
            % non-uniform sampling 
            paraml21.t_length = settings.perslength;
            paraml21.verbose = 0;
            paraml21.M = M;
            paraml21.L = L;
            paraml21.fs = fs;
            paraml21.method = thres;
            if strcmp(thres,'gl')
                f1.prox = @(x,T) myerbproxl21_gl(x,tau*T,paraml21);
            else
                f1.prox = @(x,T) myerbproxl21(x,tau*T,paraml21);
            end
            objnorm = @(x) 0.5*norm((mask.*Psit(x)) - depl_sig)^2 + tau*myl21norm(x,paraml21);
        else
            tmp = max(1,round(settings.perslength*1e-3 / (a/fs) ));
            if strcmp(thres,'gl')
                % tmp must be divisible by L
                if (tmp ~= 1)
                    tmp2 = alldiv(L/a);
                    [~,tmp2idx] = min(abs(tmp2 - tmp));
                    tmp = tmp2(tmp2idx(1));
                end
                paraml21.overlap = 0;
            else
                paraml21.overlap = 1;
            end
            if isfield(settings,'percussionpersistence')
                if settings.percussionpersistence
                    mygroupvec = ones(tmp,tmp);
                else
                    mygroupvec = ones(1,tmp);
                end
            else
                mygroupvec = ones(1,tmp);
            end
            paraml21.group = mygroupvec;
            paraml21.verbose = 0;
            paraml21.method = thres;
            f1.prox = @(x,T) myproxl21(x,tau*T,paraml21);
            objnorm = @(x) 0.5*norm((mask.*Psit(x)) - depl_sig)^2 + tau*myl21norm(x,paraml21);
        end
    otherwise
        error('only lasso method available');
end

if plotting
    dynrange = 60;
    figure(fignum);
    
    switch trans
        case 'sqrt'
            subplot(131); plotfilterbank(mat2c(Psi(orig_sig),a(:,2)),a,fc,fs,dynrange,'nocolorbar'); ylim([0,size(a,1)/2]);
            subplot(132); plotfilterbank(mat2c(Psi(depl_sig),a(:,2)),a,fc,fs,dynrange,'nocolorbar'); ylim([0,size(a,1)/2]);
        case 'mywav'
            subplot(131); plotcoeff(mat2c(Psi(orig_sig),L./am),L,fs,dynrange);
            subplot(132); plotcoeff(mat2c(Psi(depl_sig),L./am),L,fs,dynrange);
        case 'gab'
            subplot(131); plotdgtreal(Psi(orig_sig),a,M,fs,dynrange,'nocolorbar');
            subplot(132); plotdgtreal(Psi(depl_sig),a,M,fs,dynrange,'nocolorbar');
        case {'erb', 'wav'}
            subplot(131); plotnsgtf(mat2c(Psi(orig_sig),M),shift,fs,2,dynrange);
            subplot(132); plotnsgtf(mat2c(Psi(depl_sig),M),shift,fs,2,dynrange);
    end
    subplot(131); title('Original signal');
    subplot(132); title('Depleted signal');
end

if synthesis
    % synthesis model:
    % setting the function f2 (l2 norm)
    f2.grad = @(x)   Psi( mask.*(mask.*Psit(x) - depl_sig));
    f2.prox = @(x,T) Psi( Psit(x) .* (1-  mask )+ mask.* depl_sig );
    f2.eval = @(x) eps;
else
    % analysis model:
    f2.prox = @(x,T) x-mask.*x + depl_sig;
    f2.grad = @(x) mask.*(mask.*x - depl_sig);
end

% setting different parameter for the simulation
param.verbose=2; % display parameter
param.maxit=200; % maximum iteration
param.tol=1e-3; % tolerance to stop iterating

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings for the second experiment
% param.maxit = 500;
% param.tol = 1e-5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time step
param.gamma = 1;
param.lambda = 0.99;

% initialization
iter = 1;
rel_norm_vec = zeros(1,param.maxit);
rel_norm = 1;

if synthesis
    % synthethis init:
    sol = Psi(depl_sig);
else
    % analysis init:
    sol = depl_sig;
end
prev_sol = sol;

u_n = sol;
x_n = [];
obj = [];
tn = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_old = sol;
y = 0;
x_bar = sol;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while 1

    switch algo
        case 'declipping' 
            g1 = - Psi( depl_sig - (mask).*Psit(sol) );
            x_n = f1.prox(sol + (g1),param.gamma);
            sol = x_n + 0.9*(x_n - u_n);
            u_n = x_n; % updates
        case 'douglas-rachford'
            if warmstart
                param.gamma = wfactor(floor(iter/(param.maxit/numintervals))+1);
            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if synthesis || approximal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            x_n = f1.prox(2*sol-u_n,param.gamma);
            u_n = u_n + param.lambda*(x_n - sol);
            sol = f2.prox(u_n,param.gamma);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
            else
            
            CP.K = Psi;
            CP.K_adj = Psit;
            switch thres
                case {'lasso','ew'}
                    CP.prox_f = @(x,T) myproxl1(x,tau*T,paraml1);
                otherwise
                    CP.prox_f = f1.prox;
            end
            CP.prox_g = f2.prox;
            
            % CP_solver.sigma multiplied by 10 for the second experiment,
            % in order not to stuck at zero solution
            CP_solver.sigma = 10*tau*param.gamma;
            CP_solver.tau = 1/CP_solver.sigma;
            CP_solver.theta = 1;
            
            arg = y + CP_solver.sigma*CP.K(x_bar);
            y = arg - CP_solver.sigma*CP.prox_f(arg/CP_solver.sigma, 1/CP_solver.sigma);
            x_n = CP.prox_g(x_old - CP_solver.tau*CP.K_adj(y), CP_solver.tau);
            x_bar = x_n + CP_solver.theta*(x_n-x_old);
            x_old = x_n;
            sol = x_n;

            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 
        case 'forward-backward'
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if synthesis || approximal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
            x_n = f1.prox(u_n-param.gamma*f2.grad(u_n), param.gamma);
            tn1 = (1 + sqrt(1 + 4*tn^2))/2;
            u_n=x_n+(tn-1)/tn1*(x_n-sol);
            
            % updates
            tn = tn1;
            sol = x_n;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
            else
            
            CP.f = @(x) norm(x,1);
            CP.g = @(x) 1/(2*param.gamma)*norm(u_n-param.gamma*f2.grad(u_n) - x)^2;
            CP.K = Psi;
            CP.K_adj = Psit;
            switch thres
                case {'lasso','ew'}
                    CP.prox_f = @(x,T) myproxl1(x,tau*T,paraml1);
                otherwise
                    CP.prox_f = f1.prox;
            end
            CP.prox_g = @(y,T) 1/(1+T/param.gamma)*(y + T/param.gamma * (u_n-param.gamma*f2.grad(u_n)));
            CP.dim = length(orig_sig);
            
            CP_solver.sigma = 1;
            CP_solver.tau = 1;
            CP_solver.theta = 1;
            CP_solver.maxit = 25;
            
            [ x_n, ~, ~ ] = ChambollePock(CP,CP_solver);
            tn1 = (1 + sqrt(1 + 4*tn^2))/2;
            u_n=x_n+(tn-1)/tn1*(x_n-sol);
            
            %updates
            tn = tn1;
            sol = x_n;
            
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
       
    end
    
    % plotting
    if plotting == 1
        figure(fignum);
        subplot(133);
        switch trans
            case 'sqrt'
                plotfilterbank(mat2c(sol,a(:,2)),a,fc,fs,dynrange,'nocolorbar');ylim([0,size(a,1)/2]);
            case 'mywav'
                plotcoeff(mat2c(sol,L./am),L,fs,dynrange);
            case 'gab'
                if synthesis
                    % synthesis version:
                    plotdgtreal(sol,a,M,fs,dynrange,'nocolorbar');
                else
                    % analysis version:
                    plotdgtreal(Psi(sol),a,M,fs,dynrange,'nocolorbar');
                end
            case {'erb', 'wav'}
                plotnsgtf(mat2c(sol,M),shift,fs,2,dynrange);
                % plotnsgtf(sol,shift,fs,2,dynrange);
        end     
        title('Reconstruction...');              
        pause(eps);
    elseif plotting == 2
        if isfield(settings,'snr_m')
            fprintf('i: %3d   rel_norm: %6.4f - SNR:%5.2f\n',iter,rel_norm,settings.snr_m(real(Psit(sol))));
        else
            fprintf('i: %3d   rel_norm: %6.4f\n',iter,rel_norm);
        end
    end

    % convergence test
    rel_norm = norm(sol(:) - prev_sol(:))/norm(sol(:));
    rel_norm_vec(iter) = rel_norm;
    if ((abs(rel_norm) < param.tol) && iter > 2) 
        break;
    elseif (iter >= param.maxit)
        break;
    end

    if verbose
        obj(iter) = objnorm(sol);
        fprintf('norm = %1.6f\n', obj(iter));
    end
    
    prev_sol = sol;
    iter = iter + 1;

end
    
if plotting
    dynrange = 60;
    figure(fignum);
    
    subplot(133);
    switch trans
        case 'sqrt'
            plotfilterbank(mat2c(sol,a(:,2)),a,fc,fs,dynrange,'nocolorbar');ylim([0,size(a,1)/2]);
        case 'mywav'
            plotcoeff(mat2c(sol,L./am),L,fs,dynrange);
        case 'gab'
            if synthesis
                % synthesis version:
                plotdgtreal(sol,a,M,fs,dynrange,'nocolorbar');
            else
                % analysis version:
                plotdgtreal(Psi(sol),a,M,fs,dynrange,'nocolorbar');
            end
        case {'erb', 'wav'}
            plotnsgtf(mat2c(sol,M),shift,fs,2,dynrange);
    end     
    title('Reconstruction...');
end
    
rel_norm_vec(isnan(rel_norm_vec)) = 0;
rel_norm_vec = rel_norm_vec(logical(rel_norm_vec));
if synthesis
    % synthesis version:
    sol = real(Psit(sol));
end