clc;
clear;
close all;

array = [-2 0 0; -1 0 0; 0 0 0; 1 0 0; 2 0 0];
directions = [30, 0; 35, 0; 90, 0];
Z = patten(array);
plot2d3d(Z,[0:180],0,'gain in dB','initial pattern');

%% Task 2
S = spv(array, directions);
S_mu = mean(S);
S_mean_removal = S - S_mu;
Rmm = S_mean_removal' * S_mean_removal / (size(S,1));
Rmm_uncorrelated = eye(size(S,2));

% 40 dB additive isotropic noise
sigma2=0.0001;

Rxx_theoretical = S * Rmm_uncorrelated * S' + sigma2 * eye(5, 5);

%% Task 3
load Xaudio;  % X_au
load Ximage;  % X_im

soundsc(real(X_au(2,:)), 11025);
displayimage(X_im(2,:),image_size, 201, 'The received signal at the 2nd antenna');

Rxx_au= X_au* X_au'/length( X_au(1,:))
Rxx_im= X_im* X_im'/length( X_im(1,:))

%% Task 4
% Forget 3 source presentation as well as their directions and power
directions = [];
Rmm = [];
S = [];
sigma2 = [];

%% Task 5 - Detection Problem
% only Rxx,theoretical is known
Rxx_Eigen = eig(Rxx_theoretical)
Rxx_au_Eigen = eig(Rxx_au)
Rxx_im_Eigen = eig(Rxx_im)

%% Task 6 - Estimation Problem - Conventional Approach
% The source at 90 degree is desired and known
% Two sources at 30 and 35 degree are interference and unknown

% The desired source SPV
Sd = spv(array, [90, 0]);

%  optimum weight by Wiener-Hopf solution
gain_factor = 1;
wopt = gain_factor * Rxx_theoretical \ Sd

Z=patten(array, wopt);
plot2d3d(Z, [0:180], 0, 'gain in dB', 'W-H array pattern')

%% Task 7
% 10 dB additive isotropic noise
sigma2=0.1;
directions = [30, 0; 35, 0; 90, 0];
% Step 2
S = spv(array, directions);
Rmm_uncorrelated = eye(size(S,2));
Rxx_theoretical = S * Rmm_uncorrelated * S' + sigma2 * eye(5, 5);
% Setp 4
directions = [];
Rmm = [];
S = [];
sigma2 = [];
% step 5.1
Rxx_Eigen = eig(Rxx_theoretical)
% step 6
% The desired source SPV
Sd = spv(array, [90, 0]);

%  optimum weight by Wiener-Hopf solution
gain_factor = 1;
wopt = gain_factor * Rxx_theoretical \ Sd

Z=patten(array, wopt);
plot2d3d(Z, [0:180], 0, 'gain in dB', 'W-H array pattern')

%% Task 9 MUSIC
numofsources = 3;
Z=music(array, Rxx_theoretical, numofsources);
plot2d3d(Z,[0:180],0,'dB', 'MuSIC spectrum');

%% Task 10 MUSIC - Rxx_au, Rxx_im
numofsources = 3;
Z=music(array, Rxx_au, numofsources);
plot2d3d(Z,[0:180],0,'dB', 'MuSIC spectrum (R_{xx-au})');

numofsources = 3;
Z=music(array, Rxx_im, numofsources);
plot2d3d(Z,[0:180],0,'dB', 'MuSIC spectrum (R_{xx-im})');

%% Task 11 MUSIC - Spatial Smoothing 
Rmm_coherent = [1 1 0; 1 1 0; 0 0 1];
directions = [30, 0; 35, 0; 90, 0];
S = spv(array, directions);
sigma2 = 0.0001;
Rxx_theoretical = S * Rmm_coherent * S' + sigma2 * eye(5, 5);
Rxx_Eigen = eig(Rxx_theoretical)

Z=music(array, Rxx_theoretical, numofsources);
plot2d3d(Z,[0:180],0,'dB', 'MuSIC spectrum (Before Spatial Smoothing)');
%%

% % spatial
% % M = L - p + 1 is the number of subarrays. where L is number of sensors - 5, 
% % p is the number of sources - 3.
% M = length(array) - numofsources + 1;
% % The number of independent signals
% sizeSubarray = length(array) - M + 1;
% subarray = cell(M, 1);
% spvSubarray = cell(M, 1);
% covSum = 0;
% for ii = 1: M
%     subarray{ii} = array(ii: (ii + sizeSubarray - 1), :);
%     spvSubarray{ii} = spv(subarray{ii}, directions);
%     covSum = covSum + spvSubarray{ii} * Rmm_coherent * spvSubarray{ii}';
% end
% Rmm_sp_smoothed = covSum / M;
% Rxx_sp_smoothed = S * Rmm_sp_smoothed * S' + sigma2 * eye(5, 5);
% Z=music(array, Rxx_sp_smoothed, numofsources);
% plot2d3d(Z,[0:180],0,'dB', 'MuSIC spectrum (R_{xx-im})');

Rmm_sp_smoothed = spsmooth(Rxx_theoretical,3)
Rxx_sp_smoothed = S * Rmm_sp_smoothed * S' + sigma2 * eye(5, 5);
Z=music(array, Rxx_sp_smoothed, numofsources);
plot2d3d(Z,[0:180],0,'dB', 'MuSIC spectrum (After Spatial Smoothing)');

%% Task 12
%% 12.1
% The desired source SPV
Sd = spv(array, [90, 0]);

%  optimum weight by Wiener-Hopf solution
gain_factor = 1;
wopt_au = gain_factor * Rxx_au \ Sd

Z=patten(array, wopt_au);
plot2d3d(Z, [0:180], 0, 'gain in dB', 'W-H array pattern')

yt_au = wopt_au' * X_au;
soundsc(real(yt_au), 11025);

%%  12.2
gain_factor = 1;
Sd = spv(array, [90, 0]);
wopt_im = gain_factor * Rxx_im \ Sd
yt=wopt_im'* X_im; 
displayimage(yt, image_size, 202,'The received signal at o/p of W-H beamformer');

%% 12.3 construction of a multi-beam beamformer receiving 3 signals
S = spv(array, directions);
Rmm_uncorrelated = eye(size(S,2));
Rxx_theoretical = S * Rmm_uncorrelated * S' + sigma2 * eye(5, 5);
Z = music(array, Rxx_theoretical, numofsources);
plot2d3d(Z, [0:180], 0, 'gain in dB', 'W-H array pattern')
[~, azimuth_angles] = maxk(Z, numofsources);
doa = [(azimuth_angles-1)', zeros(size(azimuth_angles'))]


wsup = superresolution(array, doa(1,:), doa);
Z=patten(array, wsup);
plot2d3d(Z, [0:180], 0, 'gain in dB', "Superresolution Beamformer array pattern (Desired direction - " + doa(1,1) + " degree)");
%%
wsup = superresolution(array, doa(2,:), doa);
Z=patten(array, wsup);
plot2d3d(Z, [0:180], 0, 'gain in dB', "Superresolution Beamformer array pattern (Desired direction - " + doa(2,1) + " degree)");
%%
wsup = superresolution(array, doa(3,:), doa);
Z=patten(array, wsup);
plot2d3d(Z, [0:180], 0, 'gain in dB', "Superresolution Beamformer array pattern (Desired direction - " + doa(3,1) + " degree)");
%%
for ii = 1:size(doa,1)
    wsup = superresolution(array, doa(ii,:), doa);
    Z=patten(array, wsup);
    plot2d3d(Z, [0:180], 0, 'gain in dB', 'Superresolution Beamformer array pattern');
    hold on;
end

%% 12.4 Rxx_im
% 90 degree
wsup = superresolution(array, doa(1,:), doa);
yt=wsup'* X_im; 
displayimage(yt/max(max(yt))*255, image_size, 124,"The received signal at o/p of Superresolution Beamformer (Desired - " + doa(1,1) + " degree)");

% 30 degree
wsup = superresolution(array, doa(2,:), doa);
yt=wsup'* X_im; 
displayimage(yt/max(max(yt))*255, image_size, 125,"The received signal at o/p of Superresolution Beamformer (Desired - " + doa(2,1) + " degree)");

% 35 degree
wsup = superresolution(array, doa(3,:), doa);
yt=wsup'* X_im; 
displayimage(yt/max(max(yt))*255, image_size, 126,"The received signal at o/p of Superresolution Beamformer (Desired - " + doa(3,1) + " degree)");

%% Task 13
num_of_snapshots = 250;
N = size(array, 1);     %N = 5
[eigenVec,eigenVal] = eig(Rxx_theoretical);


% z_snapshots = randn(num_of_snapshots,N) + 1i * randn(num_of_snapshots,N);
% for n = 1:num_of_snapshots
%     X_snapshots(:,(n-1)*length(array)+1:n*length(array)) =  eigenVec * eigenVal^(1/2) * z_snapshots(n);
% end

z_snapshots = randn(N, num_of_snapshots) + 1i * randn(N, num_of_snapshots);
X_snapshots =  eigenVec * eigenVal^(1/2) * z_snapshots;
Rxx_snap = X_snapshots * X_snapshots' / num_of_snapshots
%% AIC and MDL 
[eigenVec,eigenVal] = eig(Rxx_snap);
eigenVal = sort(abs(diag(eigenVal)), 'descend');
N = length(eigenVal);           L = num_of_snapshots;
aic = double(zeros(1, N));      mdl = double(zeros(1, N));

for k = 0:N-1
    geometric_mean = 1;
    arithmetic_mean = 0;
    
    for ll = k+1:N
        d_l = eigenVal(ll);
        geometric_mean = geometric_mean * d_l ^ (1/(N-k));
        arithmetic_mean = arithmetic_mean + d_l;
    end
    arithmetic_mean = arithmetic_mean / (N-k);

    aic(1, k+1) = -2 * ((N-k)*L) * log((geometric_mean / arithmetic_mean))...
        + 2*k*(2*N-k);
    mdl(1, k+1) = -1 * ((N-k)*L) * log((geometric_mean / arithmetic_mean))...
        + 0.5*k*(2*N-k)*log(L);
%     aic(1, k+1) = -2 * log((geometric_mean / arithmetic_mean)^((N-k)*L))...
%         + 2*k*(2*N-k);
%     mdl(1, k+1) = -log((geometric_mean / arithmetic_mean)^((N-k)*L))...
%         + 0.5*k*(2*N-k)*log(L);
end

aic
mdl

%% Task 14 adaptive beamformer (LMS and RLS)
%% LMS
array = [-2 0 0; -1 0 0; 0 0 0; 1 0 0; 2 0 0];
directions = [30, 0; 35, 0; 90, 0];
S = spv(array, directions);

num_of_iters = 2000;
weights_lms = ones(5,1);
Z_desired = [1 0 0]; % 30 degree
% Z_desired = [0 1 0];  % 35 degree
% Z_desired = [0 0 1];  % 90 degree
step_size = 0.1;
Z_error_last = Z_desired - weights_lms' * S;

for i = 1:num_of_iters
     Z_l = weights_lms' * S;
     Z_error = Z_desired - Z_l;
     disp("Iter: " + i + "; Error Norm: " + norm(Z_error));
     if (norm(Z_error) > norm(Z_error_last))
         break;
     end

     weights_lms = weights_lms + (step_size * S * conj(Z_error'));
     Z_error_last = Z_error;
end

Z=patten(array, weights_lms);
plot2d3d(Z, [0:180], 0, 'gain in dB', "Adaptive Beamformer array pattern (LMS)");

%% RLS
array = [-2 0 0; -1 0 0; 0 0 0; 1 0 0; 2 0 0];
directions = [30, 0; 35, 0; 90, 0];
S = spv(array, directions);

num_of_iters = 2000;
weights_rls = ones(5,1);
regularization_parameter = 0.7;
P_n = regularization_parameter * eye(size(S,1));
lambda = 0.01;    % forgetting factor

Z_desired = [1 0 0]; % 30 degree
% Z_desired = [0 1 0];  % 35 degree
% Z_desired = [0 0 1];  % 90 degree
Z_error_last = Z_desired - weights_rls' * S;

for n = 1: num_of_iters
    k_n = (P_n * S) / (lambda + S' * P_n * S);
    Z_l = weights_rls' * S;
    Z_error = Z_desired - Z_l;
    
    disp("Iter: " + n + "; Error Norm: " + norm(Z_error));
    if (norm(Z_error) > norm(Z_error_last))
        break;
    end
    
    weights_rls = weights_rls + k_n * Z_error';
    P_n = (lambda^-1) * P_n - (lambda^-1) * k_n * S' * P_n;
    Z_error_last = Z_error;
end

Z=patten(array, weights_rls);
plot2d3d(Z, [0:180], 0, 'gain in dB', "Adaptive Beamformer array pattern (RLS)");

function ws = superresolution(array, desiredDirection, detectedDirections)
    Sd = spv(array, desiredDirection);
    Sj = spv(array, setdiff(detectedDirections, desiredDirection, 'rows'));
    %finds the Projection operator
    ws = (eye(size(fpo(Sj))) - fpo(Sj)) * Sd;
end

function Z=music(array, Rxx_theoretical, numofsources)
% plot2d3d(Z,[0:180],0,'dB', 'MuSIC spectrum');

% Step 0: Assume M and array geometry are known
M = numofsources;
% array = [-2 0 0; -1 0 0; 0 0 0; 1 0 0; 2 0 0];
% Step 1: receive signal vector - N by 1 vector
% Step 2: find the covariance matrix    -   Rxx_theoretical
% Step 3: find the eigenvalue and eigenvector of Rxx
[eigenVec,eigenVal] = eig(Rxx_theoretical);
[~, pos] = sort(abs(diag(eigenVal)));
% Step 4: form the matrix Es with columns the eigevectors which correspond to
% the M largest eigenvalues
Es = eigenVec(:, pos(end-M+1: end));
% Step 5: find the arg of the M minima of the function
azimuths = 0:180;
elevation = 0;
norm_squared_cost = zeros(size(azimuths));
Ss = [];

for azimuth_i = azimuths
    S = spv(array, [azimuth_i, elevation]);
    Pn = eye(size(Es*Es')) - Es*Es';
    norm_squared_cost(azimuth_i+1) = S' * Pn * S;
    Ss = [Ss, S];
end

Z = abs(norm_squared_cost);
Z=-10*log10(Z);

[~, azimuth_angles] = mink(norm_squared_cost, M);
azimuth_angles-1
% Z = [azimuth_angles-1, elevation * ones(size(azimuth_angles))];
% Z = azimuth_angles-1
end
