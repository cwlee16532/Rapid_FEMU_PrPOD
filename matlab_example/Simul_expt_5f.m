function [H_meas, virtualData, model] = Simul_expt_5f(options)
%% Simul_expt_5f.m
% -------------------------------------------------------------------------
% Generate virtual sensor FRF data for the 5-story frame structure
% -------------------------------------------------------------------------
%
% Supported methods:
%   options.method = 'analytic'
%       Directly computes analytic displacement FRFs from M, C, and K.
%
%   options.method = 'dynamic'
%       Generates white-noise base excitation, performs Newmark-beta dynamic
%       analysis, and estimates acceleration FRFs by modalfrf. This block
%       preserves the research-code functionality that was previously
%       commented out. It requires Signal Processing Toolbox and NewmarkBeta.m.
%
% Output H_meas fields:
%   H_meas.FRF_X_E    : measured X-input FRF over selected frequency range
%   H_meas.FRF_Y_E    : measured Y-input FRF over selected frequency range
%   H_meas.freq_E     : selected frequency vector
%   H_meas.id_sensor  : global free-DOF indices of sensor channels
%   H_meas.id_X       : column indices for X-direction sensor channels
%   H_meas.id_Y       : column indices for Y-direction sensor channels
%
% Author: Chanwoo Lee
% Repository: PrPOD-ROM FEMU for Online Dynamic Analysis
% -------------------------------------------------------------------------

    if nargin < 1
        options = struct();
    end
    options = setDefaultOptions(options);

    %% -------------------------
    % Load and assemble system
    % -------------------------
    nodedata = xlsread(options.modelFile, options.nodeSheet);
    elementdata = xlsread(options.modelFile, options.elementSheet);

    elementdata = applyStoryParameters(elementdata, options.mu_true, ...
        options.storyBeamRows, options.sectionColumns);

    [M, K, free_dof, dofs_5fs] = assembleFiveStoryFrame(nodedata, elementdata);

    %% -------------------------
    % Modal analysis and damping
    % -------------------------
    [phi, ev] = eigs(K, M, options.nModesForModalAnalysis, options.modalShift);
    fn = sqrt(diag(ev)) / (2*pi);
    omega_n = sqrt(diag(ev));

    % Mass normalization for modal-participation diagnostics.
    for iMode = 1:size(phi,2)
        phi(:,iMode) = phi(:,iMode) / sqrt(phi(:,iMode)' * M * phi(:,iMode));
    end

    switch lower(options.dampingModel)
        case 'rayleigh'
            rayleighCoeff = (0.5 * [1/(fn(1)*2*pi), fn(1)*2*pi; ...
                                    1/(fn(2)*2*pi), fn(2)*2*pi]) \ [options.xi; options.xi];
            C = rayleighCoeff(1)*M + rayleighCoeff(2)*K;
        case 'modal'
            % Optional modal damping model inherited from the commented code.
            zeta = options.modalDampingRatios(:);
            if numel(zeta) < size(phi,2)
                zeta(numel(zeta)+1:size(phi,2),1) = 0;
            end
            C = sparse(size(M,1), size(M,1));
            for iMode = find(zeta(:).' > 0)
                phi_i = phi(:,iMode);
                C = C + 2*zeta(iMode)*omega_n(iMode)*((M*phi_i)*(M*phi_i)');
            end
            rayleighCoeff = [NaN; NaN];
        otherwise
            error('Unknown dampingModel: %s', options.dampingModel);
    end

    Xi = diag(phi' * C * phi) ./ diag(2*sqrt((phi' * M * phi) .* (phi' * K * phi)));

    %% -------------------------
    % Load/input direction and modal participation diagnostics
    % -------------------------
    eq_direction = buildDirectionMatrix(size(M,1), dofs_5fs);
    modalParticipation = computeModalParticipation(M, phi, fn, eq_direction, options.nMajorModesToReport);

    %% -------------------------
    % Sensor channels
    % -------------------------
    Sensor_index = find(ismember(dofs_5fs(:,1), options.sensorNodes) & ...
                        ismember(dofs_5fs(:,2), options.sensorDOFs));

    %% =========================================================
    % STEP 1) Generate virtual FRF data
    % =========================================================
    switch lower(options.method)
        case 'analytic'
            [FRF_X, FRF_Y, f, analyticData] = computeAnalyticSensorFRF(M, C, K, ...
                eq_direction, Sensor_index, options);
            virtualData.analyticData = analyticData;
            good_sensor_id_X = 1:2:size(FRF_X,2);
            good_sensor_id_Y = 2:2:size(FRF_Y,2);

        case 'dynamic'
            [FRF_X, FRF_Y, f, good_sensor_id_X, good_sensor_id_Y, dynamicData] = ...
                computeDynamicEstimatedFRF(M, C, K, eq_direction, Sensor_index, options);
            virtualData.dynamicData = dynamicData;

        otherwise
            error('Unknown options.method: %s', options.method);
    end

    f_start_id = find(f > options.frequencyRangeHz(1), 1, 'first');
    f_end_id   = find(f >= options.frequencyRangeHz(2), 1, 'first');
    if isempty(f_start_id) || isempty(f_end_id)
        error('Selected frequency range is outside the generated frequency vector.');
    end

    %% -------------------------
    % Output measured FRF structure
    % -------------------------
    H_meas = struct();
    H_meas.FRF_X_E = FRF_X(f_start_id:f_end_id, good_sensor_id_X);
    H_meas.FRF_Y_E = FRF_Y(f_start_id:f_end_id, good_sensor_id_Y);
    H_meas.freq_E = f(f_start_id:f_end_id);
    H_meas.id_sensor = Sensor_index;
    H_meas.id_X = good_sensor_id_X;
    H_meas.id_Y = good_sensor_id_Y;
    H_meas.mu_true = options.mu_true;
    H_meas.frequencyRangeHz = options.frequencyRangeHz;

    %% -------------------------
    % Model/output data
    % -------------------------
    model = struct();
    model.M = M;
    model.C = C;
    model.K = K;
    model.phi = phi;
    model.fn = fn;
    model.Xi = Xi;
    model.rayleighCoeff = rayleighCoeff;
    model.free_dof = free_dof;
    model.dofs_5fs = dofs_5fs;
    model.eq_direction = eq_direction;
    model.Sensor_index = Sensor_index;
    model.modalParticipation = modalParticipation;

    virtualData.f = f;
    virtualData.FRF_X = FRF_X;
    virtualData.FRF_Y = FRF_Y;
    virtualData.f_start_id = f_start_id;
    virtualData.f_end_id = f_end_id;
    virtualData.good_sensor_id_X = good_sensor_id_X;
    virtualData.good_sensor_id_Y = good_sensor_id_Y;

    if options.makePlots
        plotVirtualFRF(H_meas, virtualData, model);
    end
end

%% =========================================================
% Local helper functions
% =========================================================
function options = setDefaultOptions(options)
    defaults = struct();
    defaults.method = 'analytic';
    defaults.modelFile = 'model_5f';
    defaults.nodeSheet = 'Sheet2';
    defaults.elementSheet = 'Sheet1';
    defaults.mu_true = [0.004 0.005 0.006 0.007 0.008];
    defaults.storyBeamRows = {9:10; 25:26; 41:42; 57:58; 73:74};
    defaults.sectionColumns = 5:6;
    defaults.nModesForModalAnalysis = 50;
    defaults.modalShift = 0.01;
    defaults.xi = 0.01;
    defaults.dampingModel = 'rayleigh';
    defaults.modalDampingRatios = 0.01*ones(10,1);
    defaults.sensorNodes = [105, 113, 121, 129, 137, 108, 116, 124, 132, 140];
    defaults.sensorDOFs = [1 2];
    defaults.frequencyVectorHz = (0.05:0.05:30).';
    defaults.frequencyRangeHz = [3 30];
    defaults.addNoise = false;
    defaults.noiseLevelPercent = 1;
    defaults.makePlots = false;
    defaults.nMajorModesToReport = 10;

    % Dynamic FRF-estimation settings.
    defaults.whiteNoiseSamples = 120000;
    defaults.whiteNoisePower = 10;
    defaults.dt = 0.001;
    defaults.fftBlockSize = 2^14;
    defaults.overlapRatio = 0.677;
    defaults.coherenceThreshold = 0.90;
    defaults.frfEstimator = 'H2';

    names = fieldnames(defaults);
    for i = 1:numel(names)
        if ~isfield(options, names{i}) || isempty(options.(names{i}))
            options.(names{i}) = defaults.(names{i});
        end
    end
end

function [FRF_X, FRF_Y, f, analyticData] = computeAnalyticSensorFRF(M, C, K, eq_direction, Sensor_index, options)
    f = options.frequencyVectorHz(:);
    nW = numel(f);
    nDOF = size(M,1);

    FRF_X_full = zeros(nDOF, nW);
    FRF_Y_full = zeros(nDOF, nW);
    FRF_XY_full = zeros(nDOF, nW);

    GF_X = eq_direction(:,1);
    GF_Y = eq_direction(:,2);
    GF_XY = GF_X + GF_Y;

    for iw = 1:nW
        omega = 2*pi*f(iw);
        A = K + 1i*omega*C - omega^2*M;
        FRF_X_full(:,iw) = A \ GF_X;
        FRF_Y_full(:,iw) = A \ GF_Y;
        FRF_XY_full(:,iw) = A \ GF_XY;
    end

    if options.addNoise
        FRF_X_full = addComplexNoise(FRF_X_full, options.noiseLevelPercent);
        FRF_Y_full = addComplexNoise(FRF_Y_full, options.noiseLevelPercent);
        FRF_XY_full = addComplexNoise(FRF_XY_full, options.noiseLevelPercent);
    end

    FRF_X = FRF_X_full(Sensor_index,:).';
    FRF_Y = FRF_Y_full(Sensor_index,:).';

    analyticData = struct();
    analyticData.FRF_X_full = FRF_X_full;
    analyticData.FRF_Y_full = FRF_Y_full;
    analyticData.FRF_XY_full = FRF_XY_full;
end

function [FRF_X, FRF_Y, f, good_sensor_id_X, good_sensor_id_Y, dynamicData] = ...
    computeDynamicEstimatedFRF(M, C, K, eq_direction, Sensor_index, options)
% Preserves the original commented virtual-experiment workflow.

    if exist('wgn', 'file') ~= 2
        error('The dynamic method requires wgn from the Communications Toolbox or a replacement noise generator.');
    end
    if exist('modalfrf', 'file') ~= 2
        error('The dynamic method requires modalfrf from the Signal Processing Toolbox.');
    end
    if exist('NewmarkBeta', 'file') ~= 2
        error('The dynamic method requires NewmarkBeta.m.');
    end

    GF_white = wgn(2, options.whiteNoiseSamples, options.whiteNoisePower, 'linear');
    GF_white_base = -(M) * (eq_direction(:,1)*GF_white(1,:) + eq_direction(:,2)*GF_white(2,:));

    Nmk = struct();
    Nmk.nStp = size(GF_white_base,2)-1;
    Nmk.dt = options.dt;
    Nmk.t = (0:Nmk.dt:Nmk.dt*Nmk.nStp).';

    nDOF = size(M,1);
    GU = zeros(nDOF, Nmk.nStp+1);
    GV = zeros(nDOF, Nmk.nStp+1);
    GA = zeros(nDOF, Nmk.nStp+1);

    [GU, GV, GA] = NewmarkBeta(M, C, K, GU, GV, GA, GF_white_base, Nmk); %#ok<ASGLU>

    GA_sensor = GA(Sensor_index,:);
    if options.addNoise
        GA_sensor_noise = addRealNoise(GA_sensor, options.noiseLevelPercent);
    else
        GA_sensor_noise = GA_sensor;
    end

    Nfft = options.fftBlockSize;
    Fs = 1 / Nmk.dt;
    window = hanning(Nfft);
    noverlap = round(Nfft*options.overlapRatio);

    FRF_X = zeros(Nfft/2+1, size(GA_sensor_noise,1));
    FRF_Y = zeros(Nfft/2+1, size(GA_sensor_noise,1));
    coh_X = zeros(Nfft/2+1, size(GA_sensor_noise,1));
    coh_Y = zeros(Nfft/2+1, size(GA_sensor_noise,1));

    for i = 1:size(GA_sensor_noise,1)
        [FRF_X(:,i), f, coh_X(:,i)] = modalfrf(-GF_white(1,:).', GA_sensor_noise(i,:).', ...
            Fs, window, noverlap, 'Sensor', 'acc', 'Estimator', options.frfEstimator);
        [FRF_Y(:,i), ~, coh_Y(:,i)] = modalfrf(-GF_white(2,:).', GA_sensor_noise(i,:).', ...
            Fs, window, noverlap, 'Sensor', 'acc', 'Estimator', options.frfEstimator);
    end

    % Scaling retained from the original research script.
    GF_X = eq_direction(:,1);
    GF_Y = eq_direction(:,2);
    FRF_X = FRF_X * sum(GF_X) / sum(M*GF_X);
    FRF_Y = FRF_Y * sum(GF_Y) / sum(M*GF_Y);

    f_start_id = find(f > options.frequencyRangeHz(1), 1, 'first');
    f_end_id = find(f >= options.frequencyRangeHz(2), 1, 'first');
    good_sensor_id_X = find(mean(coh_X(f_start_id:f_end_id,:),1) > options.coherenceThreshold);
    good_sensor_id_Y = find(mean(coh_Y(f_start_id:f_end_id,:),1) > options.coherenceThreshold);

    dynamicData = struct();
    dynamicData.GF_white = GF_white;
    dynamicData.Nmk = Nmk;
    dynamicData.GA_sensor = GA_sensor;
    dynamicData.GA_sensor_noise = GA_sensor_noise;
    dynamicData.coh_X = coh_X;
    dynamicData.coh_Y = coh_Y;
end

function modalParticipation = computeModalParticipation(M, phi, fn, eq_direction, nMajorModes)
    nMode = size(phi,2);
    nDir = size(eq_direction,2);

    L = zeros(nMode,nDir);
    gamma = zeros(nMode,nDir);
    effm = zeros(nMode,nDir);
    Rm = zeros(nMode,nDir);
    total_mass = zeros(1,nDir);
    major_mode_sorted = zeros(nMajorModes,nDir);
    major_mode_freq = zeros(nMajorModes,nDir);
    major_mode_Rm = zeros(nMajorModes,nDir);

    modalMass = diag(phi' * M * phi);
    for iDir = 1:nDir
        L(:,iDir) = phi' * M * eq_direction(:,iDir);
        gamma(:,iDir) = L(:,iDir) ./ modalMass;
        effm(:,iDir) = L(:,iDir) .* gamma(:,iDir);
        total_mass(iDir) = eq_direction(:,iDir)' * M * eq_direction(:,iDir);
        Rm(:,iDir) = effm(:,iDir) / total_mass(iDir);
        [effm_sorted, idx] = sort(effm(:,iDir), 'descend'); %#ok<ASGLU>
        nReport = min(nMajorModes, numel(idx));
        major_mode_sorted(1:nReport,iDir) = idx(1:nReport);
        major_mode_freq(1:nReport,iDir) = fn(idx(1:nReport));
        major_mode_Rm(1:nReport,iDir) = Rm(idx(1:nReport),iDir);
    end

    modalParticipation = struct();
    modalParticipation.L = L;
    modalParticipation.gamma = gamma;
    modalParticipation.effm = effm;
    modalParticipation.Rm = Rm;
    modalParticipation.total_mass = total_mass;
    modalParticipation.cumulative_Rm = sum(Rm,1);
    modalParticipation.major_mode_sorted = major_mode_sorted;
    modalParticipation.major_mode_freq = major_mode_freq;
    modalParticipation.major_mode_Rm = major_mode_Rm;
end

function Xn = addComplexNoise(X, noiseLevelPercent)
    noiseScale = noiseLevelPercent / 100;
    Xn = X + noiseScale * rms(abs(X(:))) * (randn(size(X)) + 1i*randn(size(X))) / sqrt(2);
end

function Xn = addRealNoise(X, noiseLevelPercent)
    noiseScale = noiseLevelPercent / 100;
    Xn = X + noiseScale * rms(X(:)) * randn(size(X));
end

function elementdata = applyStoryParameters(elementdata, mu, storyBeamRows, sectionColumns)
    for k = 1:numel(storyBeamRows)
        elementdata(storyBeamRows{k}, sectionColumns) = mu(k);
    end
end

function [M, K, free_dof, dofs_5fs] = assembleFiveStoryFrame(nodedata, elementdata)
    N = size(nodedata,1);
    El = size(elementdata,1);
    Fixity = nodedata(:,5:10);
    bounded_dof = find(Fixity' ~= 0);
    free_dof = setdiff(1:6*N, bounded_dof);

    dofs_5fs = zeros(N*6,2);
    for iNode = 1:N
        dofs_5fs(6*iNode-5:6*iNode,1) = iNode;
        dofs_5fs(6*iNode-5:6*iNode,2) = (1:6).';
    end
    dofs_5fs(bounded_dof,:) = [];

    KGlobal = zeros(6*N,6*N);
    MGlobal = zeros(6*N,6*N);

    for iEl = 1:El
        elementdata_i = elementdata(iEl,:);
        b = elementdata_i(5);
        d = elementdata_i(6);
        E = elementdata_i(7);
        G = elementdata_i(8);
        Iz = b*d^3/12;
        Iy = d*b^3/12;
        J = Iz + Iy;
        rho = elementdata_i(17);

        XYZ1 = nodedata(elementdata_i(2),2:4);
        XYZ2 = nodedata(elementdata_i(3),2:4);
        x1 = XYZ1(1); y1 = XYZ1(2); z1 = XYZ1(3);
        x2 = XYZ2(1); y2 = XYZ2(2); z2 = XYZ2(3);

        [stiff_local, mass_local] = StiffMass(E, G, d, b, Iy, Iz, J, rho, x1, y1, z1, x2, y2, z2);
        KGlobal = stiffassemble(KGlobal, stiff_local, elementdata_i(2), elementdata_i(3));
        MGlobal = massassemble(MGlobal, mass_local, elementdata_i(2), elementdata_i(3));
    end

    K = KGlobal(free_dof, free_dof);
    M = MGlobal(free_dof, free_dof);
end

function eq_direction = buildDirectionMatrix(nFreeDOF, dofs_5fs)
    eq_direction = zeros(nFreeDOF, 6);
    for iDir = 1:6
        dof_index = find(dofs_5fs(:,2) == iDir);
        eq_direction(dof_index, iDir) = 1;
    end
end

function plotVirtualFRF(H_meas, virtualData, model)
    figure('Name','Virtual FRF data','Position',[100 100 800 450]);
    for i = 1:min(10, numel(H_meas.id_X))
        subplot(5,2,i);
        semilogy(H_meas.freq_E, abs(H_meas.FRF_X_E(:,i)));
        hold on; grid on; box on;
        xlabel('Frequency (Hz)');
        ylabel('|H(f)|');
        title(sprintf('X FRF channel %d', i));
    end

    if isfield(model, 'fn')
        f_lines = model.fn(1:min(5,numel(model.fn)));
        axesList = findall(gcf,'Type','axes');
        for ia = 1:numel(axesList)
            axes(axesList(ia)); %#ok<LAXES>
            yLimits = ylim;
            for k = 1:numel(f_lines)
                line([f_lines(k) f_lines(k)], yLimits, 'Color', [0 0 0], ...
                    'LineStyle', '--', 'HandleVisibility','off');
            end
        end
    end

    fprintf('Selected frequency range: %.3f--%.3f Hz\n', ...
        virtualData.f(virtualData.f_start_id), virtualData.f(virtualData.f_end_id));
end
