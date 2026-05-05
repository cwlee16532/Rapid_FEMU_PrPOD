%% main_online_time_response.m
% -------------------------------------------------------------------------
% Online phase: fast time-domain response analysis using updated PrPOD-ROM
% -------------------------------------------------------------------------
% This script performs the final online dynamic response analysis after
% PrPOD-ROM-based FEMU. The updated parameter vector mu_star identified by
% main_online_FEMU_PSO.m is used to assemble the updated reduced-order model
% and compute the time-domain response efficiently.
%
% Main procedure:
%   1) Load FEMU result file containing mu_star, M_final, ROM, and options
%   2) Interpolate the congruence-aligned POD basis at mu_star
%   3) Assemble M(mu_star), C(mu_star), and K(mu_star)
%   4) Project full-order operators by Galerkin projection
%   5) Generate default white-noise base acceleration or read user input
%   6) Solve the reduced-order dynamic equation by Newmark-beta integration
%   7) Reconstruct full-order displacement, velocity, and acceleration
%
% Paper notation used in variable names:
%   mu_star       : updated parameter vector from FEMU
%   Phi_tilde     : interpolated congruence-aligned POD basis
%   M_r, C_r, K_r : reduced mass, damping, and stiffness matrices
%   u_hat_r       : reduced displacement response
%   u_hat         : reconstructed full-order displacement response
%   r, M_final    : retained and final POD ranks
%
% Required files/functions:
%   - FEMU_PSO_5F_PrPOD_result.mat
%   - model_5f.xlsx
%   - StiffMass.m, stiffassemble.m, massassemble.m
%
% Optional user input files:
%   - A .mat file containing either
%       input.time, input.acceleration
%     or variables named time and acceleration.
%   - acceleration must be nStep x 1, nStep x 2, or nStep x 3.
%     Columns are interpreted as X, Y, and Z base acceleration components.
%
% Author: Chanwoo Lee
% Repository: PrPOD-ROM FEMU for Online Dynamic Analysis
% -------------------------------------------------------------------------

clear; close all; clc;

%% -------------------------
% User-defined paths
% -------------------------
femuResultFile = 'FEMU_PSO_5F_PrPOD_result.mat';
timeResponseResultFile = 'time_response_5F_PrPOD_result.mat';

% These values are used if the FEMU result file does not contain them.
modelFile = 'model_5f';
nodeSheet = 'Sheet2';
elementSheet = 'Sheet1';

%% -------------------------
% Time-response analysis settings
% -------------------------
timeOptions = struct();
timeOptions.responseRank = 300;       % [] uses M_final from the FEMU result
timeOptions.analysisType = 'ROM_with_FOM_check';    % 'ROM' or 'ROM_with_FOM_check'
timeOptions.makePlots = true;
% Default time window for ROM-FOM comparison plot [t_min t_max] in seconds.
% Use [] to show the full duration.
% The ROM-FOM validation plot compares relative displacement by default.
timeOptions.plotTimeWindow = [0 10];
timeOptions.saveFullResponse = true;
timeOptions.outputSensorNodes = [105, 113, 121, 129, 137, 108, 116, 124, 132, 140];
timeOptions.outputSensorDOFs = [1 2];

% Newmark-beta settings. beta = 1/4 and gamma = 1/2 give the average
% acceleration method for linear systems.
timeOptions.newmark.beta = 1/4;
timeOptions.newmark.gamma = 1/2;

%% -------------------------
% External input settings
% -------------------------
% inputType options:
%   'white_noise_base_accel' : default broadband base acceleration
%   'user_base_accel'        : user-defined base acceleration time history
%   'user_force_matrix'      : user-defined physical force matrix F(t)
%   'function_handle'        : input generated from a MATLAB function handle
inputOptions = struct();
inputOptions.inputType = 'white_noise_base_accel';

% Default white-noise base-acceleration input.
inputOptions.dt = 0.01;
inputOptions.duration = 30.0;
inputOptions.randomSeed = 1;
inputOptions.whiteNoiseRMS = 1.0;
inputOptions.whiteNoiseDirections = [1 2];  % X and Y base acceleration

% Optional user-defined input file.
% Example file format:
%   input.time = (0:0.001:30).';
%   input.acceleration = [ax ay];
inputOptions.userInputFile = '';

% Optional direct user-defined arrays. These override userInputFile when
% inputType = 'user_base_accel' or 'user_force_matrix'.
inputOptions.userTime = [];
inputOptions.userBaseAcceleration = [];  % nStep x nDir, columns X/Y/Z
inputOptions.userForceMatrix = [];       % nDOF x nStep

% Optional function handle. It must return [time, inputData].
% For inputType = 'function_handle', define for example:
% inputOptions.inputFunction = @(opt) myInputFunction(opt);
inputOptions.inputFunction = [];

%% =========================================================
% STEP 1) Load FEMU result
% =========================================================
fprintf('Loading FEMU result file: %s\n', femuResultFile);
R = load(femuResultFile);

if ~isfield(R, 'mu_star')
    error('The FEMU result file must contain mu_star.');
end
mu_star = R.mu_star(:).';

% Exact parameter for full-order reference analysis.
% In the FEMU result file, mu_exact may already be saved. If not, this
% script also supports a separate file named mu_exact.mat.
if isfield(R, 'mu_exact')
    mu_exact = R.mu_exact(:).';
elseif exist('mu_exact.mat', 'file') == 2
    S_exact = load('mu_exact.mat');
    if isfield(S_exact, 'mu_exact')
        mu_exact = S_exact.mu_exact(:).';
    else
        error('mu_exact.mat must contain the variable mu_exact.');
    end
else
    warning('mu_exact was not found. Full-order check will use mu_star.');
    mu_exact = mu_star;
end

if isfield(R, 'M_final') && isempty(timeOptions.responseRank)
    r = R.M_final;
elseif isempty(timeOptions.responseRank)
    r = 300;
else
    r = timeOptions.responseRank;
end
M_final = r;

if isfield(R, 'ROM')
    ROM = R.ROM;
else
    error('The FEMU result file must contain ROM.');
end

if isfield(R, 'objectiveOptions') && isfield(R.objectiveOptions, 'romOptions')
    romOptions = R.objectiveOptions.romOptions;
else
    romOptions = struct();
end
romOptions = setDefaultROMOptions(romOptions, modelFile, nodeSheet, elementSheet);

fprintf('Updated parameter mu_star = %s\n', mat2str(mu_star, 8));
fprintf('Exact parameter mu_exact = %s\n', mat2str(mu_exact, 8));
fprintf('Time-response ROM rank r = %d\n', r);

%% =========================================================
% STEP 2) Assemble updated reduced-order model
% =========================================================
fprintf('Assembling updated PrPOD-ROM...\n');
tic;
[updatedROM, updatedModel] = assembleUpdatedPrPODROM(mu_star, ROM, r, romOptions);
fprintf('Updated ROM assembled in %.3f s.\n', toc);

%% =========================================================
% STEP 3) Generate external input
% =========================================================
[inputData, inputOptions] = generateDynamicInput(inputOptions, updatedModel);

fprintf('Input type: %s\n', inputOptions.inputType);
fprintf('Time step dt = %.6f s, nStep = %d, duration = %.3f s\n', ...
    inputData.dt, inputData.nStep, inputData.time(end));

%% =========================================================
% STEP 4) Reduced-order dynamic analysis
% =========================================================
fprintf('Running reduced-order Newmark-beta dynamic analysis...\n');
tic;

F_r_time = updatedROM.Phi_tilde' * inputData.forceMatrix;

[u_hat_r, v_hat_r, a_hat_r] = newmarkBetaLinear( ...
    updatedROM.M_r, updatedROM.C_r, updatedROM.K_r, ...
    F_r_time, inputData.dt, timeOptions.newmark, ...
    'ROM dynamic analysis');

elapsedROM = toc;
fprintf('ROM time integration completed in %.3f s.\n', elapsedROM);

% Reconstruct full-order responses.
% The interpolated POD basis can be complex-valued. The physical time-domain
% response is real-valued, so the real part is retained for output.
u_hat = real(updatedROM.Phi_tilde * u_hat_r);
v_hat = real(updatedROM.Phi_tilde * v_hat_r);
a_hat_relative = real(updatedROM.Phi_tilde * a_hat_r);

% Absolute acceleration is obtained by adding the imposed base acceleration
% components to the corresponding translational DOFs.
a_hat_absolute = addBaseAccelerationToRelativeResponse( ...
    a_hat_relative, inputData.baseAcceleration, updatedModel.eq_direction);

%% =========================================================
% STEP 5) Optional full-order check
% =========================================================
fullOrderResult = [];
elapsedFOM = NaN;

if strcmpi(timeOptions.analysisType, 'ROM_with_FOM_check')
    fprintf('Running optional full-order Newmark-beta check using exact parameters...\n');
    fprintf('Exact parameter mu_exact = %s\n', mat2str(mu_exact, 8));

    % Assemble the full-order reference model using exact parameters,
    % not the FEMU-updated parameters.
    [~, exactModel] = assembleUpdatedPrPODROM(mu_exact, ROM, r, romOptions);

    % Regenerate the same input type for the exact model. This is important
    % for base-acceleration input because f(t) = -M R a_g(t) depends on M.
    [inputDataExact, ~] = generateDynamicInput(inputOptions, exactModel);

    tic;
    [u_FOM, v_FOM, a_FOM_relative] = newmarkBetaLinear( ...
        exactModel.M, exactModel.C, exactModel.K, ...
        inputDataExact.forceMatrix, inputDataExact.dt, timeOptions.newmark, ...
        'Full-order dynamic analysis');
    elapsedFOM = toc;

    a_FOM_absolute = addBaseAccelerationToRelativeResponse( ...
        a_FOM_relative, inputDataExact.baseAcceleration, exactModel.eq_direction);

    fullOrderResult = struct();
    fullOrderResult.mu_exact = mu_exact;
    fullOrderResult.u = real(u_FOM);
    fullOrderResult.v = real(v_FOM);
    fullOrderResult.a_relative = real(a_FOM_relative);
    fullOrderResult.a_absolute = real(a_FOM_absolute);
    fullOrderResult.inputData = inputDataExact;
    fullOrderResult.elapsedFOM = elapsedFOM;

    fprintf('FOM time integration completed in %.3f s.\n', elapsedFOM);
end

%% =========================================================
% STEP 6) Extract selected sensor responses
% =========================================================
outputIndex = find(ismember(updatedModel.dofs_5fs(:,1), timeOptions.outputSensorNodes) & ...
                   ismember(updatedModel.dofs_5fs(:,2), timeOptions.outputSensorDOFs));

sensorResponse = struct();
sensorResponse.outputIndex = outputIndex;
sensorResponse.outputDOFs = updatedModel.dofs_5fs(outputIndex,:);
sensorResponse.u = u_hat(outputIndex,:);
sensorResponse.v = v_hat(outputIndex,:);
sensorResponse.a_relative = a_hat_relative(outputIndex,:);
sensorResponse.a_absolute = a_hat_absolute(outputIndex,:);

%% =========================================================
% STEP 7) Plot and save results
% =========================================================
if timeOptions.makePlots
    plotTimeResponse(inputData, sensorResponse, fullOrderResult, outputIndex, timeOptions.plotTimeWindow);
end

result = struct();
result.mu_star = mu_star;
result.mu_exact = mu_exact;
result.M_final = M_final;
result.r = r;
result.time = inputData.time;
result.inputData = inputData;
result.updatedROM = updatedROM;
result.updatedModel = stripLargeFields(updatedModel, timeOptions.saveFullResponse);
result.sensorResponse = sensorResponse;
result.elapsedROM = elapsedROM;
result.elapsedFOM = elapsedFOM;

if timeOptions.saveFullResponse
    result.u_hat = u_hat;
    result.v_hat = v_hat;
    result.a_hat_relative = a_hat_relative;
    result.a_hat_absolute = a_hat_absolute;
end

if ~isempty(fullOrderResult)
    result.fullOrderResult = fullOrderResult;
end

save(timeResponseResultFile, 'result', 'timeOptions', 'inputOptions', '-v7.3');

fprintf('\nOnline time-response analysis completed successfully.\n');
fprintf('ROM elapsed time: %.3f s\n', elapsedROM);
if ~isnan(elapsedFOM)
    fprintf('FOM elapsed time: %.3f s\n', elapsedFOM);
    fprintf('Speed-up ratio FOM/ROM: %.2f\n', elapsedFOM / elapsedROM);
end
fprintf('Result file: %s\n', timeResponseResultFile);

%% =========================================================
% Local helper functions
% =========================================================
function romOptions = setDefaultROMOptions(romOptions, modelFile, nodeSheet, elementSheet)
    defaults = struct();
    defaults.modelFile = modelFile;
    defaults.nodeSheet = nodeSheet;
    defaults.elementSheet = elementSheet;
    defaults.storyBeamRows = {9:10; 25:26; 41:42; 57:58; 73:74};
    defaults.sectionColumns = 5:6;
    defaults.xi = 0.01;
    defaults.dampingModel = 'rayleigh';
    defaults.modalDampingRatios = 0.01*ones(10,1);
    defaults.nModesForDamping = 5;
    defaults.modalShift = 0.01;
    defaults.podInterpolationMethod = 'linear';
    defaults.podExtrapolationMethod = 'nearest';

    names = fieldnames(defaults);
    for i = 1:numel(names)
        if ~isfield(romOptions, names{i}) || isempty(romOptions.(names{i}))
            romOptions.(names{i}) = defaults.(names{i});
        end
    end
end

function [updatedROM, model] = assembleUpdatedPrPODROM(mu_star, ROM, r, opt)
    mu_star = mu_star(:).';
    mu_star = enforceParameterBoundsIfAvailable(mu_star, ROM);

    Phi_tilde = interpolateBasisAtParameter(ROM, mu_star, r, opt);

    nodedata = xlsread(opt.modelFile, opt.nodeSheet);
    elementdata = xlsread(opt.modelFile, opt.elementSheet);
    elementdata = applyStoryParameters(elementdata, mu_star, opt.storyBeamRows, opt.sectionColumns);
    [M, K, free_dof, dofs_5fs] = assembleFiveStoryFrame(nodedata, elementdata);

    C = computeDampingMatrix(M, K, opt);
    eq_direction = buildDirectionMatrix(size(M,1), dofs_5fs);

    M_r = Phi_tilde' * M * Phi_tilde;
    C_r = Phi_tilde' * C * Phi_tilde;
    K_r = Phi_tilde' * K * Phi_tilde;

    updatedROM = struct();
    updatedROM.mu_star = mu_star;
    updatedROM.Phi_tilde = Phi_tilde;
    updatedROM.M_r = M_r;
    updatedROM.C_r = C_r;
    updatedROM.K_r = K_r;

    model = struct();
    model.M = M;
    model.C = C;
    model.K = K;
    model.free_dof = free_dof;
    model.dofs_5fs = dofs_5fs;
    model.eq_direction = eq_direction;
end

function Phi_tilde = interpolateBasisAtParameter(ROM, mu_star, r, opt)
    if isfield(ROM, 'useFastInterpolant') && ROM.useFastInterpolant
        Phi_interp_full = squeeze(ROM.POM_interpolant( ...
            mu_star(1), mu_star(2), mu_star(3), mu_star(4), mu_star(5)));

        if size(Phi_interp_full,1) < size(Phi_interp_full,2)
            Phi_interp_full = Phi_interp_full.';
        end

        r = min(r, size(Phi_interp_full,2));
        Phi_tilde = Phi_interp_full(:,1:r);
    else
        Phi_tilde = interpolateAlignedBasis(ROM.Phi_aligned_database, ROM.p_grid, ...
            mu_star, r, opt.podInterpolationMethod, opt.podExtrapolationMethod);
    end
end

function Phi_tilde = interpolateAlignedBasis(PhiDB, p_grid, mu_trial, r, interpolationMethod, extrapolationMethod)
    dbSize = size(PhiDB);
    if numel(dbSize) ~= 3
        error('Phi_aligned_database must have size Nop x nDOF x nPOM.');
    end

    nDOF = dbSize(2);
    nPOM = dbSize(3);
    r = min(r, nPOM);

    n1 = size(p_grid{1},1);
    n2 = size(p_grid{1},2);
    n3 = size(p_grid{1},3);
    n4 = size(p_grid{1},4);
    n5 = size(p_grid{1},5);

    % Interpolate all stored POMs first. The first r modes are extracted
    % after interpolation to preserve consistency with the offline database.
    Phi_reshaped = reshape(PhiDB, [n1, n2, n3, n4, n5, nDOF, nPOM]);

    Phi_interp_full = interpn(p_grid{1}, p_grid{2}, p_grid{3}, p_grid{4}, p_grid{5}, ...
        Phi_reshaped, mu_trial(1), mu_trial(2), mu_trial(3), mu_trial(4), mu_trial(5), interpolationMethod);

    if any(isnan(Phi_interp_full(:)))
        Phi_interp_full = interpn(p_grid{1}, p_grid{2}, p_grid{3}, p_grid{4}, p_grid{5}, ...
            Phi_reshaped, mu_trial(1), mu_trial(2), mu_trial(3), mu_trial(4), mu_trial(5), extrapolationMethod);
    end

    Phi_interp_full = squeeze(Phi_interp_full);
    if size(Phi_interp_full,1) ~= nDOF
        Phi_interp_full = reshape(Phi_interp_full, [nDOF, nPOM]);
    end

    Phi_tilde = Phi_interp_full(:,1:r);
end

function mu = enforceParameterBoundsIfAvailable(mu, ROM)
    if isfield(ROM, 'p_grid')
        lb = [min(ROM.p_grid{1}(:)), min(ROM.p_grid{2}(:)), min(ROM.p_grid{3}(:)), ...
              min(ROM.p_grid{4}(:)), min(ROM.p_grid{5}(:))];
        ub = [max(ROM.p_grid{1}(:)), max(ROM.p_grid{2}(:)), max(ROM.p_grid{3}(:)), ...
              max(ROM.p_grid{4}(:)), max(ROM.p_grid{5}(:))];
        mu = max(min(mu, ub), lb);
    end
end

function C = computeDampingMatrix(M, K, opt)
    switch lower(opt.dampingModel)
        case 'rayleigh'
            [~, ev] = eigs(K, M, opt.nModesForDamping, opt.modalShift);
            fn = sqrt(diag(ev)) / (2*pi);
            rayleighCoeff = (0.5 * [1/(fn(1)*2*pi), fn(1)*2*pi; ...
                                    1/(fn(2)*2*pi), fn(2)*2*pi]) \ [opt.xi; opt.xi];
            C = rayleighCoeff(1)*M + rayleighCoeff(2)*K;

        case 'modal'
            [phi, ev] = eigs(K, M, opt.nModesForDamping, opt.modalShift);
            omega_n = sqrt(diag(ev));
            for iMode = 1:size(phi,2)
                phi(:,iMode) = phi(:,iMode) / sqrt(phi(:,iMode)' * M * phi(:,iMode));
            end
            zeta = opt.modalDampingRatios(:);
            if numel(zeta) < size(phi,2)
                zeta(numel(zeta)+1:size(phi,2),1) = 0;
            end
            C = sparse(size(M,1), size(M,1));
            for iMode = find(zeta(:).' > 0)
                phi_i = phi(:,iMode);
                C = C + 2*zeta(iMode)*omega_n(iMode)*((M*phi_i)*(M*phi_i)');
            end

        otherwise
            error('Unknown dampingModel: %s', opt.dampingModel);
    end
end

function [inputData, inputOptions] = generateDynamicInput(inputOptions, model)
    switch lower(inputOptions.inputType)
        case 'white_noise_base_accel'
            rng(inputOptions.randomSeed);
            time = (0:inputOptions.dt:inputOptions.duration).';
            nStep = numel(time);
            baseAcceleration = zeros(nStep, 3);
            for iDir = inputOptions.whiteNoiseDirections(:).'
                baseAcceleration(:,iDir) = inputOptions.whiteNoiseRMS * randn(nStep,1);
            end
            forceMatrix = baseAccelerationToForceMatrix(model.M, model.eq_direction, baseAcceleration);

        case 'user_base_accel'
            [time, baseAcceleration] = readUserBaseAcceleration(inputOptions);
            inputOptions.dt = time(2) - time(1);
            forceMatrix = baseAccelerationToForceMatrix(model.M, model.eq_direction, baseAcceleration);

        case 'user_force_matrix'
            [time, forceMatrix] = readUserForceMatrix(inputOptions, size(model.M,1));
            inputOptions.dt = time(2) - time(1);
            baseAcceleration = zeros(numel(time), 3);

        case 'function_handle'
            if isempty(inputOptions.inputFunction) || ~isa(inputOptions.inputFunction, 'function_handle')
                error('inputOptions.inputFunction must be a function handle.');
            end
            [time, inputDataRaw] = inputOptions.inputFunction(inputOptions);
            time = time(:);
            inputOptions.dt = time(2) - time(1);
            if size(inputDataRaw,1) == size(model.M,1)
                forceMatrix = inputDataRaw;
                baseAcceleration = zeros(numel(time), 3);
            else
                baseAcceleration = normalizeBaseAccelerationSize(inputDataRaw, numel(time));
                forceMatrix = baseAccelerationToForceMatrix(model.M, model.eq_direction, baseAcceleration);
            end

        otherwise
            error('Unknown inputType: %s', inputOptions.inputType);
    end

    inputData = struct();
    inputData.time = time(:);
    inputData.dt = inputOptions.dt;
    inputData.nStep = numel(time);
    inputData.baseAcceleration = baseAcceleration;
    inputData.forceMatrix = forceMatrix;
end

function [time, baseAcceleration] = readUserBaseAcceleration(inputOptions)
    if ~isempty(inputOptions.userTime) && ~isempty(inputOptions.userBaseAcceleration)
        time = inputOptions.userTime(:);
        baseAcceleration = inputOptions.userBaseAcceleration;
    elseif ~isempty(inputOptions.userInputFile)
        S = load(inputOptions.userInputFile);
        if isfield(S, 'input')
            time = S.input.time(:);
            baseAcceleration = S.input.acceleration;
        elseif isfield(S, 'time') && isfield(S, 'acceleration')
            time = S.time(:);
            baseAcceleration = S.acceleration;
        else
            error('User input file must contain input.time/input.acceleration or time/acceleration.');
        end
    else
        error('No user-defined base-acceleration input was provided.');
    end
    baseAcceleration = normalizeBaseAccelerationSize(baseAcceleration, numel(time));
end

function [time, forceMatrix] = readUserForceMatrix(inputOptions, nDOF)
    if ~isempty(inputOptions.userTime) && ~isempty(inputOptions.userForceMatrix)
        time = inputOptions.userTime(:);
        forceMatrix = inputOptions.userForceMatrix;
    elseif ~isempty(inputOptions.userInputFile)
        S = load(inputOptions.userInputFile);
        if isfield(S, 'input')
            time = S.input.time(:);
            forceMatrix = S.input.forceMatrix;
        elseif isfield(S, 'time') && isfield(S, 'forceMatrix')
            time = S.time(:);
            forceMatrix = S.forceMatrix;
        else
            error('User input file must contain input.time/input.forceMatrix or time/forceMatrix.');
        end
    else
        error('No user-defined force matrix was provided.');
    end

    if size(forceMatrix,1) ~= nDOF && size(forceMatrix,2) == nDOF
        forceMatrix = forceMatrix.';
    end
    if size(forceMatrix,1) ~= nDOF
        error('forceMatrix must have size nDOF x nStep.');
    end
    if size(forceMatrix,2) ~= numel(time)
        error('forceMatrix column count must match length(time).');
    end
end

function baseAcceleration = normalizeBaseAccelerationSize(baseAcceleration, nStep)
    if size(baseAcceleration,1) ~= nStep && size(baseAcceleration,2) == nStep
        baseAcceleration = baseAcceleration.';
    end
    if size(baseAcceleration,1) ~= nStep
        error('Base acceleration must have nStep rows.');
    end
    if size(baseAcceleration,2) > 3
        error('Base acceleration may have at most three columns: X, Y, Z.');
    end
    if size(baseAcceleration,2) < 3
        baseAcceleration(:,end+1:3) = 0;
    end
end

function forceMatrix = baseAccelerationToForceMatrix(M, eq_direction, baseAcceleration)
    % For base acceleration input, the equivalent inertial force is
    % f(t) = -M * R * a_g(t), where R is the influence matrix.
    R = eq_direction(:,1:3);
    forceMatrix = -M * R * baseAcceleration(:,1:3).';
end

function [U, V, A] = newmarkBetaLinear(M, C, K, F, dt, opt, waitbarTitle)
    if nargin < 7 || isempty(waitbarTitle)
        waitbarTitle = 'Newmark-beta dynamic analysis';
    end

    beta = opt.beta;
    gamma = opt.gamma;

    nDOF = size(M,1);
    nStep = size(F,2);

    U = zeros(nDOF, nStep);
    V = zeros(nDOF, nStep);
    A = zeros(nDOF, nStep);

    % Initial acceleration from equilibrium.
    A(:,1) = M \ (F(:,1) - C*V(:,1) - K*U(:,1));

    a0 = 1/(beta*dt^2);
    a1 = gamma/(beta*dt);
    a2 = 1/(beta*dt);
    a3 = 1/(2*beta) - 1;
    a4 = gamma/beta - 1;
    a5 = dt*(gamma/(2*beta) - 1);

    K_eff = K + a0*M + a1*C;

    % Progress gauge. The waitbar is updated at about 100 increments to avoid
    % slowing down the dynamic analysis while still showing progress clearly.
    hWait = waitbar(0, sprintf('%s: 0.0%%', waitbarTitle));
    cleanupObj = onCleanup(@() closeWaitbarIfValid(hWait)); %#ok<NASGU>
    updateInterval = max(1, floor((nStep-1)/100));

    for iStep = 1:nStep-1
        F_eff = F(:,iStep+1) ...
            + M*(a0*U(:,iStep) + a2*V(:,iStep) + a3*A(:,iStep)) ...
            + C*(a1*U(:,iStep) + a4*V(:,iStep) + a5*A(:,iStep));

        U(:,iStep+1) = K_eff \ F_eff;
        A(:,iStep+1) = a0*(U(:,iStep+1) - U(:,iStep)) - a2*V(:,iStep) - a3*A(:,iStep);
        V(:,iStep+1) = V(:,iStep) + dt*((1-gamma)*A(:,iStep) + gamma*A(:,iStep+1));

        if mod(iStep, updateInterval) == 0 || iStep == nStep-1
            progress = iStep/(nStep-1);
            waitbar(progress, hWait, sprintf('%s: %.1f%%', waitbarTitle, 100*progress));
            drawnow limitrate;
        end
    end

    % Remove negligible imaginary components that may arise from complex-valued
    % POD bases or numerical round-off. Time-domain outputs are real-valued.
    U = real(U);
    V = real(V);
    A = real(A);
end

function closeWaitbarIfValid(h)
    if ~isempty(h) && isvalid(h)
        close(h);
    end
end

function a_abs = addBaseAccelerationToRelativeResponse(a_rel, baseAcceleration, eq_direction)
    a_abs = real(a_rel);
    if isempty(baseAcceleration)
        return;
    end
    for iDir = 1:min(3,size(baseAcceleration,2))
        dofIndex = find(eq_direction(:,iDir) ~= 0);
        a_abs(dofIndex,:) = a_abs(dofIndex,:) + real(baseAcceleration(:,iDir)).';
    end
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

function plotTimeResponse(inputData, sensorResponse, fullOrderResult, outputIndex, plotTimeWindow)
    if isempty(outputIndex)
        warning('No output sensor response was selected for plotting.');
        return;
    end

    if nargin < 5
        plotTimeWindow = [];
    end

    iPlot = 1;
    t = inputData.time;

    figure('Name','Input base acceleration','Position',[100 100 800 320]);
    plot(t, inputData.baseAcceleration(:,1), 'DisplayName','X'); hold on;
    plot(t, inputData.baseAcceleration(:,2), 'DisplayName','Y');
    plot(t, inputData.baseAcceleration(:,3), 'DisplayName','Z');
    grid on; box on;
    xlabel('Time (s)');
    ylabel('Base acceleration');
    legend('Location','best');

    figure('Name','ROM sensor response','Position',[120 120 850 520]);
    subplot(3,1,1);
    plot(t, sensorResponse.u(iPlot,:)); grid on; box on;
    ylabel('Relative disp.');
    title(sprintf('Sensor node %d, DOF %d', sensorResponse.outputDOFs(iPlot,1), sensorResponse.outputDOFs(iPlot,2)));

    subplot(3,1,2);
    plot(t, sensorResponse.v(iPlot,:)); grid on; box on;
    ylabel('Relative vel.');

    subplot(3,1,3);
    plot(t, sensorResponse.a_absolute(iPlot,:)); grid on; box on;
    xlabel('Time (s)');
    ylabel('Absolute acc.');

    if ~isempty(fullOrderResult)
        if isempty(plotTimeWindow)
            idxPlot = true(size(t));
            windowTitle = 'full duration';
        else
            t_min = plotTimeWindow(1);
            t_max = plotTimeWindow(2);
            idxPlot = (t >= t_min) & (t <= t_max);
            windowTitle = sprintf('%.2f--%.2f s', t_min, t_max);
        end

        figure('Name','ROM-FOM displacement comparison','Position',[140 140 850 420]);
        plot(t(idxPlot), sensorResponse.u(iPlot,idxPlot), 'DisplayName','ROM'); hold on;
        plot(t(idxPlot), fullOrderResult.u(outputIndex(iPlot),idxPlot), '--', 'DisplayName','FOM');
        grid on; box on;
        xlabel('Time (s)');
        ylabel('Relative displacement');
        title(sprintf('ROM-FOM displacement comparison (%s)', windowTitle));
        legend('Location','best');
    end
end

function modelSmall = stripLargeFields(model, keepLarge)
    modelSmall = model;
    if ~keepLarge
        modelSmall.M = [];
        modelSmall.C = [];
        modelSmall.K = [];
    end
end
