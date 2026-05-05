%% main_online_FEMU_PSO.m
% -------------------------------------------------------------------------
% Online phase: PrPOD-ROM-based FEMU using Particle Swarm Optimization (PSO)
% -------------------------------------------------------------------------
% This script performs the online model-updating phase for the numerical
% validation example of the 5-story frame structure.
%
% Main procedure:
%   1) Generate virtual sensor FRF data using Simul_expt_5f.m
%   2) Load the offline parametric POD-ROM database
%   3) Define measured FRF data H_meas
%   4) Run progressive PSO-based FEMU with increasing ROM rank r
%   5) Save updated parameters, objective values, and PSO history
%
% Paper notation used in variable names:
%   mu_trial       : trial updating parameter vector
%   mu_star        : optimal updating parameter vector
%   r, M_final     : current and final retained POD ranks
%   J_mu           : FRF-based objective function
%   epsilon_1      : stagnation threshold of best objective value
%   epsilon_2_s    : CV-based particle convergence threshold at stage s
%   epsilon_3_s    : swarm cost uniformity threshold at stage s
%   T_s            : stagnation window size
%   N_p            : swarm size
%   chi            : inertia weight range
%
% Required files/functions:
%   - Simul_expt_5f.m
%   - FEA_PrPOD_5F.m
%   - costFunction_PrPOD_ROM_5F.m
%   - PODROM_5F_param0.004_0.008_POM600_inputXY.mat
%   - model_5f.xlsx
%   - StiffMass.m, stiffassemble.m, massassemble.m
%   - Optimization Toolbox: particleswarm
% -------------------------------------------------------------------------

clear; close all; clc;

%% -------------------------
% User-defined paths
% -------------------------
% offlineROMFile = 'PODROM_5F_param0.004_0.008_POM600_inputXY.mat';
offlineROMFile = 'POM_interpolant_5F_param0.004_0.008_POM300_inputXY.mat';
resultFile     = 'FEMU_PSO_5F_PrPOD_result.mat';

modelFile    = 'model_5f';
nodeSheet    = 'Sheet2';
elementSheet = 'Sheet1';

%% -------------------------
% Virtual measurement settings
% -------------------------
% method = 'analytic' directly computes analytic FRFs from M(mu), C(mu), K(mu).
% method = 'dynamic' preserves the commented research-code workflow:
% white-noise base excitation -> Newmark-beta -> modalfrf FRF estimation.
exptOptions = struct();
exptOptions.method = 'analytic';
exptOptions.modelFile = modelFile;
exptOptions.nodeSheet = nodeSheet;
exptOptions.elementSheet = elementSheet;
exptOptions.mu_true = [0.004 0.005 0.006 0.007 0.008];
% exptOptions.mu_true = 0.0065*ones(1,5);   % optional closely spaced mode case
exptOptions.frequencyRangeHz = [3 30];
exptOptions.sensorNodes = [105, 113, 121, 129, 137, 108, 116, 124, 132, 140];
exptOptions.sensorDOFs = [1 2];             % translational X and Y DOFs
exptOptions.addNoise = false;
exptOptions.noiseLevelPercent = 1;
exptOptions.makePlots = true;

fprintf('Generating virtual sensor data...\n');
tic;
[H_meas, virtualData, modelRef] = Simul_expt_5f(exptOptions); %#ok<ASGLU>
fprintf('Virtual sensor data generated in %.3f s.\n', toc);

%% -------------------------
% Load offline parametric POD-ROM database
% -------------------------
fprintf('Loading offline POD-ROM database: %s\n', offlineROMFile);
tic;
ROMDB = load(offlineROMFile);
fprintf('Offline database loaded in %.3f s.\n', toc);

ROM = struct();

if isfield(ROMDB, 'POM_G_5F')
    ROM.POM_interpolant = ROMDB.POM_G_5F;
    ROM.useFastInterpolant = true;

elseif isfield(ROMDB, 'POM_G_5f')
    ROM.POM_interpolant = ROMDB.POM_G_5f;
    ROM.useFastInterpolant = true;

elseif isfield(ROMDB, 'U_XY_final')
    ROM.Phi_aligned_database = ROMDB.U_XY_final;
    ROM.useFastInterpolant = false;

elseif isfield(ROMDB, 'Phi_aligned_database')
    ROM.Phi_aligned_database = ROMDB.Phi_aligned_database;
    ROM.useFastInterpolant = false;

else
    error('No POD basis database or interpolation object was found.');
end

if all(isfield(ROMDB, {'p1','p2','p3','p4','p5'}))
    ROM.p_grid = {ROMDB.p1, ROMDB.p2, ROMDB.p3, ROMDB.p4, ROMDB.p5};
else
    values = [0.004 0.006 0.008];
    [p1,p2,p3,p4,p5] = ndgrid(values, values, values, values, values);
    ROM.p_grid = {p1,p2,p3,p4,p5};
end

if isfield(ROMDB, 'param_space')
    ROM.mu_operating_points = ROMDB.param_space;
end

% Optional cumulative-energy check retained from the research script.
runEnergyCheck = true;
energyCheckRank = 300;
if runEnergyCheck && isfield(ROMDB, 'S_XY') && ~isempty(ROMDB.S_XY{1})
    singularValues = diag(ROMDB.S_XY{1});
    retainedEnergy = sum(singularValues(1:energyCheckRank)) / sum(singularValues);
    fprintf('Retained POD energy at r = %d: %.12f\n', energyCheckRank, retainedEnergy);
end

%% -------------------------
% FEMU / PSO settings
% -------------------------
nvars = 5;
mu_lower_bound = 0.004 * ones(1,nvars);
mu_upper_bound = 0.008 * ones(1,nvars);

% Progressive mode counts r_1 < r_2 < ... < M_final.
% progressiveRanks = 300;              % original one-stage example
progressiveRanks = [100 200 300];  % optional multi-stage PrPOD run
M_final = progressiveRanks(end);
N_stage = numel(progressiveRanks);

% PSO hyperparameters
N_p = 50;
c_1 = 1.49;
c_2 = 1.49;
maxIterations = 1000;
maxStallIterations = 20;
functionTolerance = 1e-6;

% Multi-criteria convergence thresholds from the paper notation.
epsilon_1 = 1.0e-2;       % best-cost stagnation threshold
epsilon_2_start = 0.10;   % CV threshold, early stage
epsilon_2_final = 0.03;   % CV threshold, final stage
epsilon_3_start = 1.00;   % cost-uniformity threshold, early stage
epsilon_3_final = 0.10;   % cost-uniformity threshold, final stage
T_s = 10;

% Adaptive inertia-weight upper bound: chi ~ U(0, chi_upper_s)
chi_upper_start = 1.10;
chi_upper_final = 0.10;

% Optional exact-solution particles retained from the research script.
useExactParticles = false;
mu_exact = exptOptions.mu_true;
exact_particles = repmat(mu_exact, 5, 1); %#ok<NASGU>

%% -------------------------
% Objective/ROM analysis settings
% -------------------------
objectiveOptions = struct();
objectiveOptions.objectiveMetric = 'log_deviation';
objectiveOptions.useMagnitudeOnly = true;
objectiveOptions.frfInterpolationMethod = 'linear';
objectiveOptions.makeComparisonPlots = false;

romOptions = struct();
romOptions.modelFile = modelFile;
romOptions.nodeSheet = nodeSheet;
romOptions.elementSheet = elementSheet;
romOptions.frequencyVectorHz = (0.05:0.05:30).';
romOptions.xi = 0.01;
romOptions.dampingModel = 'rayleigh';
romOptions.modalDampingRatios = 0.01*ones(10,1);
romOptions.nModesForDamping = 5;
romOptions.modalShift = 0.01;
romOptions.podInterpolationMethod = 'linear';
romOptions.podExtrapolationMethod = 'nearest';
romOptions.computeFullOrderFRF = false;  % optional FOM-vs-ROM check

objectiveOptions.romOptions = romOptions;

% Diagnostic objective evaluation
runCostFunctionTest = true;

if runCostFunctionTest
    fprintf('\nRunning diagnostic cost-function evaluations...\n');

    tic;
    J_exact = costFunction_PrPOD_ROM_5F(H_meas, mu_exact, ROM, progressiveRanks(1), objectiveOptions);
    time_J_exact = toc;

    tic;
    J_mid = costFunction_PrPOD_ROM_5F(H_meas, 0.006*ones(1,5), ROM, progressiveRanks(1), objectiveOptions);
    time_J_mid = toc;

    fprintf('Diagnostic J(mu_exact) = %.6e | elapsed time = %.4f s\n', ...
        J_exact, time_J_exact);

    fprintf('Diagnostic J(mu_mid)   = %.6e | elapsed time = %.4f s\n', ...
        J_mid, time_J_mid);
end

%% =========================================================
% STEP 1) Progressive PSO-based FEMU
% =========================================================
PSO_data = struct();
x_data = cell(N_stage,1);
fval_data = cell(N_stage,1);
exitflag_data = cell(N_stage,1);
output_data = cell(N_stage,1);
points_data = cell(N_stage,1);
elapsed_time = zeros(N_stage,1);

for s = 1:N_stage
    r = progressiveRanks(s);

    if N_stage == 1
        epsilon_2_s = epsilon_2_final;
        epsilon_3_s = epsilon_3_final;
        chi_upper_s = chi_upper_start;
    else
        stageRatio = (s-1) / (N_stage-1);
        epsilon_2_s = epsilon_2_start - stageRatio*(epsilon_2_start - epsilon_2_final);
        epsilon_3_s = epsilon_3_start - stageRatio*(epsilon_3_start - epsilon_3_final);
        chi_upper_s = chi_upper_start - stageRatio*(chi_upper_start - chi_upper_final);
    end

    fprintf('\n=========================================================\n');
    fprintf('PrPOD stage %d / %d: r = %d, M_final = %d\n', s, N_stage, r, M_final);
    fprintf('epsilon_2^(s) = %.4g, epsilon_3^(s) = %.4g, chi in [0, %.4g]\n', ...
        epsilon_2_s, epsilon_3_s, chi_upper_s);
    fprintf('=========================================================\n');

    outputFunction = @(optimValues,state) psoOutputFunction_PrPOD( ...
        optimValues, state, s, N_stage, epsilon_1, epsilon_2_s, epsilon_3_s, T_s);

    fun = @(mu_trial) costFunction_PrPOD_ROM_5F(H_meas, mu_trial, ROM, r, objectiveOptions);

    if s == 1
        options = optimoptions('particleswarm', ...
            'CreationFcn', 'pswcreationuniform', ...
            'FunctionTolerance', functionTolerance, ...
            'InertiaRange', [0 chi_upper_s], ...
            'MaxIterations', maxIterations, ...
            'MaxStallIterations', maxStallIterations, ...
            'SelfAdjustmentWeight', c_1, ...
            'SocialAdjustmentWeight', c_2, ...
            'SwarmSize', N_p, ...
            'Display', 'iter', ...
            'OutputFcn', outputFunction);
        if useExactParticles
            options = optimoptions(options, 'InitialPoints', exact_particles);
        end
    else
        options = optimoptions('particleswarm', ...
            'InitialPoints', points_data{s-1}.X, ...
            'FunctionTolerance', functionTolerance, ...
            'InertiaRange', [0 chi_upper_s], ...
            'MaxIterations', maxIterations, ...
            'MaxStallIterations', maxStallIterations, ...
            'SelfAdjustmentWeight', c_1, ...
            'SocialAdjustmentWeight', c_2, ...
            'SwarmSize', N_p, ...
            'Display', 'iter', ...
            'OutputFcn', outputFunction);
    end

    global PSO_HISTORY_PRPOD
    PSO_HISTORY_PRPOD = [];

    tic;
    [mu_star_s, J_star_s, exitflag_s, output_s, points_s] = ...
        particleswarm(fun, nvars, mu_lower_bound, mu_upper_bound, options);
    elapsed_time(s) = toc;

    x_data{s} = mu_star_s;
    fval_data{s} = J_star_s;
    exitflag_data{s} = exitflag_s;
    output_data{s} = output_s;
    points_data{s} = points_s;

    stageName = sprintf('pso_history_stage%d_rank%d', s, r);
    PSO_data.(stageName) = PSO_HISTORY_PRPOD;

    fprintf('Stage %d completed in %.3f s.\n', s, elapsed_time(s));
    fprintf('mu_star = %s\n', mat2str(mu_star_s, 8));
    fprintf('J(mu_star) = %.6e\n', J_star_s);
end

mu_star = x_data{end};
J_mu_star = fval_data{end};

%% =========================================================
% STEP 2) Optional visualization of PSO convergence
% =========================================================
makePSOPlots = true;
if makePSOPlots
    plotPSOHistory(PSO_data, elapsed_time, mu_lower_bound, mu_upper_bound);
end

%% -------------------------
% Save results
% -------------------------
save(resultFile, ...
    'mu_star', 'J_mu_star', 'x_data', 'fval_data', 'exitflag_data', ...
    'output_data', 'points_data', 'elapsed_time', 'PSO_data', ...
    'progressiveRanks', 'M_final', 'H_meas', 'exptOptions', ...
    'objectiveOptions', 'ROM', 'mu_exact', '-v7.3');

fprintf('\nOnline FEMU completed successfully.\n');
fprintf('Final mu_star: %s\n', mat2str(mu_star, 8));
fprintf('Final objective J(mu_star): %.6e\n', J_mu_star);
fprintf('Result file: %s\n', resultFile);

%% =========================================================
% Local helper functions
% =========================================================
function stop = psoOutputFunction_PrPOD(optimValues, state, stageIndex, nStage, epsilon_1, epsilon_2_s, epsilon_3_s, T_s)
% Stores PSO history and implements the three convergence criteria:
%   (i)   stagnation of best cost value, epsilon_1, over T_s iterations
%   (ii)  CV-based particle convergence, epsilon_2^(s)
%   (iii) uniformity of swarm performance, epsilon_3^(s)

    stop = false;
    global PSO_HISTORY_PRPOD

    if strcmp(state, 'init')
        PSO_HISTORY_PRPOD = struct( ...
            'iteration', [], ...
            'bestfval', [], ...
            'bestx', [], ...
            'meanfval', [], ...
            'score', [], ...
            'swarm', [], ...
            'std_swarm', [], ...
            'CV_mu', [], ...
            'epsilon_1', epsilon_1, ...
            'epsilon_2_s', epsilon_2_s, ...
            'epsilon_3_s', epsilon_3_s, ...
            'T_s', T_s);
        return;
    end

    if strcmp(state, 'iter')
        PSO_HISTORY_PRPOD.iteration(end+1) = optimValues.iteration;
        PSO_HISTORY_PRPOD.bestfval(end+1) = optimValues.bestfval;
        PSO_HISTORY_PRPOD.bestx(:,end+1) = optimValues.bestx(:);
        PSO_HISTORY_PRPOD.meanfval(end+1) = optimValues.meanfval;
        PSO_HISTORY_PRPOD.score{end+1} = optimValues.swarmfvals;
        PSO_HISTORY_PRPOD.swarm{end+1} = optimValues.swarm;
        PSO_HISTORY_PRPOD.std_swarm(end+1,:) = std(optimValues.swarm, 0, 1);
        PSO_HISTORY_PRPOD.CV_mu(end+1,:) = std(optimValues.swarm, 0, 1) ./ (abs(optimValues.bestx(:)).' + eps);

        J_best = optimValues.bestfval;
        J_mean = optimValues.meanfval;
        CV_mu_current = PSO_HISTORY_PRPOD.CV_mu(end,:);

        isCVConverged = all(CV_mu_current < epsilon_2_s);
        isUniform = abs(J_mean - J_best) / (abs(J_best) + eps) < epsilon_3_s;

        isStagnated = false;
        if stageIndex == nStage && numel(PSO_HISTORY_PRPOD.bestfval) >= T_s + 1
            J_now = PSO_HISTORY_PRPOD.bestfval(end);
            J_old = PSO_HISTORY_PRPOD.bestfval(end - T_s);
            relativeChange = abs(J_now - J_old) / (abs(J_old) + eps);
            isStagnated = relativeChange < epsilon_1;
        end

        if stageIndex < nStage
            if isCVConverged || isUniform
                stop = true;
                fprintf('Stopping stage %d: epsilon_2 or epsilon_3 criterion satisfied.\n', stageIndex);
            end
        else
            if isStagnated && (isCVConverged || isUniform)
                stop = true;
                fprintf('Stopping final stage: epsilon_1 and epsilon_2/epsilon_3 criteria satisfied.\n');
            end
        end
    end
end

function plotPSOHistory(PSO_data, elapsed_time, mu_lower_bound, mu_upper_bound)

    historyNames = fieldnames(PSO_data);
    if isempty(historyNames)
        return;
    end

    figure('Name','PSO parameter evolution','Position',[100 100 850 420]);
    hold on; grid on; box on;

    markerList = {'+','o','*','x','diamond'};

    % 🔥 핵심: j별 고정 색상
    colorList = lines(5);  % 5 parameters → distinct colors

    cumulativeTime = 0;

    for s = 1:numel(historyNames)
        hist_s = PSO_data.(historyNames{s});
        if ~isfield(hist_s, 'swarm') || isempty(hist_s.swarm)
            continue;
        end

        nIter = numel(hist_s.swarm);

        for iter = 1:nIter
            yTime = cumulativeTime + elapsed_time(s)*iter/nIter;
            swarm = hist_s.swarm{iter};

            for j = 1:min(5,size(swarm,2))

                if s == 1 && iter == 1
                    plot(swarm(:,j), yTime*ones(size(swarm,1),1), ...
                        'LineStyle','none', ...
                        'Marker',markerList{j}, ...
                        'Color',colorList(j,:), ...
                        'DisplayName',sprintf('\\mu_%d',j));
                else
                    plot(swarm(:,j), yTime*ones(size(swarm,1),1), ...
                        'LineStyle','none', ...
                        'Marker',markerList{j}, ...
                        'Color',colorList(j,:), ...
                        'HandleVisibility','off');
                end

            end
        end

        cumulativeTime = cumulativeTime + elapsed_time(s);
    end

    xlim([min(mu_lower_bound)-0.0002, max(mu_upper_bound)+0.0002]);
    xlabel('Updating parameter, \mu');
    ylabel('Cumulative computation time (s)');
    legend('Location','best');

    f_lines = min(mu_lower_bound):0.001:max(mu_upper_bound);
    yLimits = ylim;

    for k = 1:numel(f_lines)
        line([f_lines(k) f_lines(k)], yLimits, ...
            'Color', [0 0 0], ...
            'LineStyle', '--', ...
            'HandleVisibility','off');
    end
end
