%% main_offline_build_ROM.m
% -------------------------------------------------------------------------
% Offline phase: build parametric frequency POD-ROM for the 5-story frame
% -------------------------------------------------------------------------
% This script constructs the offline ROM database used in the PrPOD-ROM
% FEMU framework.
%
% Main procedure:
%   1) Define operating points for five story-wise beam-section parameters
%   2) Assemble full-order M and K matrices for each operating point
%   3) Compute Rayleigh damping from the first two natural frequencies
%   4) Compute analytic FRFs over the target frequency range
%   5) Extract frequency-domain POD bases using SVD
%   6) Align POD bases to a common reference basis by Procrustes rotation
%   7) Save the offline ROM database for online FEMU and dynamic analysis
%
% Required files/functions:
%   - 5f_model.xlsx
%       Sheet1: element data
%       Sheet2: node data
%   - StiffMass.m
%   - stiffassemble.m
%   - massassemble.m
%
% Author: Chanwoo Lee
% Repository: https://github.com/cwlee16532/PrPOD_ROM_FE_model_updating_for_online_dynamic_analysis
% -------------------------------------------------------------------------

clear; close all; clc;

%% -------------------------
% User-defined settings
% -------------------------
modelFile = 'model_5f';   % Excel model file without extension is allowed
nodeSheet = 'Sheet2';        % node data sheet
elementSheet = 'Sheet1';     % element data sheet

% Operating points for the five updating parameters.
% Each parameter controls the beam section dimensions of one story.
paramValues = [0.004, 0.006, 0.008]; % paramValues = [0.004, 0.008];
[p1, p2, p3, p4, p5] = ndgrid(paramValues, paramValues, paramValues, paramValues, paramValues);
param_space = [p1(:), p2(:), p3(:), p4(:), p5(:)];

% Element rows corresponding to story-wise beam section parameters.
% Columns 5 and 6 of elementdata are replaced by each parameter value.
storyBeamRows = {
    9:10;    % 1F beam section area/dimension parameter
    25:26;   % 2F beam section area/dimension parameter
    41:42;   % 3F beam section area/dimension parameter
    57:58;   % 4F beam section area/dimension parameter
    73:74    % 5F beam section area/dimension parameter
};
sectionColumns = 5:6;

% Modal and damping settings
nModesForModalAnalysis = 50;
modalShift = 0.01;
xi = 0.01;                  % target damping ratio for Rayleigh damping

% Frequency grid for FRF snapshots
wVec = (0.05:0.05:30).';    % Hz
% Alternative logarithmic grid used in verification studies:
% log_wVec = linspace(log10(3), log10(30), 200);
% wVec = 10.^log_wVec(:);

% FRF input options
computeFRF_X  = false;      % keep optional X-only input FRF calculation
computeFRF_Y  = false;      % keep optional Y-only input FRF calculation
computeFRF_XY = true;       % default database uses combined X+Y input
includeDampingInFRF = true; % false: undamped FRF

% POD settings
maxPODRank = 300;            % [] keeps all available singular vectors
referenceOperatingPoint = 1;

% Optional verification block
runVerification = true;
param_test = [0.0075, 0.0045, 0.0075, 0.0045, 0.0075];
verificationRank = 200;
verificationOutDOF = 811;

% Output
outputFile = 'PODROM_5F_param0.004_0.008_POM300_inputXY.mat';
% load PODROM_5F_param0.004_0.008_POM300_inputXY.mat

% Preallocation
% -------------------------
nOp = size(param_space, 1);
U_X  = cell(nOp, 1);
U_Y  = cell(nOp, 1);
U_XY = cell(nOp, 1);
S_X  = cell(nOp, 1);
S_Y  = cell(nOp, 1);
S_XY = cell(nOp, 1);

U_X_final  = [];
U_Y_final  = [];
U_XY_final = [];

modelInfo = struct();
modelInfo.modelFile = modelFile;
modelInfo.nodeSheet = nodeSheet;
modelInfo.elementSheet = elementSheet;
modelInfo.paramValues = paramValues;
modelInfo.param_space = param_space;
modelInfo.storyBeamRows = storyBeamRows;
modelInfo.sectionColumns = sectionColumns;
modelInfo.wVec = wVec;
modelInfo.xi = xi;
modelInfo.includeDampingInFRF = includeDampingInFRF;
modelInfo.referenceOperatingPoint = referenceOperatingPoint;

%% =========================================================
% STEP 1) Generate POD bases at all operating points
% =========================================================
h = waitbar(0, 'Generating frequency POD-ROM database...');
cleanupObj = onCleanup(@() closeWaitbarIfValid(h));

for p = 1:nOp
    fprintf('Operating point %d / %d\n', p, nOp);

    %% Load and assemble full-order system
    nodedata = xlsread(modelFile, nodeSheet);
    elementdata = xlsread(modelFile, elementSheet);

    elementdata = applyStoryParameters(elementdata, param_space(p,:), storyBeamRows, sectionColumns);

    [M, K, free_dof, dofs_5fs] = assembleFiveStoryFrame(nodedata, elementdata);

    %% Modal analysis and Rayleigh damping
    [phi, ev] = eigs(K, M, nModesForModalAnalysis, modalShift);
    fn = sqrt(diag(ev)) / (2*pi);

    rayleighCoeff = (0.5 * [1/(fn(1)*2*pi), fn(1)*2*pi; ...
                            1/(fn(2)*2*pi), fn(2)*2*pi]) \ [xi; xi];
    D_ray = rayleighCoeff(1)*M + rayleighCoeff(2)*K;

    Xi = diag(phi' * D_ray * phi) ./ diag(2*sqrt((phi' * M * phi) .* (phi' * K * phi)));

    %% Ground-motion/load direction vectors
    eq_direction = buildDirectionMatrix(size(M,1), dofs_5fs);
    GF_X = eq_direction(:,1);
    GF_Y = eq_direction(:,2);

    %% Compute analytic FRFs
    if computeFRF_X
        FRF_X = computeAnalyticFRF(M, D_ray, K, GF_X, wVec, includeDampingInFRF);
        [U_X{p}, S_X{p}] = computePODBases(FRF_X, maxPODRank); 
    end

    if computeFRF_Y
        FRF_Y = computeAnalyticFRF(M, D_ray, K, GF_Y, wVec, includeDampingInFRF);
        [U_Y{p}, S_Y{p}] = computePODBases(FRF_Y, maxPODRank); 
    end

    if computeFRF_XY
        FRF_XY = computeAnalyticFRF(M, D_ray, K, GF_X + GF_Y, wVec, includeDampingInFRF);
        [U_XY{p}, S_XY{p}] = computePODBases(FRF_XY, maxPODRank); 
    end

    if p == 1
        modelInfo.nFreeDOF = size(M,1);
        modelInfo.free_dof = free_dof;
        modelInfo.dofs_5fs = dofs_5fs;
        modelInfo.rayleighCoeff_ref = rayleighCoeff;
        modelInfo.fn_ref = fn;
    end

    waitbar(p/nOp, h);
end

%% =========================================================
% STEP 2) Congruence / Procrustes basis alignment
% =========================================================
if computeFRF_X
    U_X_final = alignPODBasesToReference(U_X, referenceOperatingPoint);
end

if computeFRF_Y
    U_Y_final = alignPODBasesToReference(U_Y, referenceOperatingPoint);
end

if computeFRF_XY
    U_XY_final = alignPODBasesToReference(U_XY, referenceOperatingPoint);
end

%% =========================================================
% STEP 3) Optional POD-ROM adaptation verification
% =========================================================
if runVerification
    if ~computeFRF_XY || isempty(U_XY_final)
        error('Verification currently uses U_XY_final. Set computeFRF_XY = true.');
    end

    verificationResult = verifyInterpolatedROM( ...
        modelFile, nodeSheet, elementSheet, param_test, storyBeamRows, sectionColumns, ...
        p1, p2, p3, p4, p5, U_XY_final, wVec, xi, nModesForModalAnalysis, ...
        modalShift, verificationRank, verificationOutDOF, includeDampingInFRF); 
end

%% -------------------------
% Save offline database
% -------------------------
save(outputFile, ...
    'U_X', 'U_Y', 'U_XY', ...
    'S_X', 'S_Y', 'S_XY', ...
    'U_X_final', 'U_Y_final', 'U_XY_final', ...
    'param_space', 'paramValues', 'p1', 'p2', 'p3', 'p4', 'p5', ...
    'wVec', 'modelInfo', '-v7.3');

fprintf('Offline ROM database saved: %s\n', outputFile);

%% -------------------------
% Save fast interpolation object for online phase
% -------------------------
% This file is used to accelerate repeated basis interpolation:
%   Phi_tilde(mu_trial) = POM_G_5F(mu_1, mu_2, mu_3, mu_4, mu_5)
%
% Note:
%   The full number of stored POMs is preserved.
%   The online code should extract the first r modes after interpolation.

if computeFRF_XY
    U_reshaped = reshape(U_XY_final, ...
        [size(p1), size(U_XY_final,2), size(U_XY_final,3)]);

    POM_G_5F = griddedInterpolant( ...
        p1, p2, p3, p4, p5, ...
        U_reshaped, ...
        'linear', 'nearest');

    interpolationFile = 'POM_interpolant_5F_param0.004_0.008_POM300_inputXY.mat';

    save(interpolationFile, ...
        'POM_G_5F', ...
        'paramValues', 'param_space', ...
        'p1', 'p2', 'p3', 'p4', 'p5', ...
        'wVec', 'modelInfo', ...
        '-v7.3');

    fprintf('Fast POM interpolation object saved: %s\n', interpolationFile);
end
%% =========================================================
% Local helper functions
% =========================================================
function elementdata = applyStoryParameters(elementdata, param, storyBeamRows, sectionColumns)
    for k = 1:numel(storyBeamRows)
        elementdata(storyBeamRows{k}, sectionColumns) = param(k);
    end
end

function [M, K, free_dof, dofs_5fs] = assembleFiveStoryFrame(nodedata, elementdata)
    N  = size(nodedata, 1);
    El = size(elementdata, 1);

    Fixity = nodedata(:,5:10);
    bounded_dof = find(Fixity' ~= 0);
    free_dof = setdiff(1:6*N, bounded_dof);

    dofs_5fs = zeros(N*6, 2);
    for i = 1:N
        dofs_5fs(6*i-5:6*i,1) = i;
        dofs_5fs(6*i-5:6*i,2) = (1:6).';
    end
    dofs_5fs(bounded_dof,:) = [];

    KGlobal = zeros(6*N, 6*N);
    MGlobal = zeros(6*N, 6*N);

    for e = 1:El
        elementdata_i = elementdata(e,:);

        node1 = elementdata_i(2);
        node2 = elementdata_i(3);
        b = elementdata_i(5);
        d = elementdata_i(6);
        E = elementdata_i(7);
        G = elementdata_i(8);
        Iz = b*d^3/12;
        Iy = d*b^3/12;
        J = Iz + Iy;
        rho = elementdata_i(17);

        XYZ1 = nodedata(node1,2:4);
        XYZ2 = nodedata(node2,2:4);

        x1 = XYZ1(1); y1 = XYZ1(2); z1 = XYZ1(3);
        x2 = XYZ2(1); y2 = XYZ2(2); z2 = XYZ2(3);

        [stiff_local, mass_local] = StiffMass(E, G, d, b, Iy, Iz, J, rho, x1, y1, z1, x2, y2, z2);

        KGlobal = stiffassemble(KGlobal, stiff_local, node1, node2);
        MGlobal = massassemble(MGlobal, mass_local, node1, node2);
    end

    K = KGlobal(free_dof, free_dof);
    M = MGlobal(free_dof, free_dof);
end

function eq_direction = buildDirectionMatrix(nDOF, dofs_5fs)
    eq_direction = zeros(nDOF, 6);
    for i = 1:6
        row = find(dofs_5fs(:,2) == i);
        eq_direction(row, i) = 1;
    end
end

function FRF = computeAnalyticFRF(M, C, K, forceVector, wVec, includeDamping)
    nW = numel(wVec);
    nDOF = size(M, 1);
    FRF = complex(zeros(nDOF, nW));

    for iw = 1:nW
        w = 2*pi*wVec(iw);
        if includeDamping
            dynamicStiffness = K + 1i*w*C - w^2*M;
        else
            dynamicStiffness = K - w^2*M;
        end
        FRF(:,iw) = dynamicStiffness \ forceVector;
    end
end

function [U, S] = computePODBases(snapshotMatrix, maxPODRank)
    if isempty(maxPODRank)
        r = min(size(snapshotMatrix));
    else
        r = min(maxPODRank, min(size(snapshotMatrix)));
    end

    % svds is efficient for large sparse/rectangular matrices. If it fails,
    % fall back to economy-size svd for portability.
    try
        [U, S, ~] = svds(snapshotMatrix, r, 'largest');
    catch
        [Ufull, Sfull, ~] = svd(snapshotMatrix, 'econ');
        U = Ufull(:,1:r);
        S = Sfull(1:r,1:r);
    end
end

function U_final = alignPODBasesToReference(U_cell, refIndex)
    nOp = numel(U_cell);
    refBasis = U_cell{refIndex};
    [nDOF, r] = size(refBasis);
    U_final = complex(zeros(nOp, nDOF, r));

    for p = 1:nOp
        currentBasis = U_cell{p};
        crossBasis = currentBasis' * refBasis;
        [U, ~, V] = svd(crossBasis, 'econ');
        POM_trans = U * V';
        U_final(p,:,:) = currentBasis * POM_trans;
    end
end

function verificationResult = verifyInterpolatedROM( ...
    modelFile, nodeSheet, elementSheet, param_test, storyBeamRows, sectionColumns, ...
    p1, p2, p3, p4, p5, U_XY_final, wVec, xi, nModesForModalAnalysis, ...
    modalShift, rDOF, outDOF, includeDampingInFRF)

    fprintf('Running interpolated ROM verification...\n');

    U_reshaped = reshape(U_XY_final, [size(p1), size(U_XY_final,2), size(U_XY_final,3)]);
    POM_interp = squeeze(interpn(p1, p2, p3, p4, p5, U_reshaped, ...
        param_test(1), param_test(2), param_test(3), param_test(4), param_test(5), 'linear'));

    POM_interp = POM_interp(:,1:rDOF);

    nodedata = xlsread(modelFile, nodeSheet);
    elementdata = xlsread(modelFile, elementSheet);
    elementdata = applyStoryParameters(elementdata, param_test, storyBeamRows, sectionColumns);
    [M, K, ~, dofs_5fs] = assembleFiveStoryFrame(nodedata, elementdata);

    [phi, ev] = eigs(K, M, nModesForModalAnalysis, modalShift);
    fn = sqrt(diag(ev))/(2*pi);
    rayleighCoeff = (0.5 * [1/(fn(1)*2*pi), fn(1)*2*pi; ...
                            1/(fn(2)*2*pi), fn(2)*2*pi]) \ [xi; xi];
    D_ray = rayleighCoeff(1)*M + rayleighCoeff(2)*K;

    eq_direction = buildDirectionMatrix(size(M,1), dofs_5fs);
    GF_X = eq_direction(:,1);
    GF_Y = eq_direction(:,2);
    GF_XY = GF_X + GF_Y;

    FRF_FOM = computeAnalyticFRF(M, D_ray, K, GF_XY, wVec, includeDampingInFRF);

    rM = POM_interp' * M * POM_interp;
    rD = POM_interp' * D_ray * POM_interp;
    rK = POM_interp' * K * POM_interp;
    rF = POM_interp' * GF_XY;

    rFRF = complex(zeros(rDOF, numel(wVec)));
    for iw = 1:numel(wVec)
        w = 2*pi*wVec(iw);
        if includeDampingInFRF
            rFRF(:,iw) = (rK + 1i*w*rD - w^2*rM) \ rF;
        else
            rFRF(:,iw) = (rK - w^2*rM) \ rF;
        end
    end
    FRF_ROM = POM_interp * rFRF;

    verificationResult.param_test = param_test;
    verificationResult.outDOF = outDOF;
    verificationResult.magFOM = abs(FRF_FOM(outDOF,:));
    verificationResult.magROM = abs(FRF_ROM(outDOF,:));
    verificationResult.phaseFOM = angle(FRF_FOM(outDOF,:));
    verificationResult.phaseROM = angle(FRF_ROM(outDOF,:));

    figure('Name', 'Interpolated POD-ROM verification', 'Position', [100, 100, 800, 300]);
    subplot(2,1,1)
    semilogy(wVec, verificationResult.magFOM); hold on;
    semilogy(wVec, verificationResult.magROM, '--');
    xlabel('Frequency (Hz)'); ylabel('Magnitude');
    legend('FOM', 'POD-ROM'); grid on; box on;

    subplot(2,1,2)
    plot(wVec, verificationResult.phaseFOM); hold on;
    plot(wVec, verificationResult.phaseROM, '--');
    xlabel('Frequency (Hz)'); ylabel('Phase');
    legend('FOM', 'POD-ROM'); grid on; box on;
end

function closeWaitbarIfValid(h)
    if ~isempty(h) && isvalid(h)
        close(h);
    end
end
