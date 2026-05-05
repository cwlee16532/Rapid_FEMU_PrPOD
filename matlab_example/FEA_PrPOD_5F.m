function [FRF_X_hat, FRF_Y_hat, freq_ROM, details] = FEA_PrPOD_5F(mu_trial, ROM, r, opt)
%% FEA_PrPOD_5F.m
% -------------------------------------------------------------------------
% ROM-based FRF analysis for the 5-story frame structure
% -------------------------------------------------------------------------
% This function assembles M(mu), C(mu), and K(mu), interpolates the aligned
% POD basis Phi_tilde(mu), performs Galerkin projection, and computes the
% reconstructed FRFs for X- and Y-direction unit inputs.
%
% Paper notation preserved in variable names:
%   mu_trial   : trial parameter vector
%   Phi_tilde  : interpolated congruence-aligned POD basis
%   M_r,C_r,K_r: reduced mass, damping, and stiffness matrices
%   U_hat      : reconstructed full-order frequency response
% -------------------------------------------------------------------------

    if nargin < 4
        opt = struct();
    end
    opt = setDefaultFEAOptions(opt);

    mu_trial = mu_trial(:).';
    mu_trial = enforceParameterBounds(mu_trial, ROM.p_grid);

    %% Interpolate aligned POD basis Phi_tilde(mu_trial)
    if isfield(ROM, 'useFastInterpolant') && ROM.useFastInterpolant
        Phi_interp_full = squeeze(ROM.POM_interpolant( ...
            mu_trial(1), mu_trial(2), mu_trial(3), mu_trial(4), mu_trial(5)));
    
        Phi_tilde = Phi_interp_full(:, 1:r);
    else
        Phi_tilde = interpolateAlignedBasis( ...
            ROM.Phi_aligned_database, ROM.p_grid, mu_trial, r, ...
            romOptions.podInterpolationMethod, romOptions.podExtrapolationMethod);
    end

    %% Assemble full-order system matrices M(mu), K(mu)
    nodedata = xlsread(opt.modelFile, opt.nodeSheet);
    elementdata = xlsread(opt.modelFile, opt.elementSheet);
    elementdata = applyStoryParameters(elementdata, mu_trial, opt.storyBeamRows, opt.sectionColumns);
    [M, K, free_dof, dofs_5fs] = assembleFiveStoryFrame(nodedata, elementdata); %#ok<ASGLU>

    %% Damping matrix C(mu)
    switch lower(opt.dampingModel)
        case 'rayleigh'
            [phi, ev] = eigs(K, M, opt.nModesForDamping, opt.modalShift);
            fn = sqrt(diag(ev)) / (2*pi);
            rayleighCoeff = (0.5 * [1/(fn(1)*2*pi), fn(1)*2*pi; ...
                                    1/(fn(2)*2*pi), fn(2)*2*pi]) \ [opt.xi; opt.xi];
            C = rayleighCoeff(1)*M + rayleighCoeff(2)*K;

        case 'modal'
            [phi, ev] = eigs(K, M, opt.nModesForDamping, opt.modalShift);
            fn = sqrt(diag(ev)) / (2*pi); %#ok<NASGU>
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
            rayleighCoeff = [NaN; NaN];

        otherwise
            error('Unknown dampingModel: %s', opt.dampingModel);
    end

    %% Load vectors and Galerkin projection
    eq_direction = buildDirectionMatrix(size(M,1), dofs_5fs);
    F_X = eq_direction(:,1);
    F_Y = eq_direction(:,2);

    M_r = Phi_tilde' * M * Phi_tilde;
    C_r = Phi_tilde' * C * Phi_tilde;
    K_r = Phi_tilde' * K * Phi_tilde;
    F_r_X = Phi_tilde' * F_X;
    F_r_Y = Phi_tilde' * F_Y;

    %% Reduced FRF solution and full-space reconstruction
    freq_ROM = opt.frequencyVectorHz(:);
    nFreq = numel(freq_ROM);
    nDOF = size(M,1);
    FRF_X_hat = zeros(nDOF, nFreq);
    FRF_Y_hat = zeros(nDOF, nFreq);

    for iw = 1:nFreq
        omega = 2*pi*freq_ROM(iw);
        A_r = K_r + 1i*omega*C_r - omega^2*M_r;
        U_r_X = A_r \ F_r_X;
        U_r_Y = A_r \ F_r_Y;
        FRF_X_hat(:,iw) = Phi_tilde * U_r_X;
        FRF_Y_hat(:,iw) = Phi_tilde * U_r_Y;
    end

    %% Optional full-order FRF check retained from the original code
    if opt.computeFullOrderFRF
        FRF_X_full = zeros(nDOF, nFreq);
        FRF_Y_full = zeros(nDOF, nFreq);
        for iw = 1:nFreq
            omega = 2*pi*freq_ROM(iw);
            A = K + 1i*omega*C - omega^2*M;
            FRF_X_full(:,iw) = A \ F_X;
            FRF_Y_full(:,iw) = A \ F_Y;
        end
    else
        FRF_X_full = [];
        FRF_Y_full = [];
    end

    if nargout > 3
        details = struct();
        details.mu_trial = mu_trial;
        details.Phi_tilde = Phi_tilde;
        details.M_r = M_r;
        details.C_r = C_r;
        details.K_r = K_r;
        details.M = M;
        details.C = C;
        details.K = K;
        details.rayleighCoeff = rayleighCoeff;
        details.eq_direction = eq_direction;
        details.FRF_X_full = FRF_X_full;
        details.FRF_Y_full = FRF_Y_full;
    end
end

function opt = setDefaultFEAOptions(opt)
    defaults = struct();
    defaults.modelFile = 'model_5f';
    defaults.nodeSheet = 'Sheet2';
    defaults.elementSheet = 'Sheet1';
    defaults.storyBeamRows = {9:10; 25:26; 41:42; 57:58; 73:74};
    defaults.sectionColumns = 5:6;
    defaults.frequencyVectorHz = (0.05:0.05:30).';
    defaults.xi = 0.01;
    defaults.dampingModel = 'rayleigh';
    defaults.modalDampingRatios = 0.01*ones(10,1);
    defaults.nModesForDamping = 5;
    defaults.modalShift = 0.01;
    defaults.podInterpolationMethod = 'linear';
    defaults.podExtrapolationMethod = 'nearest';
    defaults.computeFullOrderFRF = false;

    names = fieldnames(defaults);
    for i = 1:numel(names)
        if ~isfield(opt, names{i}) || isempty(opt.(names{i}))
            opt.(names{i}) = defaults.(names{i});
        end
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

    % Interpolate the full aligned basis first, preserving all nPOM modes.
    Phi_reshaped = reshape(PhiDB, [n1, n2, n3, n4, n5, nDOF, nPOM]);

    Phi_interp_full = interpn(p_grid{1}, p_grid{2}, p_grid{3}, p_grid{4}, p_grid{5}, ...
        Phi_reshaped, ...
        mu_trial(1), mu_trial(2), mu_trial(3), mu_trial(4), mu_trial(5), ...
        interpolationMethod);

    if any(isnan(Phi_interp_full(:)))
        Phi_interp_full = interpn(p_grid{1}, p_grid{2}, p_grid{3}, p_grid{4}, p_grid{5}, ...
            Phi_reshaped, ...
            mu_trial(1), mu_trial(2), mu_trial(3), mu_trial(4), mu_trial(5), ...
            extrapolationMethod);
    end

    Phi_interp_full = squeeze(Phi_interp_full);

    if size(Phi_interp_full,1) ~= nDOF
        Phi_interp_full = reshape(Phi_interp_full, [nDOF, nPOM]);
    end

    % Extract only the first r modes after interpolation.
    Phi_tilde = Phi_interp_full(:, 1:r);

end

function mu = enforceParameterBounds(mu, p_grid)
    lb = [min(p_grid{1}(:)), min(p_grid{2}(:)), min(p_grid{3}(:)), min(p_grid{4}(:)), min(p_grid{5}(:))];
    ub = [max(p_grid{1}(:)), max(p_grid{2}(:)), max(p_grid{3}(:)), max(p_grid{4}(:)), max(p_grid{5}(:))];
    mu = max(min(mu, ub), lb);
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
