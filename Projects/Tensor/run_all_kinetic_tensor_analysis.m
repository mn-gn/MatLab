%% run_all_kinetic_tensor_analysis.m
% Final entry point for the tensor/open-generator workflow.
%
% Default run:
%   1) sequential-like q-t mock data -> shifted-time tensor -> closure rotation
%   2) tensor-derived manifold -> sparse open first-order generator fit
%   3) Debye branch and parallel mock data -> open-generator validation
%
% The lower-level method scripts are kept in analysis_modules.

clear; clc; close all;

rootDir = fileparts(mfilename("fullpath"));
moduleDir = fullfile(rootDir, "analysis_modules");
resultRoot = fullfile(rootDir, "analysis_results");
matlabRoot = fullfile(rootDir, "..", "..");
sharedFunctionDir = fullfile(matlabRoot, "functions");
mockDataDir = fullfile(matlabRoot, "mock data");

if exist(sharedFunctionDir, "dir") == 7
    addpath(sharedFunctionDir);
end
addpath(moduleDir);

cfg = struct();
cfg.runSequentialTensor = true;
cfg.runSequentialOpenGenerator = true;
cfg.runDebyeBranchParallelValidation = true;

cfg.qFile = fullfile(mockDataDir, "q.csv");
cfg.tFile = fullfile(mockDataDir, "t_81points.csv");
cfg.dataFile = fullfile(mockDataDir, "q_delta_S_noise_10%_81points.csv");

cfg.sequentialTensorOutDir = fullfile(resultRoot, "01_sequential_tensor_modes");
cfg.sequentialGeneratorOutDir = fullfile(resultRoot, "02_sequential_sparse_open_generator");
cfg.branchParallelOutDir = fullfile(resultRoot, "03_debye_branch_parallel_validation");

cfg.numPC = 3;
cfg.numTau = 3;
cfg.numStart = 100;
cfg.cpRank = 3;

ensure_folder(resultRoot);
write_run_manifest(resultRoot, cfg);

fprintf("\n=== Tensor/open-generator integrated analysis ===\n");
fprintf("Root: %s\n", rootDir);
fprintf("Results: %s\n", resultRoot);

if cfg.runSequentialTensor
    fprintf("\n[1/3] Building shifted-time tensor and closure-rotated modes...\n");
    tensorOverride = struct();
    tensorOverride.qFile = cfg.qFile;
    tensorOverride.tFile = cfg.tFile;
    tensorOverride.dataFile = cfg.dataFile;
    tensorOverride.outDir = cfg.sequentialTensorOutDir;
    tensorOverride.numPC = cfg.numPC;
    tensorOverride.numTau = cfg.numTau;
    tensorOverride.numStart = cfg.numStart;
    tensorOverride.cpRank = cfg.cpRank;
    run_analysis_module(fullfile(moduleDir, "03_find_tensor_modes_with_closure_rotation.m"), tensorOverride);
    cd(rootDir);
end

if cfg.runSequentialOpenGenerator
    fprintf("\n[2/3] Fitting sparse open generator from tensor-derived manifold...\n");
    generatorOverride = struct();
    generatorOverride.resultFile = fullfile(cfg.sequentialTensorOutDir, "workspace_result.mat");
    generatorOverride.outDir = cfg.sequentialGeneratorOutDir;
    if exist(generatorOverride.resultFile, "file") ~= 2
        error("Sequential tensor result is missing. Run step 1 first or set generatorOverride.resultFile to an existing workspace_result.mat file.");
    end
    run_analysis_module(fullfile(moduleDir, "04_fit_sparse_open_generator_from_tensor.m"), generatorOverride);
    cd(rootDir);
end

if cfg.runDebyeBranchParallelValidation
    fprintf("\n[3/3] Validating branch and parallel Debye mock data...\n");
    debyeOverride = struct();
    debyeOverride.outDir = cfg.branchParallelOutDir;
    run_analysis_module(fullfile(moduleDir, "05_validate_debye_branch_parallel_open_generator.m"), debyeOverride);
    cd(rootDir);
end

fprintf("\nDone. Open this folder for all outputs:\n%s\n", resultRoot);

function run_analysis_module(modulePath, override)
    if exist(modulePath, "file") ~= 2
        error("Required module is missing: %s", modulePath);
    end
    cfg_override = override; %#ok<NASGU>
    run(modulePath);
end

function ensure_folder(folderPath)
    if exist(folderPath, "dir") ~= 7
        mkdir(folderPath);
    end
end

function write_run_manifest(resultRoot, cfg)
    ensure_folder(resultRoot);
    manifestFile = fullfile(resultRoot, "run_manifest.txt");
    fid = fopen(manifestFile, "w");
    if fid < 0
        warning("Could not write run manifest: %s", manifestFile);
        return;
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, "Tensor/open-generator integrated analysis\n");
    fprintf(fid, "Created: %s\n\n", datestr(now));
    fprintf(fid, "Sequential input q: %s\n", cfg.qFile);
    fprintf(fid, "Sequential input t: %s\n", cfg.tFile);
    fprintf(fid, "Sequential input data: %s\n\n", cfg.dataFile);
    fprintf(fid, "01 sequential tensor modes: %s\n", cfg.sequentialTensorOutDir);
    fprintf(fid, "02 sequential sparse open generator: %s\n", cfg.sequentialGeneratorOutDir);
    fprintf(fid, "03 Debye branch/parallel validation: %s\n", cfg.branchParallelOutDir);
end
