clear
clc
close all

%% Load the selected validation profile

% PFAL_MAIN_VALID_CONFIG can point to any of the file, web-test, or custom
% profiles. Keeping selection here makes the script convenient to run while
% leaving data addresses and scientific settings in JSON.
cfg_valid = UTILS.LOAD_MAIN_VALID_CONFIG;

%% Validate LD2 results against TEM and tandem measurements

% The returned tables remain in the workspace for immediate inspection.
% Stable manuscript figures and timestamped audit files are written according
% to the outputs section of the selected configuration.
validation_result = UTILS.RUN_MAIN_VALIDATION(cfg_valid);
validation_metrics = validation_result.metrics;
validation_predictions = validation_result.predictions;
experimental_relation_summary = ...
    validation_result.experimental_relation_summary;
