function [data, dataset_src] = LOAD_SCATTER_LD1_DATASET(dataset_cfg)
%LOAD_SCATTER_LD1_DATASET Load the LD1 source library used by main_scatter.
%   DATA = UTILS.LOAD_SCATTER_LD1_DATASET(DATASET_CFG) loads the MAT-file
%   variable declared in the dataset section of the main_scatter config.
%
%   [DATA, DATASET_SRC] = UTILS.LOAD_SCATTER_LD1_DATASET(...) also returns
%   the normalized dataset metadata that was used during the load.

% The dataset section is required because the script now receives its LD1
% source library information from JSON rather than hard-coded paths.
if nargin < 1 || ~isstruct(dataset_cfg)
    error('PFAL:LOAD_SCATTER_LD1_DATASET:MissingDatasetConfig', ...
        'The dataset configuration must be provided as a struct.');
end

% Validate the three fields needed to identify and load the MAT payload.
dataset_src.id = require_text_field(dataset_cfg, 'id', 'dataset identifier');
dataset_src.file = require_text_field(dataset_cfg, 'file', 'dataset file path');
dataset_src.variable = require_text_field(dataset_cfg, 'variable', 'dataset variable name');

% The config loader attaches config_dir so relative dataset paths can be
% resolved from the same folder as the JSON file.
if isfield(dataset_cfg, 'config_dir') && ~isempty(dataset_cfg.config_dir)
    dataset_src.config_dir = char(dataset_cfg.config_dir);
else
    dataset_src.config_dir = fullfile(repo_root(), 'config');
end

% Allow either absolute paths or paths written relative to the config file.
dataset_src.file = resolve_dataset_path(dataset_src.file, dataset_src.config_dir);

% Confirm the target MAT file exists before attempting to import data.
if ~isfile(dataset_src.file)
    error('PFAL:LOAD_SCATTER_LD1_DATASET:MissingFileOnDisk', ...
        'Dataset file for "%s" was not found: %s', dataset_src.id, dataset_src.file);
end

loaded = load(dataset_src.file, dataset_src.variable); % import only the requested variable

% Return the same object that the original hard-coded load statement used.
if ~isfield(loaded, dataset_src.variable)
    error('PFAL:LOAD_SCATTER_LD1_DATASET:MissingVariableInFile', ...
        'Variable "%s" was not found in dataset file: %s', ...
        dataset_src.variable, dataset_src.file);
end

data = loaded.(dataset_src.variable);

end

function text_value = require_text_field(src, field_name, field_label)
%REQUIRE_TEXT_FIELD Read and validate a non-empty text field from a struct.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_SCATTER_LD1_DATASET:MissingField', ...
        'The dataset configuration is missing the %s.', field_label);
end

text_value = char(src.(field_name));

end

function root_dir = repo_root()
%REPO_ROOT Resolve the repository root from the current package location.

utils_dir = fileparts(mfilename('fullpath')); % .../+UTILS
root_dir = fileparts(utils_dir); % repository root above the package folder

end

function dataset_file = resolve_dataset_path(dataset_file, config_dir)
%RESOLVE_DATASET_PATH Expand relative paths from the config directory.

if is_absolute_path(dataset_file)
    return % already fully resolved
end

dataset_file = fullfile(config_dir, dataset_file);

end

function tf = is_absolute_path(path_value)
%IS_ABSOLUTE_PATH Detect Windows drive, UNC, and POSIX absolute paths.

tf = ~isempty(regexp(path_value, '^[A-Za-z]:[\\/]', 'once')) || ...
    strncmp(path_value, '\\', 2) || ...
    strncmp(path_value, '//', 2) || ...
    strncmp(path_value, '/', 1); % cover Windows drive letters, UNC paths, and POSIX roots

end
