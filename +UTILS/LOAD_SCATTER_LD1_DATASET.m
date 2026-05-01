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
dataset_src.url = optional_text_field(dataset_cfg, 'url');
dataset_src.doi = optional_text_field(dataset_cfg, 'doi');
dataset_src.sha256 = optional_text_field(dataset_cfg, 'sha256');
dataset_src.cache_file = optional_text_field(dataset_cfg, 'cache_file');

% The config loader attaches config_dir so relative dataset paths can be
% resolved from the same folder as the JSON file.
if isfield(dataset_cfg, 'config_dir') && ~isempty(dataset_cfg.config_dir)
    dataset_src.config_dir = char(dataset_cfg.config_dir);
else
    dataset_src.config_dir = fullfile(repo_root(), 'config');
end

% Allow either absolute paths or paths written relative to the config file.
dataset_src.file = resolve_dataset_path(dataset_src.file, dataset_src.config_dir);
if ~isempty(dataset_src.cache_file)
    dataset_src.cache_file = resolve_dataset_path(dataset_src.cache_file, ...
        dataset_src.config_dir);
end

if ~isfile(dataset_src.file)
    dataset_src.file = resolve_missing_dataset_file(dataset_src);
end

verify_dataset_checksum(dataset_src.file, dataset_src.sha256, dataset_src.id);
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

function text_value = optional_text_field(src, field_name)
%OPTIONAL_TEXT_FIELD Read optional text config fields.

text_value = '';
if isfield(src, field_name) && ~isempty(src.(field_name))
    text_value = char(src.(field_name));
end

end

function root_dir = repo_root()
%REPO_ROOT Resolve the repository root from the current package location.

utils_dir = fileparts(mfilename('fullpath')); % .../+UTILS
root_dir = fileparts(utils_dir); % repository root above the package folder

end

function dataset_file = resolve_missing_dataset_file(dataset_src)
%RESOLVE_MISSING_DATASET_FILE Download a missing MAT file when URL is set.

if isempty(dataset_src.url)
    error('PFAL:LOAD_SCATTER_LD1_DATASET:MissingFileOnDisk', ...
        ['Dataset file for "%s" was not found: %s\n' ...
        'Add the file locally or set dataset.url in the scatter config.'], ...
        dataset_src.id, dataset_src.file);
end

if isempty(dataset_src.cache_file)
    dataset_file = default_cache_file(dataset_src.file);
else
    dataset_file = dataset_src.cache_file;
end

ensure_parent_folder(dataset_file);

try
    websave(dataset_file, dataset_src.url);
catch err
    error('PFAL:LOAD_SCATTER_LD1_DATASET:DownloadFailed', ...
        'Could not download dataset "%s" from %s: %s', ...
        dataset_src.id, dataset_src.url, err.message);
end

end

function dataset_file = default_cache_file(configured_file)
%DEFAULT_CACHE_FILE Put remote MAT downloads in data/main_scatter by default.

[~, file_name, file_ext] = fileparts(configured_file);
dataset_file = fullfile(repo_root(), 'data', 'main_scatter', ...
    strcat(file_name, file_ext));

end

function ensure_parent_folder(file_path)
%ENSURE_PARENT_FOLDER Create the destination folder before downloading.

parent_folder = fileparts(file_path);
if ~isfolder(parent_folder)
    mkdir(parent_folder);
end

end

function verify_dataset_checksum(file_path, expected_sha256, dataset_id)
%VERIFY_DATASET_CHECKSUM Validate downloaded/local data when configured.

if isempty(expected_sha256)
    return
end

actual_sha256 = file_sha256(file_path);
if ~strcmpi(actual_sha256, expected_sha256)
    error('PFAL:LOAD_SCATTER_LD1_DATASET:ChecksumMismatch', ...
        ['SHA-256 mismatch for dataset "%s".\nExpected: %s\nActual:   %s'], ...
        dataset_id, expected_sha256, actual_sha256);
end

end

function hash_text = file_sha256(file_path)
%FILE_SHA256 Compute SHA-256 for a file without external dependencies.

digest = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(file_path, 'r');
if fid < 0
    error('PFAL:LOAD_SCATTER_LD1_DATASET:HashReadFailed', ...
        'Could not open dataset file for hashing: %s', file_path);
end
cleanup = onCleanup(@() fclose(fid));

while ~feof(fid)
    chunk = fread(fid, 1024 * 1024, '*uint8');
    if ~isempty(chunk)
        digest.update(typecast(chunk, 'int8'));
    end
end

hash_uint8 = typecast(digest.digest(), 'uint8');
hash_text = lower(reshape(dec2hex(hash_uint8, 2).', 1, []));

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
