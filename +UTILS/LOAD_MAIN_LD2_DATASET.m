function [data, dataset_src] = LOAD_MAIN_LD2_DATASET(dataset_cfg)
%LOAD_MAIN_LD2_DATASET Load the scaled aggregate library used by main_LD2.
%   DATA = UTILS.LOAD_MAIN_LD2_DATASET(DATASET_CFG) loads the MAT-file
%   variable declared in the dataset section of the main_LD2 config.
%
%   [DATA, DATASET_SRC] = UTILS.LOAD_MAIN_LD2_DATASET(...) also returns the
%   normalized dataset metadata that was used during the load.

if nargin < 1 || ~isstruct(dataset_cfg)
    error('PFAL:LOAD_MAIN_LD2_DATASET:MissingDatasetConfig', ...
        'The LD2 dataset configuration must be provided as a struct.');
end

dataset_src.id = require_text_field(dataset_cfg, 'id', 'dataset identifier');
dataset_src.file = require_text_field(dataset_cfg, 'file', 'dataset file path');
dataset_src.variable = require_text_field(dataset_cfg, 'variable', 'dataset variable name');
dataset_src.url = optional_text_field(dataset_cfg, 'url');
dataset_src.doi = optional_text_field(dataset_cfg, 'doi');
dataset_src.sha256 = optional_text_field(dataset_cfg, 'sha256');
dataset_src.cache_file = optional_text_field(dataset_cfg, 'cache_file');

if isfield(dataset_cfg, 'config_dir') && ~isempty(dataset_cfg.config_dir)
    dataset_src.config_dir = char(dataset_cfg.config_dir);
else
    dataset_src.config_dir = fullfile(repo_root(), 'config');
end

dataset_src.file = resolve_dataset_path(dataset_src.file, dataset_src.config_dir);
if ~isempty(dataset_src.cache_file)
    dataset_src.cache_file = resolve_dataset_path(dataset_src.cache_file, ...
        dataset_src.config_dir);
end

if ~isfile(dataset_src.file)
    dataset_src.file = resolve_missing_dataset_file(dataset_src);
end

verify_dataset_checksum(dataset_src.file, dataset_src.sha256, dataset_src.id);
loaded = load(dataset_src.file, dataset_src.variable);

if ~isfield(loaded, dataset_src.variable)
    error('PFAL:LOAD_MAIN_LD2_DATASET:MissingVariableInFile', ...
        'Variable "%s" was not found in LD2 dataset file: %s', ...
        dataset_src.variable, dataset_src.file);
end

data = loaded.(dataset_src.variable);

end

function text_value = require_text_field(src, field_name, field_label)
%REQUIRE_TEXT_FIELD Read and validate a non-empty text field from a struct.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_LD2_DATASET:MissingField', ...
        'The LD2 dataset configuration is missing the %s.', field_label);
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

function dataset_file = resolve_missing_dataset_file(dataset_src)
%RESOLVE_MISSING_DATASET_FILE Download a missing MAT file when URL is set.

if isempty(dataset_src.url)
    error('PFAL:LOAD_MAIN_LD2_DATASET:MissingFileOnDisk', ...
        ['LD2 dataset file for "%s" was not found: %s\n' ...
        'Run main_scatter_v8 or set dataset.url in the LD2 config.'], ...
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
    error('PFAL:LOAD_MAIN_LD2_DATASET:DownloadFailed', ...
        'Could not download LD2 dataset "%s" from %s: %s', ...
        dataset_src.id, dataset_src.url, err.message);
end

end

function dataset_file = default_cache_file(configured_file)
%DEFAULT_CACHE_FILE Put remote MAT downloads in data/main_ld2 by default.

[~, file_name, file_ext] = fileparts(configured_file);
dataset_file = fullfile(repo_root(), 'data', 'main_ld2', ...
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
    error('PFAL:LOAD_MAIN_LD2_DATASET:ChecksumMismatch', ...
        'SHA-256 mismatch for LD2 dataset "%s".\nExpected: %s\nActual:   %s', ...
        dataset_id, expected_sha256, actual_sha256);
end

end

function hash_text = file_sha256(file_path)
%FILE_SHA256 Compute SHA-256 for a file without external dependencies.

digest = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(file_path, 'r');
if fid < 0
    error('PFAL:LOAD_MAIN_LD2_DATASET:HashReadFailed', ...
        'Could not open LD2 dataset file for hashing: %s', file_path);
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
    dataset_file = normalize_path(dataset_file);
    return
end

dataset_file = normalize_path(fullfile(config_dir, dataset_file));

end

function tf = is_absolute_path(path_value)
%IS_ABSOLUTE_PATH Detect Windows drive, UNC, and POSIX absolute paths.

tf = ~isempty(regexp(path_value, '^[A-Za-z]:[\\/]', 'once')) || ...
    strncmp(path_value, '\\', 2) || ...
    strncmp(path_value, '//', 2) || ...
    strncmp(path_value, '/', 1);

end

function path_value = normalize_path(path_value)
%NORMALIZE_PATH Collapse relative path segments when MATLAB has Java.

try
    path_value = char(java.io.File(path_value).getCanonicalPath());
catch
    path_value = char(path_value);
end

end

function root_dir = repo_root()
%REPO_ROOT Resolve the repository root from the current package location.

utils_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(utils_dir);

end
