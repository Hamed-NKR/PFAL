function [loaded, source] = LOAD_VALIDATION_SOURCE(source_cfg)
%LOAD_VALIDATION_SOURCE Resolve, verify, and load one validation data source.
%   The loader prefers the configured local file, then an existing cache,
%   and finally the configured URL. Only explicitly listed MAT variables are
%   loaded, which keeps the effective-density workspace out of memory.

if nargin < 1 || ~isstruct(source_cfg)
    error('PFAL:LOAD_VALIDATION_SOURCE:MissingConfig', ...
        'A validation source configuration struct is required.');
end

source = source_cfg;
required_fields = {'id', 'file', 'variables'};
for i = 1:numel(required_fields)
    field_name = required_fields{i};
    if ~isfield(source, field_name) || isempty(source.(field_name))
        error('PFAL:LOAD_VALIDATION_SOURCE:MissingField', ...
            'Validation source config is missing "%s".', field_name);
    end
end

source.id = char(source.id);
source.file = char(source.file);
source.url = optional_text(source, 'url');
source.doi = optional_text(source, 'doi');
source.sha256 = optional_text(source, 'sha256');
source.cache_file = optional_text(source, 'cache_file');
source.variables = normalize_variable_names(source.variables);

resolved_file = source.file;
if ~isfile(resolved_file) && ~isempty(source.cache_file) && ...
        isfile(source.cache_file)
    resolved_file = source.cache_file;
end

if ~isfile(resolved_file)
    resolved_file = download_source(source);
end

if ~isempty(source.sha256)
    actual_hash = UTILS.FILE_SHA256(resolved_file);
    if ~strcmpi(actual_hash, source.sha256)
        error('PFAL:LOAD_VALIDATION_SOURCE:ChecksumMismatch', ...
            ['SHA-256 mismatch for validation source "%s".\n' ...
             'Expected: %s\nActual:   %s'], ...
            source.id, source.sha256, actual_hash);
    end
end

loaded = load(resolved_file, source.variables{:});
for i = 1:numel(source.variables)
    variable_name = source.variables{i};
    if ~isfield(loaded, variable_name)
        error('PFAL:LOAD_VALIDATION_SOURCE:MissingVariable', ...
            'Variable "%s" was not found in validation source "%s": %s', ...
            variable_name, source.id, resolved_file);
    end
end

source.resolved_file = resolved_file;
source.loaded_variables = source.variables;

end

function variable_names = normalize_variable_names(value)
%NORMALIZE_VARIABLE_NAMES Convert decoded JSON text arrays to cell strings.

if ischar(value)
    variable_names = {value};
elseif isstring(value)
    variable_names = cellstr(value(:));
elseif iscell(value)
    variable_names = cellfun(@char, value(:), 'UniformOutput', false);
else
    error('PFAL:LOAD_VALIDATION_SOURCE:InvalidVariables', ...
        'source.variables must be text or an array of text values.');
end

end

function resolved_file = download_source(source)
%DOWNLOAD_SOURCE Fetch a missing source into its configured cache location.

if isempty(source.url)
    error('PFAL:LOAD_VALIDATION_SOURCE:MissingFile', ...
        ['Validation source "%s" was not found locally: %s\n' ...
         'Provide the file or configure a URL and cache_file.'], ...
        source.id, source.file);
end

if isempty(source.cache_file)
    resolved_file = source.file;
else
    resolved_file = source.cache_file;
end

parent_dir = fileparts(resolved_file);
if ~isempty(parent_dir) && ~isfolder(parent_dir)
    mkdir(parent_dir);
end

try
    websave(resolved_file, source.url);
catch err
    error('PFAL:LOAD_VALIDATION_SOURCE:DownloadFailed', ...
        'Could not download validation source "%s" from %s: %s', ...
        source.id, source.url, err.message);
end

end

function value = optional_text(source, field_name)
%OPTIONAL_TEXT Normalize optional source metadata to a character vector.

value = '';
if isfield(source, field_name) && ~isempty(source.(field_name))
    value = char(source.(field_name));
end

end
