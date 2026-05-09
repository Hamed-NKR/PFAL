function [params_ud, params_const] = LD2_PARAMS_FROM_CONFIG(params_cfg)
%LD2_PARAMS_FROM_CONFIG Build LD2 transport parameter tables from JSON data.
%   [PARAMS_UD, PARAMS_CONST] = UTILS.LD2_PARAMS_FROM_CONFIG(PARAMS_CFG)
%   converts the transport.user_defined JSON rows into the table shape used
%   by TRANSP.INIT_DOM, PAR.INIT_LOC, and the transport routines.

if nargin < 1 || ~isstruct(params_cfg)
    error('PFAL:LD2_PARAMS_FROM_CONFIG:MissingParams', ...
        'transport.user_defined must be a struct array.');
end

params_cfg = params_cfg(:);
n_rows = numel(params_cfg);

Name = cell(n_rows, 1);
Value = zeros(n_rows, 1);
Unit = cell(n_rows, 1);
Description = cell(n_rows, 1);

for i = 1 : n_rows
    Name{i} = require_text_field(params_cfg(i), 'name', ...
        sprintf('transport.user_defined(%d).name', i));
    Value(i) = require_scalar_field(params_cfg(i), 'value', ...
        sprintf('transport.user_defined(%d).value', i));
    Unit{i} = require_text_field(params_cfg(i), 'unit', ...
        sprintf('transport.user_defined(%d).unit', i));
    Description{i} = require_text_field(params_cfg(i), 'description', ...
        sprintf('transport.user_defined(%d).description', i));
end

expected_names = {'volf'; 'dom_size(1)'; 'dom_size(2)'; 'dom_size(3)'; ...
    'n_par'; 'n_pp(1)'; 'n_pp(2)'; 'd_pp(1)'; 'd_pp(2)'; 'd_pp(3)'; ...
    'temp_f'; 'v_f(1)'; 'v_f(2)'; 'v_f(3)'; 'p_f'};
if n_rows ~= numel(expected_names) || any(~strcmp(Name, expected_names))
    error('PFAL:LD2_PARAMS_FROM_CONFIG:InvalidParamOrder', ...
        ['transport.user_defined must contain the same 15 rows as inputs/LD2_Params.txt ' ...
        'in the same order.']);
end

params_ud = table(Name, Value, Unit, Description);

Name = {'rho_bc'; 'M_air'; 'kb'; 'Na'; 'Ru'};
Value = [1.86e3; 28.97e-3; 1.381e-23; 6.022e23; 8.314];
Unit = {'kg/m3'; 'kg/mol'; 'j/k'; 'mol^-1'; 'j/mol.k'};
Description = {'Black Carbon bulk density'; 'Air molar mass'; ...
    'Boltzmann constant'; 'Avogadro constant'; 'Universal gas constant'};
params_const = table(Name, Value, Unit, Description);

end

function out = require_text_field(src, field_name, label)
%REQUIRE_TEXT_FIELD Read and validate a non-empty text field.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LD2_PARAMS_FROM_CONFIG:MissingText', ...
        'The LD2 config is missing %s.', label);
end

out = char(src.(field_name));

end

function out = require_scalar_field(src, field_name, label)
%REQUIRE_SCALAR_FIELD Read and validate a numeric scalar field.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LD2_PARAMS_FROM_CONFIG:MissingScalar', ...
        'The LD2 config is missing %s.', label);
end

out = double(src.(field_name));
if ~isscalar(out) || ~isfinite(out)
    error('PFAL:LD2_PARAMS_FROM_CONFIG:InvalidScalar', ...
        'The LD2 config field %s must be a finite numeric scalar.', label);
end

end
