function [params_domain, params_const, fl] = LD2_PARAMS_FROM_CONFIG(domain_cfg, fluid_cfg)
%LD2_PARAMS_FROM_CONFIG Build LD2 runtime parameters from JSON config.
%   [PARAMS_DOMAIN, PARAMS_CONST, FL] = UTILS.LD2_PARAMS_FROM_CONFIG(...)
%   converts the compact LD2 config into the table/struct shapes used by
%   PAR.INIT_LOC, TRANSP.PBC, TRANSP.MOBIL, and TRANSP.MARCH.

if nargin < 1 || ~isstruct(domain_cfg)
    error('PFAL:LD2_PARAMS_FROM_CONFIG:MissingDomain', ...
        'transport.domain must be a struct.');
end
if nargin < 2 || ~isstruct(fluid_cfg)
    error('PFAL:LD2_PARAMS_FROM_CONFIG:MissingFluid', ...
        'transport.fluid must be a struct.');
end

volume_fraction = require_scalar_field(domain_cfg, 'volume_fraction', ...
    'transport.domain.volume_fraction');
domain_size = require_vector_field(domain_cfg, 'size', 3, ...
    'transport.domain.size');
temperature = require_scalar_field(fluid_cfg, 'temperature', ...
    'transport.fluid.temperature');

Name = {'volf'; 'dom_size(1)'; 'dom_size(2)'; 'dom_size(3)'};
Value = [volume_fraction; domain_size(:)];
Unit = {'[-]'; '[m]'; '[m]'; '[m]'};
Description = {'Particle volume fraction'; 'Computational domain length'; ...
    'Computational domain width'; 'Computational domain height'};
params_domain = table(Name, Value, Unit, Description);

Name = {'rho_bc'; 'M_air'; 'kb'; 'Na'; 'Ru'};
Value = [1.86e3; 28.97e-3; 1.381e-23; 6.022e23; 8.314];
Unit = {'kg/m3'; 'kg/mol'; 'j/k'; 'mol^-1'; 'j/mol.k'};
Description = {'Black Carbon bulk density'; 'Air molar mass'; ...
    'Boltzmann constant'; 'Avogadro constant'; 'Universal gas constant'};
params_const = table(Name, Value, Unit, Description);

fl = struct('size', domain_size(:), 'temp', temperature, 'mu', [], ...
    'lambda', []);

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

function out = require_vector_field(src, field_name, expected_len, label)
%REQUIRE_VECTOR_FIELD Read and validate a numeric vector field.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LD2_PARAMS_FROM_CONFIG:MissingVector', ...
        'The LD2 config is missing %s.', label);
end

out = double(src.(field_name));
out = reshape(out, 1, []);
if numel(out) ~= expected_len || any(~isfinite(out))
    error('PFAL:LD2_PARAMS_FROM_CONFIG:InvalidVector', ...
        'The LD2 config field %s must be a finite numeric vector of length %d.', ...
        label, expected_len);
end

end
