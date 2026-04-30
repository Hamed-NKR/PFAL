function TEXTBAR(i, f_back)
%TEXTBAR Print a single-line text progress bar in the command window.
%
%   UTILS.TEXTBAR() or UTILS.TEXTBAR(0) initializes a text progress bar.
%
%   UTILS.TEXTBAR(PCT) updates the bar using a fractional completion
%   between 0 and 1.
%
%   UTILS.TEXTBAR([I, N]) updates the bar using the current loop index I
%   and total count N.
%
%   UTILS.TEXTBAR([I1, N1, I2, N2, ...]) combines nested loop counters
%   into a single overall progress bar.
%
%   UTILS.TEXTBAR(..., F_BACK) controls whether the previous line is
%   overwritten. The default is 1. When F_BACK is 0, the bar starts on a
%   fresh line.

% Parse inputs and apply defaults.
if ~exist('i', 'var') || isempty(i)
    i = 0;
end

if ~exist('f_back', 'var') || isempty(f_back)
    f_back = 1;
end

% Convert the caller input to a fractional completion value.
if i(1) < 0
    pct = 0;
elseif length(i) == 1
    pct = i;
    i(2) = 0;
elseif length(i) == 2
    pct = i(1) / i(2);
else
    n_total = prod(i(2:2:end));
    i_global = i(1);
    for ii = 3 : 2 : length(i)
        i_global = i_global + (i(ii) - 1) * prod(i(2:2:(ii - 1)));
    end
    i = [max(i_global, 0), n_total];
    pct = i(1) / i(2);
end

% The first draw should start a new line rather than overwriting text.
if pct == 0
    f_back = 0;
end

% Clamp the progress fraction to a valid range.
pct = max(0, min(1, pct));

% Format bar dimensions.
n_dot = 20; % number of bar slots
n_xtra = 10; % extra padding for percentage text
n_frac = 2 * length(num2str(i(2))) + 1; % width of the trailing count text
n_str = n_dot + n_xtra + n_frac;

% Move to the start of the current line when overwriting an existing bar.
if f_back
    str_prefix = sprintf('\r');
else
    str_prefix = '';
end

% Format the percentage label.
str_pct = num2str(100 * pct, '%.0f');
str_pct = [repmat(' ', [1, 3 - length(str_pct)]), str_pct];

% Format the bar body using plain ASCII characters for terminal stability.
n_fill = floor(pct * n_dot);
if pct >= 1
    str_fill = repmat('=', [1, n_dot]);
    str_empty = '';
elseif n_fill > 0
    str_fill = [repmat('=', [1, max(n_fill - 1, 0)]), '>'];
    str_empty = repmat(' ', [1, n_dot - n_fill]);
else
    str_fill = '';
    str_empty = repmat(' ', [1, n_dot]);
end

% Format the trailing counter when total count information is available.
if i(2) == 0
    str_frac = '';
else
    str_frac = [num2str(floor(i(1)), '%.0f'), '/', num2str(i(2), '%.0f')];
end
str_frac = [str_frac, repmat(' ', [1, max(n_frac - length(str_frac), 0)])];

% Assemble the full output line.
str_out = [' ', str_pct, '%%', ' |', str_fill, str_empty, '| ', str_frac];
str_out = [str_out, repmat(' ', [1, max(n_str - length(str_out), 0)])];

% Print the progress update. Only terminate with a newline when finished.
if pct < (1 - eps)
    fprintf('%s%s', str_prefix, str_out);
else
    fprintf('%s%s [DONE]\n', str_prefix, str_out);
end

end
