function RASTERIZE_PDF(pdf_file, png_file, resolution_dpi)
%RASTERIZE_PDF Create a PNG that matches a vector PDF exactly.
%   UTILS.RASTERIZE_PDF(PDF_FILE, PNG_FILE, RESOLUTION_DPI) invokes
%   pdftocairo to render the first and only page of a publication figure.
%   Rendering from the final PDF prevents MATLAB's raster renderer from
%   substituting fonts that are embedded correctly in the vector output.

pdf_file = char(pdf_file);
png_file = char(png_file);
resolution_dpi = double(resolution_dpi);
if ~isfile(pdf_file)
    error('PFAL:RASTERIZE_PDF:MissingPDF', ...
        'The source PDF does not exist: %s', pdf_file);
end
if ~isscalar(resolution_dpi) || ~isfinite(resolution_dpi) || ...
        resolution_dpi <= 0 || resolution_dpi ~= round(resolution_dpi)
    error('PFAL:RASTERIZE_PDF:InvalidResolution', ...
        'The PNG resolution must be a positive integer in dots per inch.');
end
[png_dir, png_name, png_extension] = fileparts(png_file);
if ~strcmpi(png_extension, '.png')
    error('PFAL:RASTERIZE_PDF:InvalidOutputExtension', ...
        'The output file must use the .png extension: %s', png_file);
end
if isempty(png_dir)
    png_dir = pwd;
end
if ~isfolder(png_dir)
    mkdir(png_dir);
end
if contains(pdf_file, '"') || contains(png_dir, '"')
    error('PFAL:RASTERIZE_PDF:UnsupportedPath', ...
        'PDF rasterization paths cannot contain double-quote characters.');
end

% A temporary destination prevents a failed renderer invocation from being
% mistaken for a successful export when an older PNG already exists.
temporary_prefix = tempname(png_dir);
temporary_png = [temporary_prefix, '.png'];
temporary_cleanup = onCleanup(@() delete_if_present(temporary_png)); %#ok<NASGU>
command = sprintf(['pdftocairo -png -singlefile -r %d ', ...
    '"%s" "%s"'], resolution_dpi, pdf_file, temporary_prefix);
[status, output] = system(command);
if status ~= 0 || ~isfile(temporary_png)
    error('PFAL:RASTERIZE_PDF:RendererFailed', ...
        ['pdftocairo could not create the PNG. Install Poppler or a ', ...
        'TeX distribution that provides pdftocairo, or set png_source ', ...
        'to "figure". Renderer output:\n%s'], strtrim(output));
end

final_png = fullfile(png_dir, [png_name, '.png']);
[moved, message] = movefile(temporary_png, final_png, 'f');
if ~moved
    error('PFAL:RASTERIZE_PDF:MoveFailed', ...
        'Could not replace the PNG export "%s": %s', final_png, message);
end

end

function delete_if_present(file_path)
%DELETE_IF_PRESENT Remove only the temporary raster produced by this call.

if isfile(file_path)
    delete(file_path)
end

end
