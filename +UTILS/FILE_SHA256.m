function hash_text = FILE_SHA256(file_path)
%FILE_SHA256 Compute a lowercase SHA-256 digest for a local file.
%   HASH_TEXT = UTILS.FILE_SHA256(FILE_PATH) reads the file in bounded
%   chunks so large experimental artifacts do not need to fit in memory.

if nargin < 1 || isempty(file_path) || ~isfile(file_path)
    error('PFAL:FILE_SHA256:MissingFile', ...
        'Cannot compute SHA-256 because the file does not exist: %s', ...
        char(string(file_path)));
end

digest = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(file_path, 'r');
if fid < 0
    error('PFAL:FILE_SHA256:OpenFailed', ...
        'Could not open the file for hashing: %s', file_path);
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
