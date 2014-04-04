function [ H ] = loadmatrix( CHAR )
%loadmatrix: H = loadmatrix('CHAR')
%   'CHAR' is the matrix name,
%   ie matrixA, matrixB.
%   currently 'CHAR' = 'A' or 'B'.

if strcmp(CHAR, 'A') || strcmp(CHAR, 'B')
    dim = 25600;
end

path=strcat('/net/data1/nhqm2014/matrix', CHAR)
pathr=strcat(path, '_real.bin')
pathi=strcat(path, '_imag.bin')


fid_r = fopen(pathr);
fid_i = fopen(pathi);
H = fread(fid_r, [dim dim], 'double') + 1i*fread(fid_i, [dim dim], 'double');
fclose(fid_r);
fclose(fid_i);

end

