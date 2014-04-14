
tic
fid = fopen('/net/data2/riklund/jacobi/smallBinary.obj', 'rb');
dim = fread(fid,[1 1], 'int');
U = cell(dim,1);
D = zeros(dim,1);
readtime = 0;
tic

for i=1:dim
    offset = (8+8)*(i-1);  % real + imag double = 8+8 bytes
    fseek( fid, offset, 'cof' );

    data = fread(fid, [2 dim-i+1], 'double');

    D(i) = data(1,1) + 1i*data(2,1);
    t1=tic;
    U{i} = data(1,2:end)' + 1i*data(2,2:end)';
    readtime = readtime + toc(t1);
end

reso_states = fread(fid, [2 1], 'int');
check = fread(fid, [1 1], 'long');
fclose(fid);

ttot=toc
