
fid = fopen('smallBinary.obj', 'rb');
dim = fread(fid,[1 1], 'int');
H = fread(fid, [2 dim^2], 'double');
H = H(1,:) + 1i*H(2,:);
H = reshape(H,dim,dim);


reso_states = fread(fid, [2 1], 'int');
check = fread(fid, [1 1], 'long');
fclose(fid);
