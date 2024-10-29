function y = frconv(h,u)

h = h(:).';
u = u(:).';

N = length(u);
M = length(h);

NN = N+M-1;
NNN = 2.^nextpow2(NN);

H = fft([h zeros(1,NNN-M)]);
U = fft([u zeros(1,NNN-N)]);

y = ifft(H.*U);
y = y(1:NN);