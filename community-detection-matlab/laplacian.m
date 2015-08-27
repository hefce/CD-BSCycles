
function L=laplacian(kernelcb)
D=diag(sum(kernelcb));
Dsqrt = sqrt(D);
L= Dsqrt \ kernelcb /Dsqrt;
