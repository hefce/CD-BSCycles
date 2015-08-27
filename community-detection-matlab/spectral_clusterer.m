% Ng, A., Jordan, M., and Weiss, Y. (2002). On spectral clustering: analysis and an algorithm. In T. Dietterich,
% S. Becker, and Z. Ghahramani (Eds.), Advances in Neural Information Processing Systems 14 
% (pp. 849  856). MIT Press.

% Asad Ali
% GIK Institute of Engineering Sciences & Technology, Pakistan
% Email: asad_82@yahoo.com

% CONCEPT: Introduced the normalization process of simrix matrix(D-1/2 A D-1/2), 
% eigenvectors orthonormal conversion and clustering by kmeans 
function [IDX,C] = spectral_clusterer(simrix,k)


% compute the degree matrix
for i=1:size(simrix,1)
    D(i,i) = sum(simrix(i,:));
end

% compute the normalized laplacian / simrix matrix (method 1)
%NL1 = D^(-1/2) .* L .* D^(-1/2);
for i=1:size(simrix,1)
    for j=1:size(simrix,2)
        NL1(i,j) = simrix(i,j) / (sqrt(D(i,i)) * sqrt(D(j,j)));  
    end
end

% compute the normalized laplacian (method 2)  eye command is used to
% obtain the identity matrix of size m x n
% NL2 = eye(size(simrix,1),size(simrix,2)) - (D^(-1/2) .* simrix .* D^(-1/2));

% perform the eigen value decomposition
[eigVectors,eigValues] = eig(NL1);

% select k largest eigen vectors

nEigVec = eigVectors(:,(size(eigVectors,1)-(k-1)): size(eigVectors,1));

% construct the normalized matrix U from the obtained eigen vectors
for i=1:size(nEigVec,1)
    n = sqrt(sum(nEigVec(i,:).^2));    
    U(i,:) = nEigVec(i,:) ./ n; 
end

% perform kmeans clustering on the matrix U
[IDX,C] = cluster_point_kmeans(U,k,5,20); 

end

