function uu=ang_mapping(S,k)
%function uu=ang_mapping(S,k,use_gen)
% this is outdated always want to use the generalized version. 
  
  n=size(S,2);  

  % Need to ensure that the similarity is zero at the diagonal. (as in the paper).
  % S=makeDiagonalZero(S);

  %  Compute Laplacian L=D^-1/2 S D^-1/2
  D=diag(sum(S));
  Dsqrt = sqrt(D);
%  if use_gen
	[vv ignore_lambda indices] = myeigs_gen( S , D, k );
	%uu = DSqrt vv.
	uu=Dsqrt*vv; 
	
%   else
% 	%	L=S_to_L(S); 
% 	L= Dsqrt \ S /Dsqrt; 
	
% 	%  Compute eigenvectors
	
% 	[ uu lang ] = myeigs( L, k );
% 	% [uu lang] = eigsort(L); 
% 	%uu = uu( :, 1:k ); % the first k eigen vectors 		
	
%   end
  uu=uu';  % make it k x n so that each point becomes a k-dim column vector 
  
  % normalize the uu along 2nd dimension. (so that each point is on the
  % unit sphere now. 
  uu=normalize_2nd(uu); 

