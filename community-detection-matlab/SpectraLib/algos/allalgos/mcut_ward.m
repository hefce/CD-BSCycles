function [assignment_multicut,param_string]=mcut_ward(S,k,orig_ca)

  assignment_multicut=cluster_spectral_general(S,k,'mcut_all','ward'); 
  if nargout>1
	param_string=''; 
  end
  
  
