function assignment=cluster_point_kmeans(xx,k,num_ortho,num_random) 
  global KMEANS_THRESHOLD KMEANS_MAX_ITER
  threshold=KMEANS_THRESHOLD;
  max_iter=KMEANS_MAX_ITER;
  
  [center, assignment, distortion]=kmeans_ortho_multiple(xx,k,num_ortho,num_random,threshold,max_iter);
  
  