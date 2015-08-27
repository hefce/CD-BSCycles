figure,
x1=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
y1=[1,1,1,1,1,0.99,0.7328,0.3684];
x2=x1
y2=[1,1,1,1,1,0.95,0.68,0.26]
x3=x1
y3=[1,1,1,1,1,0.99,0.73,0.36]



subplot(2,2,1)
plot(x1,y1,'-d','MarkerFaceColor', [1 0 0],'Markersize',7.34)
xlabel('\mu');ylabel('NMI');
subplot(2,2,2)
plot(x2,y2,'-d','MarkerFaceColor', [1 0 0],'Markersize',7.34)
xlabel('\mu');ylabel('Jaccard Index');
subplot(2,2,3)
plot(x3,y3,'-d','MarkerFaceColor', [1 0 0],'Markersize',7.34)
xlabel('\mu');ylabel('Adjusted Rand Index');

