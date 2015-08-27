figure;
x=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
subplot(2,2,1)
plot(x,[1,1,1,1,1,0.99,0.7328,0.3684],'-d','MarkerFaceColor', [0 0 0])
%title('Subplot 1: sin(x)')
xlabel('\mu','FontWeight','bold');ylabel('NMI','FontWeight','bold');
text(1.2,3.7, 'exponent','FontSize',8, 'FontWeight','bold');
set(gca,'fontsize',12)
subplot(2,2,2)
plot(x,[1,1,1,1,1,0.95,0.68,0.26],'-d','MarkerFaceColor', [0 0 0])
%title('Subplot 2: sin(2x)');
xlabel('\mu','FontWeight','bold');ylabel('Jaccard Index','FontWeight','bold');

subplot(2,2,3)
plot(x,[1,1,1,1,1,0.99,0.73,0.36],'-d','MarkerFaceColor', [0 0 0])
%title('Subplot 3: sin(4x)');
xlabel('\mu','FontWeight','bold');ylabel('Adjusted Rand Index','FontWeight','bold');

