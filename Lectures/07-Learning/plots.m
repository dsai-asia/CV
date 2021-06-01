
x=1:100;
x1=1:50;
x2=51:100;
y1 = exp(-x/20);

clearplot;
plot(x,y1,'b-;Training set;');
hold on;
title('Training error vs. test error','fontsize',30);
xlabel('Model complexity','fontsize',20);
ylabel('Classification error','fontsize',20)'
axis off;

y21 = 0.5+((50-x1)/1.3).^2/2200;
y22 = 0.5+((50-x2)/1.3).^2/6000;
y2=[y21,y22];
plot(x,y2,'r-;Test set;');

axis off;

print('error.eps', '-FHelvetica:30');

