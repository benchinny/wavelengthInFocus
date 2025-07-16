%%

s = 0:0.01:6;

G_fast = 8;

G_slow = 5;

tau_fast = 2.5;

tau_slow = 100;

figure;
hold on;
plot(s,G_fast./(1+tau_fast.*s),'-','LineWidth',2);
plot(s,G_slow./(1+tau_slow.*s),'-','LineWidth',2);
plot([0 6],[1 1],'-','LineWidth',2);
axis square;
set(gca,'FontSize',15);
xlabel('Frequency');
ylabel('Amplitude');
legend('Fast integrator','Slow integrator','Proportional');

t = 0:0.001:1;

figure;
hold on;
plot(t,sin(2.*pi.*2.*t),'-','LineWidth',2);
plot(t,sin(2.*pi.*2.*t+pi),'-','LineWidth',2);
plot([0 1],[0 0],'-','LineWidth',2);
axis square;
set(gca,'FontSize',15);
xlabel('Time (s)');
ylabel('Accommodation (D)');
legend('Stimulus','Response','Doing nothing');

