%%

lambda = (400:700)./1000;

deltaR = ARCrefErrorFromLCA(lambda);

figure;
plot(lambda,deltaR,'k-');