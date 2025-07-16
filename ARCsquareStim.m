%%

frqCpd = 15;
contrast = 1;
rgbAll = [0.555 0 1];
k0 = 1;
gammaR = 2.4;
gammaG = 2.6;
gammaB = 2.2;

acuStimOrig1 = ARC2Dgabor(smpPos(256,256),[],0,0,[frqCpd 3*frqCpd 5*frqCpd 7*frqCpd], ...
               [contrast contrast/3 contrast/5 contrast/7],-15,90,0.2,0.2, ...
               [rgbAll(k0,1)^gammaR rgbAll(k0,2)^gammaG rgbAll(k0,3)^gammaB],1,1,0,0);

acuStimOrig1(:,:,1) = acuStimOrig1(:,:,1).^(1/gammaR);
acuStimOrig1(:,:,2) = acuStimOrig1(:,:,2).^(1/gammaG);
acuStimOrig1(:,:,3) = acuStimOrig1(:,:,3).^(1/gammaB);