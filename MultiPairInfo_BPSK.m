clc; clear variables; close all;

N = 5*10^5;

Pt = 0:2:40;			
pt = (10^-3)*db2pow(Pt);	

BW = 10^6;			
No = -174 + 10*log10(BW);				%
no = (10^-3)*db2pow(No);	 

d1 = 500; d2 = 200; d3 = 70; d4 = 500; d5 = 200; d6 = 70;	%Distances
a1 = 0.90; a2 = 0.08; a3 = 0.02; a4 = 0.90; a5 = 0.08; a6 = 0.02;	%Power allocation coefficients

eta = 4;	%Path loss exponent

%Generate Rayleigh fading channel for the six users
h1 = sqrt(d1^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h2 = sqrt(d2^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h3 = sqrt(d3^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h4 = sqrt(d4^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h5 = sqrt(d5^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h6 = sqrt(d6^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);

%Generate noise samples for the six users
n1 = sqrt(no)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
n2 = sqrt(no)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
n3 = sqrt(no)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
n4 = sqrt(no)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
n5 = sqrt(no)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
n6 = sqrt(no)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);

%Generate random binary message data for the six users
x1 = randi([0 1],N,1);
x2 = randi([0 1],N,1);
x3 = randi([0 1],N,1);
x4 = randi([0 1],N,1);
x5 = randi([0 1],N,1);
x6 = randi([0 1],N,1);

figure(1)
subplot(6,1,1);
stem(x1, 'linewidth',3), grid on;
title('  Information before Transmiting - Channel 1 ');
axis([ 0 20 0 1.5]);
%plot(tt,y_in,'linewidth',3), grid on;

subplot(6,1,2);
stem(x2, 'linewidth',3), grid on;
title('  Information before Transmiting - Channel 2 ');
axis([ 0 20 0 1.5]);

subplot(6,1,3);
stem(x3, 'linewidth',3), grid on;
title('  Information before Transmiting - Channel 3 ');
axis([ 0 20 0 1.5]);

subplot(6,1,4);
stem(x4, 'linewidth',3), grid on;
title('  Information before Transmiting - Channel 4 ');
axis([ 0 20 0 1.5]);

subplot(6,1,5);
stem(x5, 'linewidth',3), grid on;
title('  Information before Transmiting - Channel 5 ');
axis([ 0 20 0 1.5]);

subplot(6,1,6);
stem(x6, 'linewidth',3), grid on;
title('  Information before Transmiting - Channel 6 ');
axis([ 0 20 0 1.5]);


%Create QPSKModulator and QPSKDemodulator objects
QPSKmod = comm.BPSKModulator;%('BitInput',true); 
%display(QPSKmod)
QPSKdemod = comm.BPSKDemodulator;%('BitOutput',true); 

%Perform QPSK modulation
xmod1 = step(QPSKmod, x1);
xmod2 = step(QPSKmod, x2);
xmod3 = step(QPSKmod, x3);
xmod4 = step(QPSKmod, x4);
xmod5 = step(QPSKmod, x5);
xmod6 = step(QPSKmod, x6);

x = sqrt(a1)*xmod1 + sqrt(a2)*xmod2 + sqrt(a3)*xmod3;
%disp(size(x))
x1_send = sqrt(a1)*xmod1;
x2_send = sqrt(a2)*xmod2;
x3_send = sqrt(a3)*xmod3;
x_2 = sqrt(a4)*xmod4 + sqrt(a5)*xmod5 + sqrt(a6)*xmod6;
x4_send = sqrt(a4)*xmod4;
x5_send = sqrt(a5)*xmod5;
x6_send = sqrt(a6)*xmod6;

for u = 1:length(Pt)
	
%Received signals
   y1 = sqrt(pt(u))*x.*h1 + n1;	%At user 1
   y2 = sqrt(pt(u))*x.*h2 + n2;	%At user 2
   y3 = sqrt(pt(u))*x.*h3 + n3;	%At user 3
   
   eq1 = y1./h1;
   eq2 = y2./h2;
   eq3 = y3./h3;
   
%Decode at user 1 (Direct decoding)
   dec1 = step(QPSKdemod, eq1);
   
%Decode at user 2
   dec12 = step(QPSKdemod, eq2);		
   dec12_remod = step(QPSKmod, dec12);		
   rem2 = eq2 - sqrt(a1*pt(u))*dec12_remod;	%SIC to remove U1's data
   dec2 = step(QPSKdemod, rem2);		
   
%Decode at user 3
   dec13 = step(QPSKdemod, eq3);		%Direct demodulation to get U1's data
   dec13_remod = step(QPSKmod, dec13);		%Remodulation of U1's data
   rem31 = eq3 - sqrt(a1*pt(u))*dec12_remod;	%SIC to remove U1's data
   dec23 = step(QPSKdemod, rem31);		%Direct demodulation of remaining signal to get U2's data
   dec23_remod = step(QPSKmod, dec23);		%Remodulation of U2's data
   rem3 = rem31 - sqrt(a2*pt(u))*dec23_remod;	%SIC to remove U2's data
   dec3 = step(QPSKdemod, rem3);		%Demodulate remaining signal to get U3's data
   
%BER calculation
   ber1(u) = biterr(dec1, x1)/N;
   ber2(u) = biterr(dec2, x2)/N;
   ber3(u) = biterr(dec3, x3)/N;
end

for u = 1:length(Pt)
	
%Received signals
   y4 = sqrt(pt(u))*x_2.*h4 + n4;	%At user 1
   y5 = sqrt(pt(u))*x_2.*h5 + n5;	%At user 2
   y6 = sqrt(pt(u))*x_2.*h6 + n6;	%At user 3
   
   eq4 = y4./h4;
   eq5 = y5./h5;
   eq6 = y6./h6;
   
%Decode at user 1 (Direct decoding)
   dec4 = step(QPSKdemod, eq4);
   
%Decode at user 2
   dec45 = step(QPSKdemod, eq5);		
   dec45_remod = step(QPSKmod, dec45);		
   rem5 = eq5 - sqrt(a4*pt(u))*dec45_remod;	%SIC to remove U1's data
   dec5 = step(QPSKdemod, rem5);		
   
%Decode at user 3
   dec46 = step(QPSKdemod, eq6);		%Direct demodulation to get U1's data
   dec46_remod = step(QPSKmod, dec46);		%Remodulation of U1's data
   rem64 = eq6 - sqrt(a4*pt(u))*dec45_remod;	%SIC to remove U1's data
   dec56 = step(QPSKdemod, rem64);		%Direct demodulation of remaining signal to get U2's data
   dec56_remod = step(QPSKmod, dec56);		%Remodulation of U2's data
   rem6 = rem64 - sqrt(a5*pt(u))*dec56_remod;	%SIC to remove U2's data
   dec6 = step(QPSKdemod, rem6);		%Demodulate remaining signal to get U3's data
   
%BER calculation
   ber4(u) = biterr(dec4, x4)/N;
   ber5(u) = biterr(dec5, x5)/N;
   ber6(u) = biterr(dec6, x6)/N;
end

C1 = xor(dec1,dec4);
C2 = xor(dec2,dec5);
C3 = xor(dec3,dec6);
%release(QPSKmod);
QPSKmod = comm.BPSKModulator;%('BitInput',true);
%class(QPSKmod)
QPSKdemod = comm.BPSKDemodulator;%('BitOutput',true); 

cmod1 = step(QPSKmod, C1);
cmod2 = step(QPSKmod, C2);
cmod3 = step(QPSKmod, C3);
%disp(size(C))

C = sqrt(a1)*cmod1 + sqrt(a2)*cmod2 + sqrt(a3)*cmod3;

for u = 1:length(Pt)
	
%Received signals
   y1_ = sqrt(pt(u))*C.*h1 + n1;	%At user 1
   %display(y1)
   y2_ = sqrt(pt(u))*C.*h2 + n2;	%At user 2
   y3_ = sqrt(pt(u))*C.*h3 + n3;	%At user 3
   
   eq1_ = y1./h1;
   eq2_ = y2./h2;
   eq3_ = y3./h3;
   
%Decode at user 1 (Direct decoding)
   dec1_ = step(QPSKdemod, eq1_);
   received_1 = or(dec1_,x1);
   
%Decode at user 2
   dec12_ = step(QPSKdemod, eq2_);	
   release(QPSKmod);
   dec12_remod_ = step(QPSKmod, dec12_);		
   rem2_ = eq2_ - sqrt(a1*pt(u))*dec12_remod_;	%SIC to remove U1's data
   dec2_ = step(QPSKdemod, rem2_);
   received_2 = or(dec2_,x2);
   
%Decode at user 3
   dec13_ = step(QPSKdemod, eq3_);		%Direct demodulation to get U1's data
   dec13_remod_ = step(QPSKmod, dec13_);		%Remodulation of U1's data
   rem31_ = eq3_ - sqrt(a1*pt(u))*dec12_remod_;	%SIC to remove U1's data
   dec23_ = step(QPSKdemod, rem31_);		%Direct demodulation of remaining signal to get U2's data
   dec23_remod_ = step(QPSKmod, dec23_);		%Remodulation of U2's data
   rem3_ = rem31_ - sqrt(a2*pt(u))*dec23_remod_;	%SIC to remove U2's data
   dec3_ = step(QPSKdemod, rem3_);		%Demodulate remaining signal to get U3's data
   received_3 = or(dec3_,x3);
   
%BER calculation
   ber1_(u) = biterr(dec1_, C1)/N;
   ber2_(u) = biterr(dec2_, C2)/N;
   ber3_(u) = biterr(dec3_, C3)/N;
end

for u = 1:length(Pt)
	
%Received signals
   y4_ = sqrt(pt(u))*C.*h4 + n4;	%At user 1
   y5_ = sqrt(pt(u))*C.*h5 + n5;	%At user 2
   y6_ = sqrt(pt(u))*C.*h6 + n6;	%At user 3
   
   eq4_ = y4./h4;
   eq5_ = y5./h5;
   eq6_ = y6./h6;
   
%Decode at user 1 (Direct decoding)
   dec4_ = step(QPSKdemod, eq4_);
   received_4 = or(dec4_,x4);
   
%Decode at user 2
   dec45_ = step(QPSKdemod, eq5_);	
   release(QPSKmod);
   dec45_remod_ = step(QPSKmod, dec45_);		
   rem5_ = eq5_ - sqrt(a4*pt(u))*dec45_remod_;	%SIC to remove U1's data
   dec5_ = step(QPSKdemod, rem5_);
   received_5 = or(dec5_,x5);
   
%Decode at user 3
   dec46_ = step(QPSKdemod, eq6_);		%Direct demodulation to get U1's data
   dec46_remod_ = step(QPSKmod, dec46_);		%Remodulation of U1's data
   rem64_ = eq6_ - sqrt(a4*pt(u))*dec45_remod_;	%SIC to remove U1's data
   dec56_ = step(QPSKdemod, rem64_);		%Direct demodulation of remaining signal to get U2's data
   dec56_remod_ = step(QPSKmod, dec56_);		%Remodulation of U2's data
   rem6_ = rem64_ - sqrt(a5*pt(u))*dec56_remod_;	%SIC to remove U2's data
   dec6_ = step(QPSKdemod, rem6_);		%Demodulate remaining signal to get U3's data
   received_6 = or(dec6_,x6);
   
%BER calculation
   ber4_(u) = biterr(dec4_, C1)/N;
   ber5_(u) = biterr(dec5_, C2)/N;
   ber6_(u) = biterr(dec6_, C3)/N;
end

figure(2)
subplot(6,1,1);
stem(received_1, 'linewidth',3), grid on;
title('  Information after Transmiting - Channel 1 ');
axis([ 0 20 0 1.5]);
%plot(tt,y_in,'linewidth',3), grid on;

subplot(6,1,2);
stem(received_2, 'linewidth',3), grid on;
title('  Information after Transmiting - Channel 2 ');
axis([ 0 20 0 1.5]);

subplot(6,1,3);
stem(received_3, 'linewidth',3), grid on;
title('  Information after Transmiting - Channel 3 ');
axis([ 0 20 0 1.5]);

subplot(6,1,4);
stem(received_4, 'linewidth',3), grid on;
title('  Information after Transmiting - Channel 4 ');
axis([ 0 20 0 1.5]);

subplot(6,1,5);
stem(received_5, 'linewidth',3), grid on;
title('  Information after Transmiting - Channel 5 ');
axis([ 0 20 0 1.5]);

subplot(6,1,6);
stem(received_6, 'linewidth',3), grid on;
title('  Information after Transmiting - Channel 6 ');
axis([ 0 20 0 1.5]);

figure(3)
semilogy(Pt, ber1, '-o', 'linewidth', 2); hold on; grid on;
semilogy(Pt, ber2, '-o', 'linewidth', 2);
semilogy(Pt, ber3, '-o', 'linewidth', 2);
semilogy(Pt, ber4, '-o', 'linewidth', 2);
semilogy(Pt, ber5, '-o', 'linewidth', 2);
semilogy(Pt, ber6, '-o', 'linewidth', 2);

xlabel('Transmit power (dBm)');
ylabel('BER');
legend('User 1 (Weakest user)', 'User 2', 'User 3 (Strongest user)','User 4 (Weakest user)', 'User 5', 'User 6 (Strongest user)');

figure(4)
semilogy(Pt, ber1_, '-o', 'linewidth', 2); hold on; grid on;
semilogy(Pt, ber2_, '-o', 'linewidth', 2);
semilogy(Pt, ber3_, '-o', 'linewidth', 2);




