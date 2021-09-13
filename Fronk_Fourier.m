%% Fronk Fourier Lab
%
%% PART 1: Single Plane Wave
close all
clear all
A = 1;
lambda = 5.5e-7;
k = (2*pi)/lambda;
x  = -2e-6:1e-9:2e-6;
psi = A*exp(1i*k*x); % Set the wave function
figure(1)
plot3(x, real(psi), imag(psi), 'LineWidth',2) % Make a 3D plot
hold on
plot3(x, real(psi), zeros(size(x))-1.1) % Plot the real part only
plot3(x, zeros(size(x))-1.1, imag(psi)) % Plot the imag part only
hold off
grid on
axis([-1e-6 1e-6 -1.1 1.1 -1.1 1.1]) % Set limits for axes
view([-125 30]) % Change the view angle
xlabel('Position', 'Rotation',-30)
ylabel('Real{\psi}', 'Rotation',10)
zlabel('Imag{\psi}')
legend('complex','real','imag');

figure(2) %Create separate figure
hold on
plot(x,real(psi))
plot(x,imag(psi),'r') %Plot real, imag, and amplitude components of psi separately
plot(x,abs(psi),'LineWidth',1)
title('Real and Imaginary components and amplitude vs Position')
xlabel('Position') 
ylabel('Real{\psi} and Imag{\psi} and A') %Set axis labels/title
legend('real','imag','amp');
hold off

%%
% Looking at these two plots, there are two important things to notice.
% Firstly, the real and imaginary components are out of phase with one
% another by pi/2. Meaning when one is at its peak, the other is at 0. This
% is significant because it explains the constant amplitude found in figure
% 2. Since whenever one is at the max amplitude the other is at 0, value of psi 
% never exceeds the set amplitude of 1. 
% This explains the constant amplitude and the constant radius of the complete complex spiral
% in figure 1. 

%% PART 2: Multiple Waves and their Sum

%Three plane waves
n = -1:1; %Create spacing to generate wave numbers 
x2  = -1.5e-5:1e-8:1.5e-5; %Set x bounds and spacing
deltak = .05*k; 
kvals = k + n*deltak;
psi1 = A*exp(1i*kvals(1)*x2);
psi2 = A*exp(1i*kvals(2)*x2); %Create deltak and calculate wave values
psi3 = A*exp(1i*kvals(3)*x2);
figure(3)
hold on
plot(x2,real(psi1),'r')
plot(x2,real(psi2),'g') %Plot individual waves
plot(x2,real(psi3),'b')
title('Real Components of three waves with wavenumber spacing of 5 percent vs Position')
xlabel('Position') 
ylabel('Real{\psi}')
hold off

sum = psi1+psi2+psi3;
amp = sqrt(real(sum).^2+imag(sum).^2); %Sum waves
figure(4)
hold on
plot(x2,real(sum),'r')
plot(x2,amp,'b') %Plot sum of the waves 
title('Sum of 3 waves vs Position') 
xlabel('Position') 
ylabel('Real{\psi} 3 sum')
hold off

%Nine plane waves
n2 = -4:4; %Spacing for the 9 waves
kvals2 = k + n2*deltak; %Wave numbers for the waves 
psimatrix = zeros(9,length(x2));
for j=1:9
    psimatrix(j,:) = A*exp(1i*kvals2(j)*x2); %Populate matrix with psi values
end
j = 0;


figure(5)
for j =1:9
  hold on
  plot(x2,real(psimatrix(j,:)),'color',rand(1,3)); %Plot each wave individually
  hold off
end
title('Real Components of 9 waves with 5 percent deltak spacing vs Position')
xlabel('Position') 
ylabel('Real{\psi}')
j = 0;
sum2 = 0;
for j=1:9
    sum2 = sum2 + psimatrix(j,:);
end
amp2 = sqrt(real(sum2).^2+imag(sum2).^2); 
figure(6)
hold on
plot(x2,real(sum2),'r') %Calculate and plot sum and amplitude of the 9 waves
plot(x2,amp2,'b')
title('Sum of 9 waves vs Position')
xlabel('Position') 
ylabel('Real{\psi} component of sum of 9 waves')
hold off

psiorig = A*exp(1i*k*x2);
figure(7)
hold on
plot(x2,real(sum2),'r')
plot(x2,real(psiorig),'b') %Plot the sum of the 9 waves and the wave with wavenumber kc
title('Sum of 9 waves and original wave vs Position')
xlabel('Position') 
ylabel('Real{\psi} component of sum of 9 waves/original wave')
hold off
j = 0;

%%
% Looking at the graphs of the 3 waves and 9 waves plotted individually,
% zooming in reveals that at certain repeating spots, the individual waves
% fall into phase with one another, and create cohesive peaks. They then
% fall out of phase again until the next occurrence of falling into phase.
% This is visible on the graphs of the sum of the waves too. At these
% positions where we see the waves falling the most into phase, the graph of
% the sum of the waves experiences its highest peaks and amplitudes. These
% sections of high peaks repeat just as the spots where the waves fall into phase 
% on the graphs of the individual waves do. The areas of peak also are more
% compressed on the graph of the sum of 9 waves, as the deltak spacing is
% the same for both the 3 wave sum and the 9 wave sum, but the number of
% waves for the 9 wave sum is increased.

%% PART 3: Wave Packets with Different Wavenumber Spacing
n2 = -4:4;
kvals2 = k + n2*deltak;
deltakvals = [k*.02 k*.05 k*.1]; %Create the values for deltak
kvalsmatrix = zeros(3,9);
for j = 1:3
    kvalsmatrix(j,:) = k+n2*deltakvals(j); %Create a matrix of the kvals
end
j = 0;
psimatrixtwospace = zeros(9,length(x2));
psimatrixfivespace = zeros(9,length(x2)); %Generate matrices of each wave fxn
psimatrixtenspace = zeros(9,length(x2));
for j=1:9
    psimatrixtwospace(j,:) = A*exp(1i*kvalsmatrix(1,j)*x2);
    psimatrixfivespace(j,:) = A*exp(1i*kvalsmatrix(2,j)*x2); %Populate wave matrices
    psimatrixtenspace(j,:) = A*exp(1i*kvalsmatrix(3,j)*x2);
end
j = 0;
sumtwospace = 0;
sumfivespace = 0;
sumtenspace = 0;
for j=1:9
    sumtwospace = sumtwospace + psimatrixtwospace(j,:);
    sumfivespace = sumfivespace + psimatrixfivespace(j,:); %Create arrays of wave sums
    sumtenspace = sumtenspace + psimatrixtenspace(j,:);
end

figure(8)
subplot(3,1,1);
plot(x2,real(sumtwospace));
title('Sum of 9 waves with wavenumber spacing of 2 percent of kcenter')
xlabel('Position') 
ylabel('Real{\psi}')
subplot(3,1,2);
plot(x2,real(sumfivespace));
title('Sum of 9 waves with wavenumber spacing of 5 percent of kcenter') %Plot all sums
xlabel('Position') 
ylabel('Real{\psi}')
subplot(3,1,3);
plot(x2,real(sumtenspace));
title('Sum of 9 waves with wavenumber spacing of 10 percent of kcenter')
xlabel('Position') 
ylabel('Real{\psi}')
j = 0;

%%
% Looking at these plots, changing the value of deltak affects the "spread"
% of the wave, so to speak. For example, the plot of the sum of 9 waves with
% a ten percent spacing has several sets of peaks with an amplitude of 10
% spread out across the entire region. On the other hand, the plot with a
% spacing of 2 percent, while it does have some areas with greater peaks,
% these are much smaller than in the plot with 10 percent spacing. Most of
% the high amplitude is contained in one central pulse in the 2 percent
% plot, and the pattern seems to be that decreasing the spacing between
% wavenumbers causes the wave to coalesce into a more and more uniform
% pulse.

%% PART 4: Different Wavenumber Spacing, Constant Spectral Width
n3 = -8:8;
deltakfive = .05*k;
deltakten = .1*k;
kvalsfivespace = k + n3*deltakfive; %Create kvals
kvalstenspace = k + n2*deltakten;

psimatrixfivespace = zeros(17,length(x2)); %Create wave matrices
psimatrixtenspace = zeros(9,length(x2));
for j=1:17
    psimatrixfivespace(j,:) = A*exp(1i*kvalsfivespace(j)*x2); %Populate wave matrix 
end
j = 0;

for j=1:9
    psimatrixtenspace(j,:) = A*exp(1i*kvalstenspace(j)*x2);%Populate other wave matrix
end
j = 0;
sumfivespace = 0;
sumtenspace = 0;
for j=1:17
    sumfivespace = sumfivespace + psimatrixfivespace(j,:); %Create sum array 
end
j = 0;
for j=1:9
    sumtenspace = sumtenspace + psimatrixtenspace(j,:); %Create other sum array 
end
%Plot sums below
figure(9)
subplot(2,1,1);
plot(x2,real(sumfivespace));
title('Sum of 17 waves with wavenumber spacing of 5 percent of kcenter')
xlabel('Position') 
ylabel('Real{\psi}')
subplot(2,1,2);
plot(x2,real(sumtenspace));
title('Sum of 9 waves with wavenumber spacing of 10 percent of kcenter')
xlabel('Position') 
ylabel('Real{\psi}')
j = 0;

%%
% Looking at these figures, there's a clear relationship between the
% component spacing dk and packet spacing dx. In these two figures, the
% component spacing is changing, and so the spacing between the packets
% themselves are changing. In the plot of 17 waves, the packets are more
% spaced out as compared to the plot of 9 waves. However, the width of the
% packets themselves does not change, as the spectral width is kept
% constant. So, by summing more waves but decreasing the component spacing,
% the spectral width is kept constant, and so the packet width is constant, 
% but the packet spacing is changed. 

%% PART 5: Modulated Plane Wave Amplitude
n2 = -4:4;
deltakvals = [k*.02 k*.05 k*.1]; %Create deltak vals 
kvalsmatrix = zeros(3,9);
for j = 1:3
    kvalsmatrix(j,:) = k+n2*deltakvals(j); %Create kvals
end
j = 0;
Avals = zeros(3,9);
for j=1:3
    Avals(j,:) = exp(-((kvalsmatrix(j,:)-k)/(2*deltakvals(j))).^2); %Create amplitudes
end
j = 0;
psimatrixtwospace = zeros(9,length(x2));
psimatrixfivespace = zeros(9,length(x2)); %Create empty wave matrices
psimatrixtenspace = zeros(9,length(x2));
for j=1:9
    psimatrixtwospace(j,:) = Avals(1,j)*exp(1i*kvalsmatrix(1,j)*x2);
    psimatrixfivespace(j,:) = Avals(2,j)*exp(1i*kvalsmatrix(2,j)*x2); %Populate matrices
    psimatrixtenspace(j,:) = Avals(3,j)*exp(1i*kvalsmatrix(3,j)*x2);
end
j = 0;
sumtwospace = 0;
sumfivespace = 0;
sumtenspace = 0;
for j=1:9
    sumtwospace = sumtwospace + psimatrixtwospace(j,:);
    sumfivespace = sumfivespace + psimatrixfivespace(j,:); %Create wave sum arrays
    sumtenspace = sumtenspace + psimatrixtenspace(j,:);
end
%Plot sums 
figure(10)
subplot(3,1,1);
plot(x2,real(sumtwospace));
title('Sum of 9 waves with wavenumber spacing of 2 percent of kcenter with modulated amplitude ')
xlabel('Position') 
ylabel('Real{\psi}')
subplot(3,1,2);
plot(x2,real(sumfivespace));
title('Sum of 9 waves with wavenumber spacing of 5 percent of kcenter with modulated amplitude')
xlabel('Position') 
ylabel('Real{\psi}')
subplot(3,1,3);
plot(x2,real(sumtenspace));
title('Sum of 9 waves with wavenumber spacing of 10 percent of kcenter with modulated amplitude')
xlabel('Position') 
ylabel('Real{\psi}')
j = 0;
%%
% These waves have a similar shape to those in part 3, but are not entirely
% the same. These waves have little to no fluctuation in between the areas
% of greatest amplitudes, unlike those in part 3. The activity of the waves
% are contained in discrete pulses, separated by areas of 0 amplitude in
% between. There is more complete destructive interference between the
% packets due to the modulation of the amplitude and the containment of the
% packets in a Gaussian envelope. 
