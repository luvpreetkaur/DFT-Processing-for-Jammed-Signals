
%DSP-EECE 520 Project
%DFT Processing for Jammed Signals
%Scenario #2 (Comb Filter)
    
    %initializing
    f_h=1.1*10^9;   %highest frequency
    f_l=1*10^9; %lowest frequency
    B=f_h-f_l;  %bandwidth
    g=5*10^6;   %gap spacing of 5Mhz
    Fs=2*B; %sampling rate
    T=1/Fs; %sample spacing
    N=2^10;
    N_zp=4*N;
    j=0;
    t=(0:N-1)*T;
    P_x=0;
    P_jam=0;
    P_z=0;
    
    %desired signal
    gauss=randn(1,N);
    h=[1 1];    %FIR filter with 2 taps
    x=filter(h,1,(gauss));
    X=fftshift(fft(x,N_zp));
    
    %jamming signal
    for f=2*((-B/2+g):g:(B/2-g))
        y=exp(1j*2*pi*f*t);
        j=j+y;
    end
    
    %required signal
    J=fftshift(fft(j,N_zp));
    r=x+j;
    R=fftshift(fft(r,N_zp));
        
    omega=(-N_zp/2:N_zp/2-1)*2*pi/N_zp;
    f=f_l+(B/2)*((omega/pi)+1);
    omega1=(-N/2:N/2-1)*2*pi/N;
    f1=f_l+(B/2)*((omega1/pi)+1);
    
    %windowing
    w=rectwin(N);
    win=transpose(w);
    r2=r.*win;
    R2=fftshift(fft(r2,N_zp));
    
    %Comb Filter Design
    M=20;
    vect=ones(1,M);
    c=1/(M+1)*vect;
    C=fftshift(fft(c,N_zp));
    z=filter(c,1,r);
    Z=fftshift(fft(z,N_zp));
    
    %plots
    figure(1);
    subplot(3,1,1);
    plot(t,20*log10(abs(x)));
    title('Desired signal (time domain)');
    subplot(3,1,2);
    plot(t,20*log10(abs(j)));
    title('Jamming noise (time domain)');
    subplot(3,1,3);
    plot(t,20*log10(abs(r)));
    title('Desired signal+jamming noise (time domain)');
    
    figure(2);
    subplot(2,1,1);
    plot(f,20*log10(abs(X)));
    title('Magnitude spectrum of desired signal');
    subplot(2,1,2);
    plot(f,20*log10(abs(J)));
    title('Magnitude spectrum of jamming signal');
    axis([10^9 1.1*10^9 -100 100]);
   
    figure(3);
    subplot(2,1,1);
    plot(f,20*log10(abs(R)),'g');
    title('Magnitude spectrum of (desired signal+jamming) without windowing');
    subplot(2,1,2);
    plot(f, 20*log10(abs(R2)),'g');
    title('Magnitude spectrum of (desired signal+jamming) with windowing');   
    
    figure(4);
    subplot(3,1,1);
    plot(omega,20*log10(abs(C)),'k');
    title('Magnitude spectrum of comb filter');
    subplot(3,1,2);
    plot(omega,unwrap(angle(C)),'k');
    title('Phase of the comb filter');
    subplot(3,1,3);
    plot((1:M)*2*pi/M,20*log10(abs(c)),'k');
    title('Comb filter (time domain)');
 
    figure(5);
    subplot(2,1,1)
    plot(f,20*log10(abs(Z)));
    title('Spectrum of filtered output');
    subplot(2,1,2)
    plot(t,20*log10(abs(z)));
    title('Time domain view of filtered output');

    %Power computations
    for k=1:1:N
        P_jam=P_jam+(1/N)*(abs(j(k))^2);
        P_x=P_x+(1/N)*(abs(x(k))^2);
        P_z=P_z+(1/N)*(abs(z(k))^2);
    end
    error=10*log10(P_x/P_z);
    disp('Error = ');
    disp(error);