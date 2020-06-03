
%DSP-EECE 520 Project
%DFT Processing for Jammed Signals
%Scenario #1 (Sinusoidal Jamming) - Part A

    %initializing
    f_h=1.1*10^9;   %highest frequency
    f_l=1*10^9; %lowest frequency
    B=f_h-f_l;  %bandwidth
    g=5*10^6;   %gap spacing of 5Mhz
    Fs=2*B; %sampling rate
    T=1/Fs; %sample spacing
    N=2^10;
    N_zp=4*N;   %zero padding
    s=0;
    t=(0:N-1)*T;
    N1=2^16;
    t1=(0:N1-1)*T;
    
    fm=2*((-B/2)+2.5*10^6+(g*randi(18)));   %unknown RF frequency
    x=(10^-3)*exp(1j*2*pi*fm*t);    %unknown RF signal
    for f=2*((-B/2+g):g:(B/2-g))
        y=exp(1j*2*pi*f*t);
        s=s+y;
    end
    
    %received signal
    r=x+s;  
    X=fftshift(fft(x,N_zp));
    S=fftshift(fft(s,N_zp));
    R=fftshift(fft(r,N_zp));
    
    omeg=(-(N_zp/2):(N_zp/2)-1)*2*pi/(N_zp);   %DT Frequency
    f2=f_l+(B/2)*((omeg/pi)+1);
    
    %Windowing
    H=hanning(length(x)); %hanning
    hann=transpose(H);
    r2=r.*hann;
    R2=fftshift(fft(r2,N_zp));
    f3=(-N/2:N/2-1)*2*pi/N;
    
    %plots
    figure(1);
    subplot(3,1,1)
    plot(f2,20*log10(abs(X)));
    xlabel('RF frequency (Hz)');
    ylabel('Magnitude(decibel)');
    title('Equivalent Low pass signal');
    axis([10^9 1.1*10^9 -100 0]);
    grid
    subplot(3,1,2);
    plot(f2,20*log10(abs(S)));
    xlabel('RF frequency (Hz)');
    ylabel('Magnitude(decibel)');
    title('Jamming signal');
    axis([10^9 1.1*10^9 -100 100]);
    grid
    subplot(3,1,3);
    plot(f2,20*log10(abs(R)));
    xlabel('RF frequency (Hz)');
    ylabel('Magnitude(decibel)');
    title('Received signal');
    axis([10^9 1.1*10^9 -100 100]);
    grid
    
    figure(2);
    subplot(2,1,1);
    plot(f3, hann,'g');
    title('Hanning Window');
    subplot(2,1,2);
    grid
    plot(f2, 20*log10(abs(R2)),'g');
    title('Detected RF signal');
    xlabel('RF frequency (Hz)');
    ylabel('Magnitude(decibel)');
    grid
    
%Scenario #1 (Sinusoidal Jamming with addition of Additive White Gausian
%Noise) - Part B
    
    %desired signal
    x1=(10^-1.5)*exp(1j*2*pi*fm*t1);
    X1=fftshift(fft(x1,4*N1));
    
    %jamming signal (gaussian noise)
    sx=size(x1);
    gn=(1/sqrt(2))*(randn(sx)+ 1j*randn(sx));
    GN=fftshift(fft(gn,4*N1));
    
    %windowing 
    w=kaiser(N1,4);
    win=transpose(w);
    W=fftshift(fft(win,2*N1));
   
    %received signal 
    z=(x1+gn).*win;
    
    %processing at different values of N
    Z=fftshift(fft(z(1:N1)));
    Z1=fftshift(fft(z(1:N1), 2*N1));
    Z2=fftshift(fft(z(1:N1), 4*N1));
    Z3=fftshift(fft(z(1:N1), 8*N1));
    Z4=fftshift(fft(z(1:N1), 16*N1));
    
    omega=(-N1/2:N1/2-1)*2*pi/(N1);
    f=f_l+(B/2)*((omega/pi)+1);
    omega1=(-(2*N1/2):(2*N1/2)-1)*2*pi/(2*N1);
    f1=f_l+(B/2)*((omega1/pi)+1);
    omega2=(-(N1*4/2):(N1/2*4)-1)*2*pi/(N1*4);
    f2=f_l+(B/2)*((omega2/pi)+1);
    omega3=(-(N1*8/2):(N1*8/2)-1)*2*pi/(N1*8);
    f3=f_l+(B/2)*((omega3/pi)+1);
    omega4=(-(N1*16/2):(N1*16/2)-1)*2*pi/(N1*16);
    f4=f_l+(B/2)*((omega4/pi)+1);
    
    %plots
    figure(3);
    subplot(2,2,1)
    plot(t1,20*log10(abs(x1)));
    xlabel('Time(sec)');
    ylabel('Magnitude (decibel)');
    title('Equivalent Low pass signal');
    grid
    subplot(2,2,2)
    plot(t1,20*log10(abs(gn)));
    xlabel('Time(sec)');
    ylabel('Magnitude (decibel)');
    title('Complex Noise');
    grid
    subplot(2,2,3)
    plot(t1,20*log10(abs(x1+gn)));
    xlabel('Time(sec)');
    ylabel('Magnitude (decibel)');
    title('Signal added with Noise');
    grid
    subplot(2,2,4)
    plot(t1,20*log10(abs(z)));
    xlabel('Time(sec)');
    ylabel('Magnitude (decibel)');
    title('Signal added with noise (after windowing)');
    grid
    
    figure(4);
    subplot(2,1,1);
    plot((-N1/2:N1/2-1),win,'g');
    title('Kaiser Window (time domain)');
    grid
    subplot(2,1,2);
    plot(omega1,10*log10(abs(W)),'g');
    title('Kaiser Window (frequency domain)');
    grid
    
    figure(5);
    subplot(2,2,1);
    plot(f1,10*log10(abs(Z1)));
    title('2X Zero padding');
    grid
    subplot(2,2,2);
    plot(f2,10*log10(abs(Z2)));
    title('4X Zero padding');
    grid
    subplot(2,2,3);
    plot(f3,10*log10(abs(Z3)));
    title('8X Zero padding');
    grid
    subplot(2,2,4);
    plot(f4,10*log10(abs(Z4)));
    title('16X Zero padding');
    grid
    
    figure(6);
    subplot(2,2,1);
    plot(f2,10*log10(abs(X1)));
    xlabel('RF frequency(Hz)');
    ylabel('Magnitude(decibel)');
    title('Signal peak');
    axis([10^9 1.1*10^9 -150 60]);
    grid
    subplot(2,2,2);
    plot(f2,10*log10(abs(GN)));
    xlabel('RF frequency(Hz)');
    ylabel('Magnitude (decibel)');
    title('AWGN (frequency domain)');
    grid
    subplot(2,2,3);
    plot(f2,10*log10(abs(X1+GN)));
    xlabel('RF range(Hz)');
    ylabel('Magnitude (decibel)');
    title('Signal added with AWG noise');
    grid
    subplot(2,2,4);
    plot(f,10*log10(abs(Z)));
    xlabel('RF range(Hz)');
    ylabel('Magnitude (decibel)');
    title('Signal added with AWG noise (after windowing)');
    grid