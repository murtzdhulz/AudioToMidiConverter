function monoSynthFM(M,Fs,b,fc,wavef1,wavef2,At,Av,Dt,Dv,St,Rt)
    y=0;
    x1=str2func(wavef1);
    x2=str2func(wavef2);
    for i=1:length(M(:,1))
        t=M(i,5):1/Fs:M(i,6);
        f= (2^((M(i,3)-69)/12))*(440);
        Ae=linspace(0,Av,At*length(t));
        De=linspace(Av,Dv,Dt*length(t));
        Se=linspace(Dv,Dv,St*length(t));
        Re=linspace(Dv,0,Rt*length(t));
        ADSR=[Ae De Se Re];
        ADSR=padarray(ADSR,[0 (length(t)-length(ADSR))],0,'post');
        a=M(i,4)/127;
        b=a*b;
        I_env = 5.*ADSR(1,:);
        y = [y (b*x2(2*fc*pi*t+I_env(1,:).*(a*x1(f*2*pi*t)))).*ADSR(1,:)];
    end
    figure(10);
    plot(ADSR);
    soundsc(y,Fs);
end