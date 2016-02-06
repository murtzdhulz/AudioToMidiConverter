function monoSynthAdd(M,Fs,b,c,wavef1,wavef2,wavef3,At,Av,Dt,Dv,St,Rt)
    y=0;
    x1=str2func(wavef1);
    x2=str2func(wavef2);
    x3=str2func(wavef3);
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
        c=a*c;
        y = [y (a*x1(f*2*pi*t)+b*x2(f*pi*t)+c*x3(f*4*pi*t)).*ADSR(1,:)];
    end
    figure(10);
    plot(ADSR);
    soundsc(y,Fs);
end