# Code for Figs 3 and 4 from
# Dynamics and bifurcations of the adaptive exponential integrate-and-fire model,
# Touboul J and Brette R. Biol Cyber (2008).
# --------------- ------------------------------
# R. Brette (2008) brette@di.ens.fr
#
from brian import *

# FIGURE PARAMETERS
rc('lines',linewidth=2)
rc('font',size=16)
rc('xtick',labelsize=16)
rc('ytick',labelsize=16)
rc('legend',fontsize=16)
rc('axes',labelsize=16,titlesize=16)
w,h=rcParamsDefault['figure.figsize']

fontsize=16

# EQUATIONS
C=281*pF # Can be fixed
gL=30*nS
taum=C/gL
EL=-70.6*mV # Same as changing I
VT=-50.4*mV
DeltaT=3*mV
Vcut=VT+5*DeltaT
b=0*nA
Vr=-70.6*mV
tauw=taum
a=gL

eqs="""
dvm/dt=(gL*(EL-vm)+gL*DeltaT*exp((vm-VT)/DeltaT)+I-w)/C : volt
dw/dt=(a*(vm-EL)-w)/tauw : amp
I : amp
"""

# Backward equations
eqs_rev="""
dvm/dt=-(gL*(EL-vm)+gL*DeltaT*exp((vm-VT)/DeltaT)+I-w)/C : volt
dw/dt=-(a*(vm-EL)-w)/tauw : amp
I : amp
"""

vnullcline=lambda v:gL*(EL-v)+gL*DeltaT*exp((v-VT)/DeltaT)+I
wnullcline=lambda v:a*(v-EL)
def vplus():
    return optimize.fsolve(lambda v:gL*(EL-v*.001)+gL*DeltaT*exp((v*.001-VT)/DeltaT)+I-a*(v*.001-EL),Vcut/mV,xtol=0.00001)*mV
I_SN=lambda :(gL+a)*(VT-EL-DeltaT+DeltaT*log(1+a/gL))
I_Hopf=lambda :(gL+a)*(VT-EL-DeltaT+DeltaT*log(1+taum/tauw))+DeltaT*gL*(a/gL-taum/tauw)

def myreset(P,spikes):
    P.vm_[spikes]=Vr
    P.w_[spikes]+=b

def pics():
    # Pics: integrator/resonator/spike/limit cycle
    figure()
    x=linspace(0,1,200)
    subplot(221)
    plot(x,cos(15*x)*exp(-2*x)) # resonator
    subplot(222)
    plot(x,exp(-4*x)) # integrator
    subplot(223)
    plot(cos(2*pi*x),sin(2*pi*x)) # limit cycle
    subplot(224)
    plot(x,(2+sin(2*pi*x**.6))**4) # spike
    savefig('pics.eps')

def bigdiagram():
    global I,a,tauw
    myclock=Clock(dt=0.01*ms)
    # Class diagram of the B-G model
    # Subthreshold dynamics
    x=linspace(0.01,4,100) # taum/tauw
    y=linspace(0.01,4,100) # a/gl
    tau1=.2
    tau2=2
    
    figure(figsize=(w,2*h))
    subplot(511)
    plot(x,y,'k')
    plot(x,.25*x*(1-1/x)**2,'k')
    plot(0*x+tau1,y,'b--')
    plot(0*x+tau2,y,'b--')
    xlabel('$\\tau_m/\\tau_w$')
    ylabel('$a/g_L$')
    ylim(0,1)
    
    Vt=15. # mV
    deltat=3. # mV
    x=linspace(0.01,3,500) # a/gl
    # I/gl
    Irh1=(1+x)*(Vt-deltat+deltat*log(1+x))
    
    n=523
    for tau in [tau1,tau2]:
        subplot(n)
        n+=1
        x2=linspace(tau,3,500)
        Irh2=(1+x2)*(Vt-deltat+deltat*log(1+tau))+deltat*(x2-tau)
        Icycle=Irh2[0]-12./25.*deltat/(tau**2+tau)*(x2-x2[0])**2
        um=1-tau-2*(x*tau)**.5
        Im=(1+x)*(Vt+deltat*(log(um)))-deltat*um
        xm=x[Im!=NaN]
        Im=Im[Im!=NaN]
        up=1-tau+2*(x*tau)**.5
        Ip=(1+x)*(Vt+deltat*(log(up)))-deltat*up
        xp=x[(Ip!=NaN) & (Ip<Irh1)]
        Ip=Ip[(Ip!=NaN) & (Ip<Irh1)]
        plot(x,Irh1,'b')
        plot(x2,Irh2,'b')
        plot(x2,Icycle,'b--')
        plot(xm,Im,'r')
        plot(xp,Ip,'g')
        ylabel('$I/g_L$')
        xlabel('$a/g_L$')
        
    # Resonator
    subplot(525)
    tauw=taum
    a=10*gL
    I=0*mV*gL
    v=linspace(EL-5*mV,EL+5*mV,200)
    plot(v/mV,wnullcline(v)/nA)
    plot(v/mV,vnullcline(v)/nA)
    xlabel('V (mV)')
    ylabel('w (nA)')
    neuron=NeuronGroup(1,model=eqs)
    neuron.I=I
    neuron.vm=EL+I/(gL+a)
    neuron.w=a*(neuron.vm-EL)
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    net=Network(neuron,mv,mw)
    net.run(10*ms)
    neuron.vm+=5*mV
    net.run(50*ms)
    plot(mv[0]/mV,mw[0]/nA)
    ylim(-1,1)
    subplot(526)
    plot(mv.times/ms,mv[0]/mV)
    xlabel('Time (ms)')
    ylabel('V (mV)')
    xlim(0,60)
    myclock.reinit()
    
    # Integrator
    subplot(527)
    tauw=taum/2
    a=.1*gL
    I=0*mV*gL
    v=linspace(EL-5*mV,EL+5*mV,200)
    plot(v/mV,wnullcline(v)/nA)
    plot(v/mV,vnullcline(v)/nA)
    xlabel('V (mV)')
    ylabel('w (nA)')
    neuron=NeuronGroup(1,model=eqs)
    neuron.I=I
    neuron.vm=EL+I/(gL+a)
    neuron.w=a*(neuron.vm-EL)
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    net=Network(neuron,mv,mw)
    net.run(10*ms)
    neuron.vm+=5*mV
    net.run(50*ms)
    plot(mv[0]/mV,mw[0]/nA)
    ylim(-.02,.02)
    subplot(528)
    plot(mv.times/ms,mv[0]/mV)
    xlabel('Time (ms)')
    ylabel('V (mV)')
    xlim(0,60)
    myclock.reinit()    

    savefig('diagram.eps')

def IFcurve():
    N=2000
    duration=5*second
    neuron=NeuronGroup(2*N,model=eqs+'a:siemens',threshold=Vcut,reset=myreset)
    neuron.a[:N]=.1*gL
    neuron.a[N:]=3*gL
    neuron.vm=EL
    neuron.w=0*amp
    neuron.I[:N]=linspace(0.4*nA,1.2*nA,N)
    neuron.I[N:]=linspace(1.8*nA,2.6*nA,N)
    counter=SpikeCounter(neuron)
    run(duration)
    
    figure()
    # Type I
    subplot(121)
    plot(neuron.I[:N]/nA,counter.count[:N]/duration)
    # Type II
    subplot(122)
    plot(neuron.I[N:]/nA,counter.count[N:]/duration)

def nullclines():
    global I
    figure()
    v=linspace(EL+5*mV,VT+10*mV,200)
    subplot(121)
    I=I_SN()-.2*nA
    plot(v/mV,wnullcline(v)/nA)
    plot(v/mV,vnullcline(v)/nA)
    xlabel('v (mV)')
    ylabel('w (nA)')
    subplot(122)
    I=I_SN()+.2*nA
    plot(v/mV,wnullcline(v)/nA)
    plot(v/mV,vnullcline(v)/nA)
    xlabel('v (mV)')
    ylabel('w (nA)')

def typeI():
    global I,a,tauw
    a=.2*gL
    tauw=taum/3.
    myclock=Clock(dt=.01*ms)
    figure()

    # Nullclines
    v0=-58.2*mV
    w0=0.5*nA

    v=linspace(EL+2*mV,VT+10*mV,200)
    subplot(221)
    I=I_SN()-.2*nA
    plot(v/mV,wnullcline(v)/nA)
    plot(v/mV,vnullcline(v)/nA)
    xlabel('v (mV)')
    ylabel('w (nA)')
    # Trajectory 1
    neuron=NeuronGroup(1,model=eqs)
    neuron.vm=v0
    neuron.w=w0
    neuron.I=I
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    Network(neuron,mv,mw).run(100*ms)
    plot(mv[0]/mV,mw[0]/nA)
    ylim(-.2,2.5)
    myclock.reinit()
    
    subplot(222)
    I=I_SN()+.2*nA
    plot(v/mV,wnullcline(v)/nA)
    plot(v/mV,vnullcline(v)/nA)
    xlabel('v (mV)')
    ylabel('w (nA)')
    # Trajectory 2
    neuron=NeuronGroup(1,model=eqs)
    neuron.vm=v0
    neuron.w=w0
    neuron.I=I
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    Network(neuron,mv,mw).run(50*ms)
    ind=mv[0]<max(v)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    ylim(-.2,2.5)
    myclock.reinit()

    # I-F curve
    subplot(212)
    N=2000
    duration=20*second
    neuron=NeuronGroup(N,model=eqs,threshold=Vcut,reset=myreset)
    neuron.vm=EL
    neuron.w=0*amp
    neuron.I=linspace(0.4*nA,1.*nA,N)
    counter=SpikeCounter(neuron)
    Network(neuron,counter).run(duration)
    plot(neuron.I/nA,counter.count/duration)
    xlabel('I (nA)')
    ylabel('F (Hz)')
    savefig('typeI.eps')
    
def typeII():
    global I,a,tauw
    a=3*gL
    tauw=2*taum
    myclock=Clock(dt=.01*ms)
    figure()

    # Nullclines
    v0=-48.5*mV
    w0=2*nA

    v=linspace(EL+2*mV,VT+15*mV,200)
    subplot(221)
    I=I_Hopf()-.2*nA
    plot(v/mV,wnullcline(v)/nA)
    plot(v/mV,vnullcline(v)/nA)
    xlabel('v (mV)')
    ylabel('w (nA)')
    # Trajectory 1
    neuron=NeuronGroup(1,model=eqs)
    neuron.vm=v0
    neuron.w=w0
    neuron.I=I
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    Network(neuron,mv,mw).run(200*ms)
    ind=mv[0]<max(v)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    ylim(1,3)
    myclock.reinit()
    
    subplot(222)
    I=I_Hopf()+.05*nA
    plot(v/mV,wnullcline(v)/nA)
    plot(v/mV,vnullcline(v)/nA)
    xlabel('v (mV)')
    ylabel('w (nA)')
    # Trajectory 2
    neuron=NeuronGroup(1,model=eqs)
    neuron.vm=v0
    neuron.w=w0
    neuron.I=I
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    Network(neuron,mv,mw).run(250*ms)
    ind=mv[0]<max(v)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    ylim(1,3)
    myclock.reinit()

    # I-F curve
    subplot(212)
    N=2000
    duration=20*second
    neuron=NeuronGroup(N,model=eqs,threshold=Vcut,reset=myreset)
    neuron.vm=EL
    neuron.w=0*amp
    neuron.I=linspace(1.6*nA,2.6*nA,N)
    counter=SpikeCounter(neuron)
    Network(neuron,counter).run(duration)
    plot(neuron.I/nA,counter.count/duration)
    xlabel('I (nA)')
    ylabel('F (Hz)')
    savefig('typeII.eps')
    
def attractor():
    # attractor and rebound
    global I,a,tauw
    
    ### Resonator type II
    figure(figsize=(w,2*h))
    subplot(421)
    a=3*gL
    tauw=2*taum
    I=I_Hopf()-.1*nA
    myclock=Clock(dt=.05*ms)
    v=linspace(EL+10*mV,VT+10*mV,200)
    # nullclines
    plot(v/mV,wnullcline(v)/nA)
    plot(v/mV,vnullcline(v)/nA)

    # Limit circle
    neuron=NeuronGroup(1,model=eqs_rev)
    neuron.vm=EL+I/(gL+a)
    neuron.w=a*(neuron.vm-EL)
    neuron.I=I
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    Network(neuron,mv,mw).run(500*ms)
    ind=(mv[0]<max(v)) & (mv[0]>min(v)) & (mw[0]<5*nA) & (mw.times>200*ms)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    myclock.reinit()
    
    # Resting point
    v0=EL
    w0=0*nA
    def myreset(P,spikes):
        P.vm_[spikes]=v0
        P.w_[spikes]=w0
    neuron=NeuronGroup(1,model=eqs,threshold=Vcut,reset=myreset)
    neuron.vm=EL+I/(gL+a)
    neuron.w=a*(neuron.vm-EL)
    neuron.I=I
    Network(neuron).run(100*ms)
    v0=neuron.vm[0]
    w0=neuron.w[0]
    myclock.reinit()
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    net=Network(neuron,mv,mw)
    net.run(30*ms)
    neuron.vm-=6*mV
    net.run(50*ms)
    neuron.I-=.4*nA
    net.run(150*ms)
    neuron.I+=.4*nA
    net.run(50*ms)
    ind=(mv[0]<0*mV) & (mv[0]>-100*mV) & (mw[0]<5*nA)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    xlim(min(v)/mV,max(v)/mV)
    
    subplot(422)
    plot(mv.times/ms,mv[0]/mV)
    myclock.reinit()
    
    # --------------------
    # Resonator type I
    # --------------------
    subplot(423)
    a=10*gL
    tauw=taum/12
    I=0*nA
    v=linspace(EL-300*mV,VT+25*mV,200)
    # nullclines
    plot(v/mV,wnullcline(v)/nA)
    plot(v/mV,vnullcline(v)/nA)

    # Separatrix (approximation)
    neuron=NeuronGroup(1,model=eqs_rev)
    neuron.vm=vplus()
    neuron.w=a*(neuron.vm-EL)+.02*nA
    neuron.I=I
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    Network(neuron,mv,mw).run(200*ms)
    ind=(mv[0]<max(v)) & (mv[0]>min(v)) & (mw[0]<100*nA)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    myclock.reinit()
    neuron=NeuronGroup(1,model=eqs_rev)
    neuron.vm=vplus()
    neuron.w=a*(neuron.vm-EL)-.02*nA
    neuron.I=I
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    Network(neuron,mv,mw).run(200*ms)
    ind=(mv[0]<max(v)) & (mv[0]>min(v)) & (mw[0]<100*nA)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    myclock.reinit()
    
    # Resting point
    v0=EL
    w0=0*nA
    def myreset(P,spikes):
        P.vm_[spikes]=v0
        P.w_[spikes]=w0
    neuron=NeuronGroup(1,model=eqs,threshold=Vcut,reset=myreset)
    neuron.vm=EL+I/(gL+a)
    neuron.w=a*(neuron.vm-EL)
    neuron.I=I
    Network(neuron).run(100*ms)
    v0=neuron.vm[0]
    w0=neuron.w[0]
    myclock.reinit()
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    net=Network(neuron,mv,mw)
    net.run(30*ms)
    neuron.vm-=6*mV
    net.run(50*ms)
    neuron.I-=.4*nA
    net.run(150*ms)
    neuron.I+=.4*nA
    net.run(50*ms)
    ind=(mv[0]<0*mV) & (mv[0]>-100*mV) & (mw[0]<5*nA)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    xlim(min(v)/mV,max(v)/mV)
    ylim(-100,100)
    
    subplot(424)
    plot(mv.times/ms,mv[0]/mV)
    myclock.reinit()
    
    # --------------------
    # Integrator
    # --------------------
    subplot(425)
    a=.2*gL
    tauw=taum/3
    I=0*nA
    v=linspace(-100*mV,VT+15*mV,200)
    # nullclines
    plot(v/mV,wnullcline(v)/nA)
    plot(v/mV,vnullcline(v)/nA)

    # Separatrix (approximation)
    neuron=NeuronGroup(1,model=eqs_rev)
    neuron.vm=vplus()
    neuron.w=a*(neuron.vm-EL)+.02*nA
    neuron.I=I
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    Network(neuron,mv,mw).run(200*ms)
    ind=(mv[0]<max(v)) & (mv[0]>min(v)) & (mw[0]<100*nA)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    myclock.reinit()
    neuron=NeuronGroup(1,model=eqs_rev)
    neuron.vm=vplus()
    neuron.w=a*(neuron.vm-EL)-.02*nA
    neuron.I=I
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    Network(neuron,mv,mw).run(200*ms)
    ind=(mv[0]<max(v)) & (mv[0]>min(v)) & (mw[0]<100*nA)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    myclock.reinit()
    
    # Resting point
    v0=EL
    w0=0*nA
    def myreset(P,spikes):
        P.vm_[spikes]=v0
        P.w_[spikes]=w0
    neuron=NeuronGroup(1,model=eqs,threshold=Vcut,reset=myreset)
    neuron.vm=EL+I/(gL+a)
    neuron.w=a*(neuron.vm-EL)
    neuron.I=I
    Network(neuron).run(100*ms)
    v0=neuron.vm[0]
    w0=neuron.w[0]
    myclock.reinit()
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    net=Network(neuron,mv,mw)
    net.run(30*ms)
    neuron.vm-=6*mV
    net.run(50*ms)
    neuron.I-=.4*nA
    net.run(150*ms)
    neuron.I+=.4*nA
    net.run(50*ms)
    ind=(mv[0]<0*mV) & (mv[0]>-100*mV) & (mw[0]<5*nA)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    xlim(min(v)/mV,max(v)/mV)
    ylim(-7,10)
    
    subplot(426)
    plot(mv.times/ms,mv[0]/mV)
    myclock.reinit()
    
    # ------------------------
    # Semi-oscillator, type II
    # ------------------------
    subplot(427)
    a=gL
    tauw=taum*10
    I=0*nA
    v=linspace(-150*mV,VT+15*mV,200)
    # nullclines
    plot(v/mV,wnullcline(v)/nA)
    plot(v/mV,vnullcline(v)/nA)
    xlabel('v (mV)')
    ylabel('w (nA)')

    # Separatrix (approximation)
    neuron=NeuronGroup(1,model=eqs_rev)
    neuron.vm=vplus()
    neuron.w=a*(neuron.vm-EL)+.02*nA
    neuron.I=I
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    Network(neuron,mv,mw).run(1000*ms)
    ind=(mv[0]<max(v)) & (mv[0]>min(v)) & (mw[0]<100*nA)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    myclock.reinit()
    neuron=NeuronGroup(1,model=eqs_rev)
    neuron.vm=vplus()
    neuron.w=a*(neuron.vm-EL)-.02*nA
    neuron.I=I
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    Network(neuron,mv,mw).run(1000*ms)
    ind=(mv[0]<max(v)) & (mv[0]>min(v)) & (mw[0]<100*nA)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    myclock.reinit()
    
    # Resting point
    v0=EL
    w0=0*nA
    def myreset(P,spikes):
        P.vm_[spikes]=v0
        P.w_[spikes]=w0
    neuron=NeuronGroup(1,model=eqs,threshold=Vcut+5*mV,reset=myreset)
    neuron.vm=EL+I/(gL+a)
    neuron.w=a*(neuron.vm-EL)
    neuron.I=I
    Network(neuron).run(100*ms)
    v0=neuron.vm[0]
    w0=neuron.w[0]
    myclock.reinit()
    mv=StateMonitor(neuron,'vm',record=0)
    mw=StateMonitor(neuron,'w',record=0)
    net=Network(neuron,mv,mw)
    net.run(30*ms)
    neuron.vm-=6*mV
    net.run(50*ms)
    neuron.I-=2.5*nA
    net.run(150*ms)
    neuron.I+=2.5*nA
    net.run(50*ms)
    ind=(mv[0]<0*mV) & (mw[0]<5*nA)
    plot(mv[0][ind]/mV,mw[0][ind]/nA)
    xlim(min(v)/mV,max(v)/mV)
    ylim(-3,6)
    
    subplot(428)
    plot(mv.times/ms,mv[0]/mV)
    xlabel('Time (ms)')
    ylabel('V (mV)')
    myclock.reinit()
    
    savefig('rebound.eps')
    
#pics()
#bigdiagram()
#IFcurve()
#nullclines()
#typeI()
#typeII()
#attractor()
show()
