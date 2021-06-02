# Code for Figs 5 and 6 from
# Dynamics and bifurcations of the adaptive exponential integrate-and-fire model,
# Touboul J and Brette R. Biol Cyber (2008).
# --------------- ------------------------------
# R. Brette (2008) brette@di.ens.fr
#

# Spike patterns
from brian import *

# FIGURE PARAMETERS
rc('lines',linewidth=2)
rc('font',size=16)
rc('xtick',labelsize=16)
rc('ytick',labelsize=16)
rc('legend',fontsize=16)
rc('axes',labelsize=16,titlesize=16)
width,height=rcParamsDefault['figure.figsize']

fontsize=16

cl=Clock(dt=0.01*ms)

# RS
duration=150*ms
Vr=-60*mV
C=281*pF
gL=30*nS
EL=-70.6*mV
VT=-50.4*mV
DeltaT=2*mV
tauw=144*ms*.5
a=4*nS
b=0.0805*nA*10
I=2*nA
Vcut=VT+5*DeltaT
name='phi'

# Burst 2
C=281*pF
gL=30*nS
EL=-70.6*mV
VT=-50.4*mV
DeltaT=2*mV
tauw=40*ms #20*ms
a=4*nS
b=0.08*nA
I=.8*nA
Vr=-48.5*mV
duration=1500*ms
name='burst2'

# Burst 3
Vr=-47.7*mV
name='burst3'

# Burst 4
Vr=-47.2*mV
name='burst4'

# Chaos
Vr=-48*mV
#I=1*nA
#duration=1000*ms
name='chaos'

def plotphi(wmin=-1*nA,wmax=3*nA,N=200):
    eqs="""
    dvm/dt=u*(gL*(EL-vm)+gL*DeltaT*exp((vm-VT)/DeltaT)+I-w)/C : volt
    dw/dt=u*(a*(vm-EL)-w)/tauw : amp
    u:1
    """
    # Computes phi(w)
    def myreset(P,spikes):
        P.u_[spikes]=0 # stop
    
    w0=linspace(wmin,wmax,N)
    neuron=NeuronGroup(N,model=eqs,threshold=Vcut,reset=myreset)
    neuron.vm=Vr
    neuron.u=1
    neuron.w=w0
    
    run(500*ms)
    phiw=neuron.w+b
    plot(w0/nA,phiw/nA,'k')
    plot(w0/nA,w0/nA,'k--')

eqs=Equations("""
dvm/dt=(gL*(EL-vm)+gL*DeltaT*exp((vm-VT)/DeltaT)+I-w)/C : volt
dw/dt=(a*(vm-EL)-w)/tauw : amp
""")

def myreset(P,spikes):
    P.vm_[spikes]=Vr
    P.w_[spikes]+=b

def vminus():
    return optimize.fsolve(lambda v:gL*(EL-v*.001)+gL*DeltaT*exp((v*.001-VT)/DeltaT)-a*(v*.001-EL),EL/mV,xtol=0.00001)*mV

neuron=NeuronGroup(1,model=eqs,threshold=Vcut,reset=myreset)
neuron.vm=vminus()
neuron.w=a*(neuron.vm-EL)

mon_v=StateMonitor(neuron,'vm',record=0,timestep=10)
mon_w=StateMonitor(neuron,'w',record=0,timestep=10)
spikes=SpikeMonitor(neuron)

run(duration)
# Plots
# Trace
vm=mon_v[0]
w=zeros(len(spikes.spikes)+1)*amp
w[0]=mon_w[0][0]
n=1
nstart=-1
for _,t in spikes.spikes:
    i=int(t/(cl.dt*10))
    vm[i]=20*mV
    #print mon_w[0]
    w[n]=mon_w[0][i+1]
    if (t>duration-300*ms) and (nstart==-1):
        nstart=n-1
    n+=1

figure(figsize=(width,.5*height))
subplot(221)
plot(mon_v.times/ms,vm/mV,'k')
ylabel('V (mV)')
xlim(duration/ms-300,duration/ms)
xticks([])
yticks([-70,-50,-20,20])

# w
subplot(223)
plot(mon_w.times/ms,mon_w[0]/nA,'k')
ylabel('w (nA)')
xticks([0,50,100,150,200,250,300])
#yticks([0,1,2])
xlim(duration/ms-300,duration/ms)

subplot(122)
plotphi(wmin=0.2*nA,wmax=0.4*nA)
w=w/nA
for i in range(nstart,len(w)-1):
    plot([w[i],w[i],w[i+1]],[w[i],w[i+1],w[i+1]],'k')
xticks([0.25,0.35])
yticks([0.25,0.35])
xlim(0.25,0.35)
ylim(0.25,0.35)
#ylabel('$w_{n+1} (nA)$')
#savefig(name+'.eps')

show()
