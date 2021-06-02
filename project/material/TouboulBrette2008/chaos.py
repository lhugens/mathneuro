# Code for Fig 7 from
# Dynamics and bifurcations of the adaptive exponential integrate-and-fire model,
# Touboul J and Brette R. Biol Cyber (2008).
# --------------- ------------------------------
# R. Brette (2008) brette@di.ens.fr
#
# Chaos in Brette-Gerstner model
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

C=281*pF
gL=30*nS
EL=-70.6*mV
VT=-50.4*mV
DeltaT=2*mV
tauw=40*ms #20*ms
a=4*nS
b=0.08*nA
I=.8*nA
Vcut=VT+5*DeltaT
N=500

eqs=Equations("""
dvm/dt=(gL*(EL-vm)+gL*DeltaT*exp((vm-VT)/DeltaT)+I-w)/C : volt
dw/dt=(a*(vm-EL)-w)/tauw : amp
Vr:volt
""")

l=[]
def myreset(P,spikes):
    P.vm_[spikes]=P.Vr_[spikes]
    P.w_[spikes]+=b
    if cl.t>3000*ms:
        l.extend(zip(P.Vr[spikes],P.w[spikes]))

def vminus():
    return optimize.fsolve(lambda v:gL*(EL-v*.001)+gL*DeltaT*exp((v*.001-VT)/DeltaT)-a*(v*.001-EL),EL/mV,xtol=0.00001)*mV

neuron=NeuronGroup(N,model=eqs,threshold=Vcut,reset=myreset)
neuron.vm=vminus()
neuron.w=a*(neuron.vm-EL)
#neuron.Vr=linspace(-49*mV,-46*mV,N)
neuron.Vr=linspace(-48.3*mV,-47.7*mV,N)

mon_v=StateMonitor(neuron,'vm',record=0)
mon_w=StateMonitor(neuron,'w',record=0)
spikes=SpikeMonitor(neuron)

run(5000*ms)

figure()
x,y=zip(*l)
plot(qarray(x)/mV,qarray(y)/nA,'.k')
xlabel('Vr (mV)')
ylabel('w (nA)')
#xticks([-49,-48,-47,-46])
#yticks([0.2,0.4,0.6])
#savefig('periodadding.eps')
savefig('periodaddingzoom.eps')
show()
