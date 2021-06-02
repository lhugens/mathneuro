from brian import *
from brian.library.IF import *
import time

cl=Clock(dt=0.1*ms)

N=1000
Vr=-70.6*mV+20*mV
C=281*pF
gL=30*nS
EL=-70.6*mV
VT=-50.4*mV
DeltaT=2*mV
tauw=.1*144*ms
a=4*nS
b=0.0805*nA

eqs="""
dvm/dt=(gL*(EL-vm)+gL*DeltaT*exp((vm-VT)/DeltaT)+I-w)/C : volt
dw/dt=(a*(vm-EL)-w)/tauw : amp
I : amp
w0 : amp
"""

# Computes phi(w)
wmin=-1*nA
wmax=3*nA
def myreset(P,spikes):
    P.vm_[spikes]=Vr
    P.w_[spikes]=rand(len(spikes))*(wmax-wmin)+wmin
    P.w0_[spikes]=P.w_[spikes]

neuron=NeuronGroup(N,model=eqs,threshold=-20*mV,reset=myreset)
neuron.vm=Vr
neuron.w=rand(N)*(wmax-wmin)+wmin
neuron.w0=neuron.w
neuron.I=1*nA

w=[]
phiw=[]
def record_spikes(spikes):
    w.extend(neuron.w0[spikes])
    phiw.extend(neuron.w[spikes])

mon=SpikeMonitor(neuron,function=record_spikes)

start_time=time.time()
run(500*ms)
print len(phiw)
print time.time()-start_time

phiw=qarray(phiw)+0.0805*nA
plot(qarray(w)/nA,phiw/nA,'.')
plot([wmin/nA,wmax/nA],[wmin/nA,wmax/nA])
show()
