# Response of Brette-Gerstner model with different parameter values
#
# See the following paper:
# Dynamics and bifurcations of the adaptive exponential integrate-and-fire model,
# Touboul J and Brette R. Biol Cyber (2008).
# --------------- ------------------------------
# R. Brette (2008) brette@di.ens.fr
#

from brian import *
import time

cl=Clock(dt=0.01*ms)

N=1
C=281*pF # Can be fixed
gL=30*nS
taum=C/gL
EL=-70.6*mV # Same as changing I
VT=-50.4*mV
DeltaT=2*mV
Vcut=VT+5*DeltaT

I0=0*nA
dv0=0*mV # voltage transient
T0=50*ms
T1=500*ms
T2=100*ms
name='rebound_burst'
if name=='latency':
    tauw=150*ms
    a=2*C/tauw # type II
    b=0*nA
    Vr=-70.6*mV
    #I=0*nA
    EL=-60*mV
    I_john=(1+a/gL)*log(1+taum/tauw)-(1+taum/tauw)
    I0=gL*DeltaT*I_john+(VT-EL)*(gL+a)-0.03*nA
    I=I0
    T0=100*ms
    T1=200*ms
    T2=20*ms
    dv0=2.5*mV
elif name=='rebound_spike':
    tauw=150*ms
    #T0=
    a=200*nS
    b=0.1*nA
    I0=0*nA
    I=-0.5*nA
    EL=-60*mV
    T1=50*ms
    VT=-54*mV
    Vr=EL
elif name=='rebound_burst':
    tauw=150*ms
    #T0=
    a=200*nS
    b=0.1*nA
    I0=0*nA
    I=-0.5*nA
    EL=-60*mV
    T1=50*ms
    VT=-54*mV
    Vr=VT+3*mV
elif name=='regular':
    tauw=144*ms
    a=4*nS
    b=0.0805*nA
    Vr=-70.6*mV
    I=1*nA
elif name=='phasic':
    tauw=144*ms
    a=4*nS
    b=0.0805*nA
    Vr=-70.6*mV
    I=.6*nA
elif name=='on_off':
    tauw=10*ms
    a=800*nS
    b=10*nA
    Vr=-70.6*mV
    I=11*nA
elif name=='oscillator??': # numerical integration problem
    tauw=14*ms
    a=8000*nS
    b=0.0805*nA
    Vr=-70.6*mV
    I=1.3*nA
elif name=='mixed':
    tauw=144*ms
    a=4*nS
    b=0.15*nA
    Vr=VT+2*mV
    I=1*nA
elif name=='fast':
    # Type II: a > C/tauw
    tauw=144*ms
    a=2*C/tauw
    b=0*nA
    Vr=-70.6*mV
    I_john=(1+a/gL)*log(1+taum/tauw)-(1+taum/tauw)
    I=gL*DeltaT*I_john+(VT-EL)*(gL+a)+0.01*nA
    I0=I-0.1*nA
    T0=200*ms
elif name=='bursting_tonic':
    tauw=20*ms
    a=4*nS
    b=0.5*nA
    Vr=VT+5*mV
    I=.8*nA
elif name=='bursting_phasic':
    tauw=144*ms
    a=4*nS
    b=0.1*nA
    Vr=VT+4*mV
    I=.6*nA
elif name=='bursting_phasic2': # big burst, then small bursts
    tauw=144*ms
    a=4*nS
    b=0.15*nA
    Vr=VT+2*mV
    I=1.5*nA

eqs=Equations("""
dvm/dt=(gL*(EL-vm)+gL*DeltaT*exp((vm-VT)/DeltaT)+I-w)/C : volt
dw/dt=(a*(vm-EL)-w)/tauw : amp
I : amp
""")

def myreset(P,spikes):
    P.vm_[spikes]=Vr
    P.w_[spikes]+=b

neuron=NeuronGroup(N,model=eqs,threshold=Vcut,reset=myreset,freeze=True)
# NB: freezing removes units from equations

# Solve for resting potential
dvm=eqs._function['vm']
dw=eqs._function['w']
v0,w0=optimize.fsolve(lambda x:array([dvm(x[0],x[1],0.),dw(x[0],x[1])]),[EL,0*nA])
#print v0*volt,w0*amp
neuron.vm=v0*volt
neuron.w=w0*amp
#neuron.w=8*nA

mon_v=StateMonitor(neuron,'vm',record=0)
mon_I=StateMonitor(neuron,'I',record=0)
mon_w=StateMonitor(neuron,'w',record=0)
spikes=SpikeMonitor(neuron)

neuron.I=I0
run(T0)
neuron.I=I
neuron.vm+=dv0
run(T1)
neuron.I=0*nA
run(T2)

# Plots
# Trace
vm=mon_v[0]
for _,t in spikes.spikes:
    i=int(t/cl.dt)
    vm[i]=20*mV
subplot(221)
plot(mon_v.times/ms,vm/mV)

# w
subplot(222)
plot(mon_w.times/ms,mon_w[0]/nA)

# I
subplot(223)
plot(mon_I.times/ms,mon_I[0]/nA)

subplot(224)
# Phase plot
vm=mon_v[0]
last_i=-1
for _,t in spikes.spikes:
    i=int(t/cl.dt)
    vm[i]=Vcut
    plot(vm[last_i+1:i+1]/mV,mon_w[0][last_i+1:i+1]/nA,'b')
    last_i=i
plot(vm[last_i+1:]/mV,mon_w[0][last_i+1:]/nA,'b')
# Spikes
for _,t in spikes.spikes:
    i=int(t/cl.dt)
    plot([vm[i]/mV,vm[i+1]/mV],[mon_w[0][i]/nA,mon_w[0][i+1]/nA],'b--')
wrange=array([min(mon_w[0]/nA),max(mon_w[0]/nA)])
plot((Vr/mV)*array([1,1]),wrange,'k') # rather from axes handle
# Nullclines
v=linspace(min(vm),max(vm),100)
plot(v/mV,(gL*(EL-v)+gL*DeltaT*exp((v-VT)/DeltaT)+I)/nA,'k') # v-nullcline
plot(v/mV,(a*(v-EL))/nA,'k') # w-nullcline
xlim(min(v/mV),max(v/mV))
ylim(wrange[0],wrange[1])

show()
