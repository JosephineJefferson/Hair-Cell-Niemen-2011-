"""
Based on Nieman 2011
Model of bullfrog saccular hair cell
Hodgkin-Huxely type
"""

using Unitful
using Plots
gr()

const k_int = 112.0u"mmol/L" #Intracellular potassium conc. (millimolar)
const k_ext = 2.0u"mmol/L" #Extracellular potassium conc. (millimolar)
const F = 96485.3329u"s*A/mol" #Faraday constant
const R = 8.314u"J/(mol*K)" #universal gas constant
const T = 295.15u"K" #Temperature in Kalvin
const dt= 1.0e-5u"s" #timestep - theirs was 1e-5

"""
    Structure modelling a bullfrog saccular hair cell
"""
struct Hair_Cell
    #state vector(s)
    x::Array{Any,1} #contains V, m_K1f, m_K1s, m_h, m_DRK, m_Ca
    BK::Array{Any,1} #contains Ca_conc, h_BKT
    PoMET::Array{Any,1} #contains open state prob.
    currents::Array{typeof(1.0u"nA"),1} #stores current current values
    #Parameters
    #
    #Reversal Potentials
    E_K::typeof(1.0u"mV")
    E_h::typeof(1.0u"mV") #for cation h-current
    E_Ca::typeof(1.0u"mV")
    E_L::typeof(1.0u"mV")
    E_MET::typeof(1.0u"mV")
    #Maximal Conductances
    g_K1::typeof(1.0u"nS") #g_K1 is one of the main control parameters
    g_h::typeof(1.0u"nS") #for cation h-current
    g_Ca::typeof(1.0u"nS")
    g_L::typeof(1.0u"nS")
    g_MET::typeof(1.0u"nS")
    #Max permeabilities of currents
    P_DRK::typeof(1.0u"L/s")#"L/s")
    P_BKS::typeof(1.0u"L/s") #calcium current
    P_BKT::typeof(1.0u"L/s") #calcium current
    #b is another key control parameter (dimensionless)
    #b controls the strength of the above currents (BKS and BKT)
    b::Float64
    #cell capacitance
    Cm::typeof(1.0u"pF")
end

"""
    Contructor. Note: gK1 and b are the main control parameters of the model.
"""
function Hair_Cell(V::typeof(1.0u"mV"),g_K1=20.0u"nS",b=0.5) # g and b from (a)
    #Reversal Potentials
    E_K=-95.0u"mV" # (a,c)
    E_h=-45.0u"mV" #for cation h-current (e)
    E_Ca=42.5u"mV" # (a,f)
    E_L=-40.0u"mV" # (f)
    E_MET=0.0u"mV" # (a)
    #Maximal Conductances
    g_K1=g_K1#g_K1 is one of the main control parameters of the model
    g_h=2.2u"nS" #for cation h-current (e)
    g_Ca=1.2u"nS" # (d,e)
    g_L=0.77u"nS" # (b)
    g_MET=0.65u"nS" # (a)
    #Max permeabilities of currents
    P_DRK=2.4e-14u"L/s" # (a)
    P_BKS=2e-13u"L/s" #calcium current (a,c)
    P_BKT=14e-13u"L/s" #calcium current (a,c)
    #b is another key control parameter (dimensionless)
    #b controls the strength of the above currents (BKS and BKT)
    b=b
    #cell capacitance
    Cm=20.0u"pF" # (a,c)
    #initial m values for state vector, set as 0 for now
    m_K1f=m_k1_inf(V)
    m_K1s=m_k1_inf(V)
    m_h=m_h_inf(V)
    m_DRK=m_DRK_inf(V)
    m_Ca=m_Ca_inf(V)
    #initial values for params of BK-current stuff (another state vector)
    Ca_conc= 0.0u"mmol/L"#find something good
    h_BKT=h_BKT_inf(V) #h for inactivation
    C0=0.6
    C1=0.1
    C2=0.1
    O2=0.1
    O3=0.1

    PoMet=0.15
    return Hair_Cell([V,m_K1f, m_K1s, m_h, m_DRK, m_Ca], [Ca_conc, h_BKT,C0,C1,C2,O2,O3],[PoMet],[0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA"],
    E_K,E_h,E_Ca,E_L,E_MET,g_K1,g_h,g_Ca,g_L,g_MET,P_DRK,P_BKS,P_BKT,b,Cm)
end

function update(cell::Hair_Cell,input=0.0u"nA")
    V=cell.x[1]
    thing=(V*F)/(R*T)
    pow=uconvert(Unitful.NoUnits,thing)
    IK1 = cell.g_K1*(V-cell.E_K)*(0.7*mk1f(V,cell.x[2])+0.3*mk1s(V,cell.x[3]))
    cell.currents[1]=IK1
    Ih = cell.g_h*(V-cell.E_h)*(3*mh(V,cell.x[4])^2*(1-mh(V,cell.x[4]))+mh(V,cell.x[4])^3)
    cell.currents[2]=Ih
    IDRK = uconvert(u"nA",cell.P_DRK*((V*F^2)/(R*T))*((k_int-k_ext*exp(-pow))/(1-exp(-pow)))*(mDRK(V,cell.x[5])^2))
    cell.currents[3]=IDRK
    ICa = cell.g_Ca*(V-cell.E_Ca)*mCa(V,cell.x[6])^3
    cell.currents[4]=ICa
    IBKS = cell.b*cell.P_BKS*((V*F^2)/(R*T))*((k_int-k_ext*exp(-pow))/(1-exp(-pow)))*(O_2(V,cell.BK[5],cell.BK[6],cell.BK[7])+O_3(V,cell.BK[6],cell.BK[7]))
    cell.currents[5]=IBKS
    IBKT = cell.b*cell.P_BKT*((V*F^2)/(R*T))*((k_int-k_ext*exp(-pow))/(1-exp(-pow)))*(O_2(V,cell.BK[5],cell.BK[6],cell.BK[7])+O_3(V,cell.BK[6],cell.BK[7]))*hBKT(V,cell.BK[2])
    cell.currents[6]=IBKT
    IL = cell.g_L*(V-cell.E_L)
    cell.currents[7]=IL
    IMET = cell.g_MET*cell.PoMET[1]*(V-cell.E_MET)
    cell.currents[8]=IMET
    ΔV=uconvert(u"mV",(-IK1-Ih-IDRK-ICa-IBKS-IBKT-IL-IMET+input)*dt/cell.Cm)

    #update things
    cell.x[1]= V + ΔV #voltage/membrane potential
    cell.x[2]= mk1f(V,cell.x[2]) #m_K1f
    cell.x[3]= mk1s(V,cell.x[3]) #m_K1s
    cell.x[4]=   mh(V,cell.x[4]) #m_h
    cell.x[5]= mDRK(V,cell.x[5]) #m_DRK
    cell.x[6]=  mCa(V,cell.x[6]) #m_Ca
    #cell.BK[1]= Ca_conc(cell.BK[1],ICa) #Intracellular calcium concentration
    cell.BK[2]= hBKT(V,cell.BK[2]) #h_BKT
    cell.BK[3]= C_0(cell.BK[4],cell.BK[5],cell.BK[6],cell.BK[7]) #C0
    cell.BK[4]= C_1(V,cell.BK[3],cell.BK[4],cell.BK[5]) #C1
    cell.BK[5]= C_2(V,cell.BK[4],cell.BK[5],cell.BK[6]) #C2
    cell.BK[6]= O_2(V,cell.BK[5],cell.BK[6],cell.BK[7]) #02
    cell.BK[7]= O_3(V,cell.BK[6],cell.BK[7]) #03
end

#equations for K1 current
mk1f(V,mk1f)=mk1f+Δm_k1f(V,mk1f)
Δm_k1f(V,mk1f)=(m_k1_inf(V)-mk1f)*dt/τ_k1f(V)
τ_k1f(V)=0.7u"ms"*exp(-(V + 120u"mV")/43.8u"mV")+0.04u"ms"

mk1s(V,mk1s)=mk1s+Δm_k1s(V,mk1s)
Δm_k1s(V,mk1s)=(m_k1_inf(V)-mk1s)*dt/τ_k1s(V)
τ_k1s(V)=14.1u"ms"*exp(-(V + 120u"mV")/28u"mV")+0.04u"ms"

m_k1_inf(V)=(1+exp((V + 110u"mV")/11u"mV"))^-1

#equations for h cation current
mh(V,mh)=mh+Δm_h(V,mh)
Δm_h(V,mh)=(m_h_inf(V)-mh)*dt/τ_h(V)
m_h_inf(V)=(1+exp((V+87u"mV")/16.7u"mV"))^(-1)
τ_h(V)=63.7u"ms"+135.7u"ms"*exp(-((V+91.4u"mV")/21.2u"mV")^2)

#equations for DRK current
mDRK(V,mDRK)= mDRK+ΔmDRK(V,mDRK)
ΔmDRK(V, mDRK) = (m_DRK_inf(V)- mDRK)*(dt/τ_DRK(V))
m_DRK_inf(V) = (1+exp((V+48.3u"mV")/4.19u"mV"))^(-1/2) #CHANGED V1/2*********************

τ_DRK(V) = (α_DRK(V)+β_DRK(V))^(-1)
α_DRK(V) = (3.2u"ms"*exp(-V/20.9u"mV")+3.0u"ms")^(-1)
β_DRK(V) = (1467.0u"ms"*exp(V/5.96u"mV")+9.0u"ms")^(-1)

#equations for ICa current
mCa(V,mCa) = mCa+Δm_Ca(V,mCa)
Δm_Ca(V, mCa) = (m_Ca_inf(V)-mCa)*(dt/τ_Ca(V))
m_Ca_inf(V) = (1+exp(-(V+55u"mV")/12.2u"mV"))^(-1)
τ_Ca(V) = 0.046u"ms" + 0.325u"ms"*(exp(-((V+77u"mV")/51.67u"mV")^2))

#equations for BK currents:
#Params:
#Used directly in state-probability equations
const k_1=300.0u"1/s" #s^-1
const k_2=5000.0u"1/s" #s^-1
const k_3=1500.0u"1/s" #s^-1
const β_c=2500.0u"1/s"#1000.0u"1/s" #s^-1 #changed to OG value *******************************

#for auxiliary equations
const α_c_0=450.0u"1/s" #s^-1
const K1_0=6.0#u"μmol/L" #micromolar
const K2_0=45.0#u"μmol/L" #micromolar
const K3_0=20.0#u"μmol/L" #micromolar
const z = 2 #sign of charge of Ca2+
const δ1 = 0.2 #fraction of the electric field experienced by Ca2 at the 1st binding site
const δ2 = 0.0 #fraction of the electric field experienced by Ca2 at the 2nd binding site
const δ3 = 0.2 #fraction of the electric field experienced by Ca2 at the 3rd binding site
const V_A = 30.0u"mV" #potential used to express the voltage dependence of ac

α_c(V)=α_c_0*exp(-V/V_A)

K1(V)=K1_0*exp(δ1*z*uconvert(Unitful.NoUnits,(V*F)/(R*T)))
K2(V)=K2_0*exp(δ2*z*uconvert(Unitful.NoUnits,(V*F)/(R*T)))
K3(V)=K3_0*exp(δ3*z*uconvert(Unitful.NoUnits,(V*F)/(R*T)))

k1_Ca(V)=k_1/(K1(V))
k2_Ca(V)=k_2/(K2(V))
k3_Ca(V)=k_3/(K3(V))

C_0(C1,C2,O2,O3) = 1.0 - (C1 + C2 + O2 + O3)
C_1(V,C0,C1,C2)= C1 + ΔC_1(V,C0,C1,C2)
ΔC_1(V,C0,C1,C2)=uconvert(Unitful.NoUnits,(k1_Ca(V)*C0+k_2*C2-(k_1+k2_Ca(V))*C1)*dt)
C_2(V,C1,C2,O2)= C2 + ΔC_2(V,C1,C2,O2)
ΔC_2(V,C1,C2,O2)=uconvert(Unitful.NoUnits,(k2_Ca(V)*C1+α_c(V)*O2-(k_2+β_c)*C2)*dt)
O_2(V,C2,O2,O3)= O2 + ΔO_2(V,C2,O2,O3)
ΔO_2(V,C2,O2,O3)=uconvert(Unitful.NoUnits,(β_c*C2+k_3*O3-(α_c(V)+k3_Ca(V))*O2)*dt)
O_3(V,O2,O3)= O3+ΔO_3(V,O2,O3)
ΔO_3(V,O2,O3)=uconvert(Unitful.NoUnits,(k3_Ca(V)*O2-k_3*O3)*dt)

#Ca_conc(Ca,ICa)=Ca+ΔCa_conc(Ca,ICa)
#ΔCa_conc(Ca,ICa)=(-0.00061u"mol/(L*ms*nA)"*ICa - 2800u"1/ms"*Ca)*dt #this smells funky
#Ca_conc(Ca,ICa)=uconvert(u"mmol/L",Ca+ΔCa_conc(Ca,ICa))
#ΔCa_conc(Ca,ICa)=(-0.00061*ICa - 2800*Ca)*dt

#inactivation of BKT current
hBKT(V,hBKT) = hBKT+Δh_BKT(V,hBKT)
Δh_BKT(V,hBKT)= (h_BKT_inf(V)-hBKT)*(dt/τ_BKT(V))
h_BKT_inf(V) = (1+exp((V+61.6u"mV")/3.65u"mV"))^(-1)
τ_BKT(V) = 2.1u"ms"+9.4u"ms"*exp(-((V+66.9u"mV")/17.7u"mV")^2)

"""
    runs for a given time to let parameters settle into equilibrium state
"""
function burnin(cell, time=5.0u"s", displayV=true,displayEachC=false,displayallC=true,printMP=true)
    t=0.0u"s":dt:time
    voltBI=Array{typeof(1.0u"mV")}(undef, length(t))
    currents_p=Array{Any,2}(undef,length(t),8)
    for i in 1:length(t)
        update(cell)
        voltBI[i]=cell.x[1]
        currents_p[i,:].=cell.currents
    end
    if displayV
        display(plot(ustrip.(t),ustrip.(voltBI), title="Membrane Potential during Burn-in"))
    end
    if displayEachC
        plotEach(ustrip.(t),ustrip.(currents_p))
    end
    if displayallC
        plotTOGETHER(ustrip.(t),ustrip.(currents_p))
    end
    if printMP
        println(string("RMP after burn-in: ",cell.x[1]))
    end
    return cell
end

"""
    Called by burnin to plot a graph of each current over time
"""
function plotEach(t,currents)
    display(plot(t,currents[:,1], title="IK1"))
    display(plot(t,currents[:,2], title="Ih"))
    display(plot(t,currents[:,3], title="IDRK"))
    display(plot(t,currents[:,4], title="ICa"))
    display(plot(t,currents[:,5], title="IBKS"))
    display(plot(t,currents[:,6], title="IBKT"))
    display(plot(t,currents[:,7], title="IL"))
    display(plot(t,currents[:,8], title="IMET"))
end

"""
    Called by burnin to plot a graph of all currents over time
"""
function plotTOGETHER(t,currents)
    plot(t,currents[:,1], label="IK1")
    plot!(t,currents[:,2], label="Ih")
    plot!(t,currents[:,3], label="IDRK")
    plot!(t,currents[:,4], label="ICa")
    plot!(t,currents[:,5], label="IBKS")
    plot!(t,currents[:,6], label="IBKT")
    plot!(t,currents[:,7], label="IL")
    display(plot!(t,currents[:,8], label="IMET"))
end

"""
   pulse waveform same length as t - called by pulseReponse
"""
function pulse(time=0.5u"s", start=2.0e-1u"s", len=2.0e-1u"s", amplitude=0.005u"nA")
    t=0.0u"s":dt:time
    u = Array{typeof(1.0u"nA")}(undef, length(t))
    u.=0.0u"nA"
    u[findall( t-> (t>=start) & ( t<start+len), t)] .= amplitude

    return u
end

"""
    Models response of cell to a pulse of given amplitude and length
"""
function pulseResponse(cell, amplitude=0.005u"nA", time=0.5u"s", start=2.0e-1u"s", len=2.0e-1u"s", displayV=true,displayEachC=false,displayallC=true,printMP=false)
    t=0.0u"s":dt:time
    input = pulse(time,start,len,amplitude) #in nA
    voltages=Array{typeof(1.0u"mV")}(undef, length(t))
    currents_p=Array{Any,2}(undef,length(t),8)

    for i in 1:length(input)
        update(cell,input[i])
        voltages[i]=helga.x[1]
        currents_p[i,:].=cell.currents
    end
    if displayV
        plot(ustrip.(t),ustrip.(voltages), title="Membrane Potential during Pulse",ylabel="Membrane Potential (mV)")
        display(plot!(twinx(),ustrip.(t),ustrip.(input),ylims=(-1e-3,ustrip(amplitude)*3),ylabel="Pulse Amplitude (nA)",linecolor=:violet,label="pulse"))
    end
    if displayEachC
        plotEach(ustrip.(t),ustrip.(currents_p))
    end
    if displayallC
        plotTOGETHER(ustrip.(t),ustrip.(currents_p))
    end
    if printMP
        println(string("RMP after pulse: ",cell.x[1]))
    end
    return cell
end

const helga = Hair_Cell(-50.0u"mV")

burnin(helga)
pulseResponse(helga,0.01u"nA")
