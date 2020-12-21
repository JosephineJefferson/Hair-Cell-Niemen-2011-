"""
Based on Nieman 2011
Model of bullfrog saccular hair cell
Hodgkin-Huxely type
"""

using Unitful
using Distributions
using Plots
using StatsBase
gr()

const k_int = 112.0u"mmol/L" #Intracellular potassium conc. (millimolar)
const k_ext = 2.0u"mmol/L" #Extracellular potassium conc. (millimolar)
const F = 96485.3329u"s*A/mol" #Faraday constant
const R = 8.314u"J/(mol*K)" #universal gas constant
const T = 295.15u"K" #Temperature in Kalvin
const kβ = 1.38e-23u"J/K" #Boltzmann constant
const zgate  = 40u"fN"# gatingforce 40 fN (Howard, Roberts & Hudspeth 1988)
const pr = 0.15   # resting/spontaneous open state probability
const max_deflect= 1000.0u"nm" #maximum kinocilium deflection
const min_deflect= -500.0u"nm" #minimum kinocilium deflection
const dt= 1.0e-5u"s" #timestep - theirs was 1e-5
const n_met_chan=50

"""
    Structure modelling a bullfrog saccular hair cell
"""
struct Hair_Cell
    #state vector(s)
    x::Array{Any,1} #contains V, m_K1f, m_K1s, m_h, m_DRK, m_Ca
    BK::Array{Any,1} #contains Ca_conc, current state probabilities, h_BKT
    PoMET::Array{Any,1} #contains open state prob and deflec due to brownian motion and last step of brownian
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
    P_BKS::typeof(1.0u"L/s") #calcium  dep current
    P_BKT::typeof(1.0u"L/s") #calcium dep current
    P_A::typeof(1.0u"L/s") #IA
    #b is another key control parameter (dimensionless)
    #b controls the strength of the above currents (BKS and BKT)
    b::Float64
    #cell capacitance
    Cm::typeof(1.0u"pF")
end

"""
    Contructor. Note: gK1 and b are the main control parameters of the model.
"""
function Hair_Cell(Vp,
                    #Reversal potentials
                    E_Kp=-95.0, E_hp=-45.0, E_Cap=42.5, E_Lp=0.0,
                    #Maximal Conductances
                    g_K1p=7.0, g_hp=2.2, g_Cap=1.2, g_Lp=0.1,
                    bp=1.0)

    #V0
    V=Vp*1.0u"mV"
    #Reversal Potentials
    E_K=E_Kp*1.0u"mV"
    E_h=E_hp*1.0u"mV"
    E_Ca=E_Cap*1.0u"mV"
    E_L=E_Lp*1.0u"mV"
    E_MET=0.0u"mV"
    #Maximal Conductances
    g_K1=g_K1p*1.0u"nS"
    g_h=g_hp*1.0u"nS"
    g_Ca=g_Cap*1.0u"nS"
    g_L=g_Lp*1.0u"nS"
    g_MET=0.65u"nS"
    #Max permeabilities of currents
    P_DRK=2.4e-14u"L/s" # (a)
    P_BKS=2e-13u"L/s" #calcium current (a,c)
    P_BKT=14e-13u"L/s" #calcium current (a,c)
    P_A=1.08e-13u"L/s"
    #b is another key control parameter (dimensionless)
    #b controls the strength of the above currents (BKS and BKT)
    b=bp
    #cell capacitance
    Cm=10.0u"pF" # (a,c)
    #initial m values for state vector, set as 0 for now
    m_K1f=m_k1_inf(V)
    m_K1s=m_k1_inf(V)
    m_h=m_h_inf(V)
    m_DRK=m_DRK_inf(V)
    m_Ca=m_Ca_inf(V)
    m_A = m_A_inf(V)
    hA1 = h_A_inf(V)
    hA2 = h_A_inf(V)
    #initial values for params of BK-current stuff (another state vector)
    Ca_conc= 0.0u"μmol/L"#find something good
    h_BKT=h_BKT_inf(V) #h for inactivation
    C0=0.9
    C1=0.025
    C2=0.025
    O2=0.025
    O3=0.025

    PoMet=0.15
    #deflec=0.0u"nm"
    brown_def=0.0u"nm"
    last_step=0.0u"nm"

    return Hair_Cell([V,m_K1f, m_K1s, m_h, m_DRK, m_Ca, m_A, hA1, hA2], [Ca_conc, h_BKT,C0,C1,C2,O2,O3],[PoMet,brown_def,last_step],[0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA"],
    E_K,E_h,E_Ca,E_L,E_MET,g_K1,g_h,g_Ca,g_L,g_MET,P_DRK,P_BKS,P_BKT,P_A,b,Cm)
end

function update(cell::Hair_Cell,input=0.0u"nA",deflection=0.0u"nm",stoch=true)
    V=cell.x[1]
    thing=(V*F)/(R*T)
    pow=uconvert(Unitful.NoUnits,thing)
    IK1 = cell.g_K1*(V-cell.E_K)*(0.7*mk1f(V,cell.x[2])+0.3*mk1s(V,cell.x[3]))
    #IK1=0.0u"nA"
    cell.currents[1]=IK1
    Ih = cell.g_h*(V-cell.E_h)*(3*mh(V,cell.x[4])^2*(1-mh(V,cell.x[4]))+mh(V,cell.x[4])^3)
    cell.currents[2]=Ih
    IDRK = uconvert(u"nA",cell.P_DRK*((V*F^2)/(R*T))*((k_int-k_ext*exp(-pow))/(1-exp(-pow)))*(mDRK(V,cell.x[5])^2))
    cell.currents[3]=IDRK
    ICa = cell.g_Ca*(V-cell.E_Ca)*mCa(V,cell.x[6])^3
    cell.currents[4]=ICa
    IBKS = cell.b*cell.P_BKS*((V*F^2)/(R*T))*((k_int-k_ext*exp(-pow))/(1-exp(-pow)))*(O_2(V,cell.BK[5],cell.BK[6],cell.BK[7],cell.BK[1])+O_3(V,cell.BK[6],cell.BK[7],cell.BK[1]))
    cell.currents[5]=IBKS
    IBKT = cell.b*cell.P_BKT*((V*F^2)/(R*T))*((k_int-k_ext*exp(-pow))/(1-exp(-pow)))*(O_2(V,cell.BK[5],cell.BK[6],cell.BK[7],cell.BK[1])+O_3(V,cell.BK[6],cell.BK[7],cell.BK[1]))*hBKT(V,cell.BK[2])
    cell.currents[6]=IBKT
    IL = cell.g_L*(V-cell.E_L)
    cell.currents[7]=IL
    deflec=deflection # deflection = imposed defection
    # now_brow_deflec=cell.PoMET[2]
    # if stoch
    #     cell.PoMET[3]=ΔX(cell.PoMET[3])
    #     cell.PoMET[2]+=cell.PoMET[3]#+now_brow_deflec #add next step of brownian to cumulative
    #     deflec+=cell.PoMET[2] #add new cumulative to deflec
    #     if cell.PoMET[2]>900.0u"nm"
    #         cell.PoMET[2]=900.0u"nm"
    #     elseif cell.PoMET[2]<-450u"nm"
    #         cell.PoMET[2]=-450u"nm"
    #     end
    # end
    cell.PoMET[1]=p_open(deflec) # calc open-state prob based on semi-brownian deflec
    IMET = cell.g_MET*cell.PoMET[1]*(V-cell.E_MET)
    cell.currents[8]=IMET
    aIA=a1(V)
    IA=cell.P_A*((V*F^2)/(R*T))*((k_int-k_ext*exp(-pow))/(1-exp(-pow)))*((mA(V,cell.x[7]))^3)*(aIA*hA1(V,cell.x[8])+(1-aIA)hA2(V,cell.x[9]))
    cell.currents[9]=IA

    ΔV=uconvert(u"mV",(-IK1-Ih-IDRK-ICa-IBKS-IBKT-IL-IMET-IA+input)*dt/cell.Cm)

    #update things
    cell.x[1]= V + ΔV #voltage/membrane potential
    cell.x[2]= mk1f(V,cell.x[2]) #m_K1f
    cell.x[3]= mk1s(V,cell.x[3]) #m_K1s
    cell.x[4]=   mh(V,cell.x[4]) #m_h
    cell.x[5]= mDRK(V,cell.x[5]) #m_DRK
    cell.x[6]=  mCa(V,cell.x[6]) #m_Ca
    cell.x[7]=   mA(V,cell.x[7]) #mA
    cell.x[8]=  hA1(V,cell.x[8]) #hA1
    cell.x[9]=  hA2(V,cell.x[9]) #hA2
    cell.BK[1]= Ca_conc(cell.BK[1],ICa) #Intracellular calcium concentration
    cell.BK[2]= hBKT(V,cell.BK[2]) #h_BKT
    cell.BK[3]= C_0(cell.BK[4],cell.BK[5],cell.BK[6],cell.BK[7]) #C0
    cell.BK[4]= C_1(V,cell.BK[3],cell.BK[4],cell.BK[5],cell.BK[1]) #C1
    cell.BK[5]= C_2(V,cell.BK[4],cell.BK[5],cell.BK[6],cell.BK[1]) #C2
    cell.BK[6]= O_2(V,cell.BK[5],cell.BK[6],cell.BK[7],cell.BK[1]) #02
    cell.BK[7]= O_3(V,cell.BK[6],cell.BK[7],cell.BK[1]) #03
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
const K1_0=6.0u"μmol/L" #micromolar
const K2_0=45.0u"μmol/L" #micromolar
const K3_0=20.0u"μmol/L" #micromolar
const z = 2 #sign of charge of Ca2+
const δ1 = 0.2 #fraction of the electric field experienced by Ca2 at the 1st binding site
const δ2 = 0.0 #fraction of the electric field experienced by Ca2 at the 2nd binding site
const δ3 = 0.2 #fraction of the electric field experienced by Ca2 at the 3rd binding site
const V_A = 30.0u"mV" #potential used to express the voltage dependence of ac

α_c(V)=α_c_0*exp(-V/V_A)

K1(V)=K1_0*exp(δ1*z*uconvert(Unitful.NoUnits,(V*F)/(R*T)))
K2(V)=K2_0*exp(δ2*z*uconvert(Unitful.NoUnits,(V*F)/(R*T)))
K3(V)=K3_0*exp(δ3*z*uconvert(Unitful.NoUnits,(V*F)/(R*T)))

k1(V)=k_1/(K1(V))
k2(V)=k_2/(K2(V))
k3(V)=k_3/(K3(V))

C_0(C1,C2,O2,O3) = 1.0 - (C1 + C2 + O2 + O3)
C_1(V,C0,C1,C2,Ca)= C1 + ΔC_1(V,C0,C1,C2,Ca)
ΔC_1(V,C0,C1,C2,Ca)=uconvert(Unitful.NoUnits,(k1(V)*Ca*C0+k_2*C2-(k_1+k2(V)*Ca)*C1)*dt)
C_2(V,C1,C2,O2,Ca)= C2 + ΔC_2(V,C1,C2,O2,Ca)
ΔC_2(V,C1,C2,O2,Ca)=uconvert(Unitful.NoUnits,(k2(V)*Ca*C1+α_c(V)*O2-(k_2+β_c)*C2)*dt)
O_2(V,C2,O2,O3,Ca)= O2 + ΔO_2(V,C2,O2,O3,Ca)
ΔO_2(V,C2,O2,O3,Ca)=uconvert(Unitful.NoUnits,(β_c*C2+k_3*O3-(α_c(V)+k3(V)*Ca)*O2)*dt)
O_3(V,O2,O3,Ca)= O3+ΔO_3(V,O2,O3,Ca)
ΔO_3(V,O2,O3,Ca)=uconvert(Unitful.NoUnits,(k3(V)*Ca*O2-k_3*O3)*dt)

Ca_conc(Ca,ICa)=Ca+ΔCa_conc(Ca,ICa)
const U=0.02 #free calcium proportion
const Cvol = 1256.637u"μm^3"#1.25u"pL"#cell volume
const ξ = 3.4e-5#proportion of cell hoarding calcium
const Ksca = 2800u"1/s" #rate constant

ΔCa_conc(Ca,ICa)=uconvert(u"μmol/L",(-(U*ICa)/(z*F*Cvol*ξ)-Ksca*Ca)*dt)

#inactivation of BKT current
hBKT(V,hBKT) = hBKT+Δh_BKT(V,hBKT)
Δh_BKT(V,hBKT)= (h_BKT_inf(V)-hBKT)*(dt/τ_BKT(V))
h_BKT_inf(V) = (1+exp((V+61.6u"mV")/3.65u"mV"))^(-1)
τ_BKT(V) = 2.1u"ms"+9.4u"ms"*exp(-(((V+66.9u"mV")/17.7u"mV")^2))

#Equations for IA current
mA(V,mA)=mA+Δm_A(V,mA)
Δm_A(V,mA)=(m_A_inf(V)-mA)*dt/τ_A(V)
m_A_inf(V)=(1+exp(-(V+61u"mV")/10.7u"mV"))^(-1)
function τ_A(V)
    if ustrip(V)<-55
        a=173u"ms"
        b=17.9u"mV"
        c=5.4u"ms"
    else
        a=0.48u"ms"
        b=-23.1u"mV"
        c=1.14u"ms"
    end
    return a*exp(V/b)+c
end

#Inactivation of A-current
hA1(V,hA1) = hA1+Δh_A1(V,hA1)
Δh_A1(V,hA1)= (h_A_inf(V)-hA1)*(dt/τ_hA1(V))

hA2(V,hA2) = hA2+Δh_A2(V,hA2)
Δh_A2(V,hA2)= (h_A_inf(V)-hA2)*(dt/300u"ms")

h_A_inf(V) = (1+exp((V+83u"mV")/3.9u"mV"))^(-1)
τ_hA1(V) = 74u"ms"+321u"ms"*exp(-((V+82u"mV")/15.5u"mV")^2)

a1(V)=(1-0.54)/(1+exp((V+21.5u"mV")/15.6u"mV")) + 0.54

#Equation for open state probability based on deflection

# solve p₀(x₀)= 1/2 (deflection when open state prob = 1/2)
const x₀ =  kβ*T*log( (1-pr)/pr)/zgate
#open state probability based on deflection
p_open(x) = 1.0/(1.0 + exp(-zgate*(x-x₀)/(kβ*T)))

#Brownian Motion of deflection (x is brownian of prev)
#getting params
const Q=1.0u"s^(-2)"#u"m^2 /s^(-1)" #thermal noise power very unsure about units
const a=exp(-dt/2u"ms") #effect of previous timestep on next
const σ=sqrt(Q*(dt/(1-a^2))) #effect of random noise on next timestep
const noise=Normal(0,1)

#actual equation
ΔX(x)=(a*x + σ*sqrt(dt)*rand(noise)*1.0u"nm")

"""
    runs for a given time to let parameters settle into equilibrium state
"""
function burnin(cell, time=5.0u"s", displayV=true,displayEachC=false,displayallC=true,printMP=true)
    t=0.0u"s":dt:time
    voltBI=Array{typeof(1.0u"mV")}(undef, length(t))
    #Ca_p=Array{typeof(1.0u"μmol/L")}(undef, length(t))
    #POPEN=Array{typeof(1.0u"nm")}(undef, length(t))
    currents_p=Array{Any,2}(undef,length(t),9)
    for i in 1:length(t)
        update(cell)
        voltBI[i]=cell.x[1]
        currents_p[i,:].=cell.currents
    #    Ca_p[i]=cell.BK[1]
    #    POPEN[i]=cell.PoMET[2]
    end
    #display(plot(ustrip.(t),ustrip.(Ca_p), title="Calcium"))
    #tot=0.0u"nm"
    #for x in POPEN
    #    tot+=x
    #end
    #println(tot)
    #println(string("average displacement: ",mean(sqrt.(POPEN[:].^2))))
    #println(string("standard dev: ",std(POPEN)))
    #display(plot(ustrip.(t),ustrip.(POPEN), title="Brownian Def"))
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
    display(plot(t,currents[:,9], title="IA"))
end

"""
    Called by burnin to plot a graph of all currents over time
"""
function plotTOGETHER(t,currents)
    plot(t,currents[:,1], label="IK1",color=:red)
    plot!(t,currents[:,2], label="Ih", color=:blue)
    plot!(t,currents[:,3], label="IDRK",color=:blue, linestyle=:dash)
    plot!(t,currents[:,4], label="ICa",color=:gold)
    plot!(t,currents[:,5], label="IBKS",color=:green)
    plot!(t,currents[:,6], label="IBKT",color=:green, linestyle=:dash)
    plot!(t,currents[:,7], label="IL",color=:gold, linestyle=:dash)
    plot!(t,currents[:,8], label="IMET",color=:black)
    display(plot!(t,currents[:,9], label="IA",color=:red, linestyle=:dash))
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
    currents_p=Array{Any,2}(undef,length(t),9)
    Ca_p=Array{typeof(1.0u"μmol/L")}(undef, length(t))

    for i in 1:length(input)
        update(cell,input[i])
        voltages[i]=helga.x[1]
        currents_p[i,:].=cell.currents
        Ca_p[i]=cell.BK[1]
    end
    display(plot(ustrip.(t),ustrip.(Ca_p), title="Calcium"))
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

function go_once()
    #v0: -20 - -60
    V0=(rand()-0.5)*40 - 40
    #E_K: -75 - -100
    E_K=(rand()-0.5)*25 - 87.5
    #E_h: -38 - -48
    E_h=(rand()-0.5)*10 - 43
    #E_Ca: -39 - -43
    E_Ca=(rand()-0.5)*4 - 41
    #E_L: 0 - -50
    E_L=(rand()-0.5)*50 - 25
    #g_K: 0 - 40
    g_K=(rand()-0.5)*40 + 20
    #g_h: 1 - 6
    g_h=(rand()-0.5)*5 + 3.5
    #g_Ca: 1-4.5
    g_Ca=(rand()-0.5)*3.5 + 2.75
    #g_L: 0.1-0.5
    g_L=(rand()-0.5)*0.4 + 0.2
    #b: 0-1
    b=rand()
    params=[V0,E_K,E_h,E_Ca,E_L,g_K,g_h,g_Ca,g_L,b]
    key= ["V0","E_K","E_h","E_Ca","E_L","g_K","g_h","g_Ca","g_L","b"]

    romy = Hair_Cell(V0,E_K,E_h,E_Ca,E_L,g_K,g_h,g_Ca,g_L,b)

    #burnin(cell, time=5.0u"s", displayV=true,displayEachC=false,displayallC=true,printMP=true)

    burnin(romy,2.0u"s",true,false,true,true)
    return ustrip(romy.x[1]),params
end

function go_lots()
    n_trials=5
    params=Array{Any,2}(undef,n_trials,10)
    voltages=zeros(n_trials)
    i=1
    while i<= n_trials
        stuff=go_once()
        params[i,:].=stuff[2]
        voltages[i]=stuff[1]
        i+=1
    end
    fpath="C:/Users/josep/OneDrive/Desktop/Paramter exploration/"
    overall=histogram(voltages,nbins=10,legend=false,ylabel="Number of Trials", xlabel="Resting Membrane Potential (mV)", title="Distribution of  RMP with Random Paramter Allocation")
    display(overall)
    savefig(overall,string(fpath,"overall"))
    v0=scatter(params[:,1],voltages[:], xlabel="Parameter Value (mV)", ylabel="RMP (mV)", title="V0",legend=false)
    display(v0)
    savefig(v0,string(fpath,"v0"))
    Ek=scatter(params[:,2],voltages[:], xlabel="Parameter Value (mV)", ylabel="RMP (mV)", title="E_K",legend=false)
    display(Ek)
    savefig(Ek,string(fpath,"EK"))
    fname=string(n_trials," trials")
    file=open(string(fpath,fname),"w")
    x=1
    while x<= n_trials
        for j in 1:10
            write(file,string(params[x,j],";"))
        end
        write(file,string(voltages[x],"\n"))
        x+=1
    end
    close(file)
end

go_lots()
