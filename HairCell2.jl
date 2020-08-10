"""
Based on Nieman 2011
Model of bullfrog saccular hair cell
Hodgkin-Huxely type
"""

using Unitful

const k_int = 112.0u"mmol/L" #Intracellular potassium conc. (millimolar)
const k_ext = 2.0u"mmol/L" #Extracellular potassium conc. (millimolar)
const F = 96485.3329u"s*A/mol" #Faraday constant
const R = 8.314u"J/(mol*K)" #universal gas constant
const T = 295.15u"K" #Temperature in Kalvin
const dt= 1e-6u"s" #timestep
const e = 2.7182818284590452353602874713527 #euler's constant


struct Hair_Cell
    #state vector(s)
    x::Array{Any,1} #contains V, m_K1f, m_K1s, m_h, m_DRK, m_Ca
    BK::Array{Any,1} #contains Ca_conc, h_BKT
    #Parameters
    #
    #Reversal Potentials
    E_K::typeof(1.0u"mV")
    E_h::typeof(1.0u"mV") #for cation h-current
    E_Ca::typeof(1.0u"mV")
    E_L::typeof(1.0u"mV")
    #Maximal Conductances
    g_K1::typeof(1.0u"nS") #g_K1 is one of the main control parameters
    g_h::typeof(1.0u"nS") #for cation h-current
    g_Ca::typeof(1.0u"nS")
    g_L::typeof(1.0u"nS")
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
Note: gK1 and b are the main control parameters of the model.
"""
function Hair_Cell(V::typeof(1.0u"mV"),g_K1=20.0u"nS",b=0.50)
    #Reversal Potentials
    E_K=-95.0u"mV"
    E_h=-45.0u"mV" #for cation h-current
    E_Ca=42.5u"mV"
    E_L=-60.0u"mV" #changed from paper bc I think they made a mistake
    #Maximal Conductances
    g_K1=g_K1#g_K1 is one of the main control parameters of the model
    g_h=2.2u"nS" #for cation h-current
    g_Ca=1.2u"nS"
    g_L=0.77u"nS" #not in paper, so took value from Canelli et al*****************
    #Max permeabilities of currents
    P_DRK=2.4e-14u"L/s"
    P_BKS=2e-13u"L/s" #calcium current
    P_BKT=14e-13u"L/s" #calcium current
    #b is another key control parameter (dimensionless)
    #b controls the strength of the above currents (BKS and BKT)
    b=b
    #cell capacitance
    Cm=10.0u"pF"
    #initial m values for state vector, set as 0 for now
    m_K1f=0
    m_K1s=0
    m_h=0
    m_DRK=0
    m_Ca=0
    #initial values for params of BK-current stuff (another state vector)
    Ca_conc= 0.0u"mmol/L"#find something good
    h_BKT=0 #h for inactivation
    C0=0.0
    C1=0.0
    C2=0.0
    O2=0.0
    O3=0.0
    return Hair_Cell([V,m_K1f, m_K1s, m_h, m_DRK, m_Ca], [Ca_conc, h_BKT,C0,C1,C2,O2,O3],
    E_K,E_h,E_Ca,E_L,g_K1,g_h,g_Ca,g_L,P_DRK,P_BKS,P_BKT,b,Cm)
end

function update(cell::Hair_Cell,input=0.0u"nA")
    V=cell.x[1]
    thing=(V*F)/(R*T)
    pow=uconvert(Unitful.NoUnits,thing)
    IK1 = cell.g_K1*(V-cell.E_K)*(0.7*mk1f(V,cell.x[2])+0.3*mk1s(V,cell.x[3]))
    Ih = cell.g_h*(V-cell.E_h)*(3*mh(V,cell.x[4])^2*(1-mh(V,cell.x[4]))+mh(V,cell.x[4])^3)
    IDRK = uconvert(u"nA",cell.P_DRK*((V*F^2)/(R*T))*(k_int-k_ext*e^-pow)/(1-e^-pow)*mDRK(V,cell.x[5])^2)
    ICa = cell.g_Ca*(V-cell.E_Ca)*mCa(V,cell.x[6])^3
    IBKS = cell.b*cell.P_BKS*((V*F^2)/(R*T))*((k_int-k_ext*e^-pow)/(1-e^-pow))*(O_2(cell.BK[1],cell.BK[5],cell.BK[6],cell.BK[7])+O_3(cell.BK[1],cell.BK[6],cell.BK[7]))
    IBKT = cell.b*cell.P_BKT*((V*F^2)/(R*T))*((k_int-k_ext*e^-pow)/(1-e^-pow))*(O_2(cell.BK[1],cell.BK[5],cell.BK[6],cell.BK[7])+O_3(cell.BK[1],cell.BK[6],cell.BK[7]))*hBKT(V,cell.BK[2])
    IL = cell.g_L*(V-cell.E_L)
    IMET = 0.0u"nA" #while we do not care
    ΔV=uconvert(u"mV",(-IK1-Ih-IDRK-ICa-IBKS-IBKT-IL-IMET)*dt/cell.Cm)

    #update things
    cell.x[1]= V + ΔV #voltage/membrane potential
    cell.x[2]= mk1f(V,cell.x[2]) #m_K1f
    cell.x[3]= mk1s(V,cell.x[3]) #m_K1s
    cell.x[4]=   mh(V,cell.x[4]) #m_h
    cell.x[5]= mDRK(V,cell.x[5]) #m_DRK
    cell.x[6]=  mCa(V,cell.x[6]) #m_Ca
    cell.BK[1]= Ca_conc(cell.BK[1],ICa) #Intracellular calcium concentration
    cell.BK[2]= hBKT(V,cell.BK[2]) #h_BKT
    cell.BK[3]= C_0(cell.BK[4],cell.BK[5],cell.BK[6],cell.BK[7]) #C0
    cell.BK[4]= C_1(cell.BK[1],cell.BK[3],cell.BK[4]) #C1
    cell.BK[5]= C_2(cell.BK[1],cell.BK[4],cell.BK[5],cell.BK[6]) #C2
    cell.BK[6]= O_2(cell.BK[1],cell.BK[5],cell.BK[6],cell.BK[7]) #02
    cell.BK[7]= O_3(cell.BK[1],cell.BK[6],cell.BK[7]) #03
end

#equations for K1 current
mk1f(V,mk1f)=ustrip(mk1f)+Δm_k1f(V,mk1f)
Δm_k1f(V,mk1f)=(m_k1_inf(V)-ustrip(mk1f))*ustrip(dt)/τ_k1f(V)
τ_k1f(V)=0.7exp(-(ustrip(V) + 120)/43.8)+0.04

mk1s(V,mk1s)=ustrip(mk1s)+Δm_k1s(V,mk1s)
Δm_k1s(V,mk1s)=(m_k1_inf(V)-ustrip(mk1s))*ustrip(dt)/τ_k1s(V)
τ_k1s(V)=14.1exp(-(ustrip(V) + 120)/28)+0.04

m_k1_inf(V)=(1+exp((ustrip(V) + 110)/11))^-1

#equations for h cation current
mh(V,mh)=ustrip(mh)+Δm_h(V,mh)
Δm_h(V,mh)=(m_h_inf(V)-ustrip(mh))*ustrip(dt)/τ_h(V)
m_h_inf(V)=(1+exp((ustrip(V)+87)/16.7))^(-1)
τ_h(V)=63.7+135.7exp(-((ustrip(V)+91.4)/21.2)^2)

#equations for DRK current
mDRK(V,mDRK)= ustrip(mDRK)+ΔmDRK(V,mDRK)
ΔmDRK(V, mDRK) = (m_DRK_inf(V)- ustrip(mDRK))*(ustrip(dt)/τ_DRK(V))
m_DRK_inf(V) = (1+exp((ustrip(V)+48.3)/4.19))^(-1/2)
#mDRK(V,mDRK)= ustrip(mDRK)+ΔmDRK(V,mDRK)
#ΔmDRK(V, mDRK) = (m_DRK_inf(V)- ustrip(mDRK))*(ustrip(dt)/τ_DRK(V))
#m_DRK_inf(V) = (1+exp((ustrip(V)+48.3)/4.19))^(-1/2)
τ_DRK(V) = (α_DRK(V)+β_DRK(V))^(-1)
α_DRK(V) = (3.2*e^(-ustrip(V)/20.9)+3)^(-1)
β_DRK(V) = (1467*e^(ustrip(V)/5.96)+9)^(-1)

#equations for ICa current
mCa(V,mCa) = ustrip(mCa)+Δm_Ca(V,mCa)
Δm_Ca(V, mCa) = (m_Ca_inf(V)-ustrip(mCa))*(ustrip(dt)/τ_Ca(V))
m_Ca_inf(V) = (1+exp(-(ustrip(V)+55)/12.2))^(-1)
τ_Ca(V) = 0.046 + 0.325(exp(-((ustrip(V)+77)/51.67)^2))

#equations for BK currents:
#Params:
const k_1=300.0#u"1/s" #s^-1
const k_2=5000.0#u"1/s" #s^-1
const k_3=1500.0#u"1/s" #s^-1
const β_c=2500.0#u"1/s" #s^-1
#below this may actually be intended to be functions
const α_c=450.0#u"1/s" #s^-1
const k1=6.0#u"μmol/L" #micromolar
const k2=45.0#u"μmol/L" #micromolar
const k3=20.0#u"μmol/L" #micromolar

C_0(C1,C2,O2,O3)= 1.0 - (C1 + C2 +O2 +O3)
C_1(Ca,C0,C1)= C1 + ΔC_1(Ca,C0,C1)
ΔC_1(Ca,C0,C1)=(k1*ustrip(Ca)*C0+k_2*C1-(k_1+k2*ustrip(Ca))*C1)*ustrip(dt)
C_2(Ca,C1,C2,O2)= C2 + ΔC_2(Ca,C1,C2,O2)
ΔC_2(Ca,C1,C2,O2)=(k2*ustrip(Ca)*C1+α_c*O2-(k_2+β_c)*C2)*ustrip(dt)
O_2(Ca,C2,O2,O3)= O2 + ΔO_2(Ca,C2,O2,O3)
ΔO_2(Ca,C2,O2,O3)=(β_c*C2+k_3*O3-(α_c+k3*ustrip(Ca))*O2)*ustrip(dt)
O_3(Ca,O2,O3)= O3+ΔO_3(Ca,O2,O3)
ΔO_3(Ca,O2,O3)=(k3*ustrip(Ca)*O2-k_3*O3)*ustrip(dt)

#really wanna give this units
Ca_conc(Ca,ICa)=ustrip(Ca)+ΔCa_conc(Ca,ICa)
ΔCa_conc(Ca,ICa)=(-0.00061*ustrip(ICa) - 2800*ustrip(Ca))*ustrip(dt) #this smells funky
#Ca_conc(Ca,ICa)=uconvert(u"mmol/L",Ca+ΔCa_conc(Ca,ICa))
#ΔCa_conc(Ca,ICa)=(-0.00061*ICa - 2800*Ca)*dt

#inactivation of BKT current
hBKT(V,hBKT) = ustrip(hBKT)+Δh_BKT(V,hBKT)
Δh_BKT(V,hBKT)= (h_BKT_inf(V)-ustrip(hBKT))*(ustrip(dt)/τ_BKT(V))
h_BKT_inf(V) = (1+exp((ustrip(V)+61.6)/3.65))^(-1)
τ_BKT(V) = 2.1+9.4*exp(-((ustrip(V)+66.9)/17.7)^2)

helga = Hair_Cell(-60.0u"mV")
println(helga.x)

for i in 1:10000
    update(helga)
    if i%100 == 0
        println(helga.x[1])
    end
end
println(helga.x)
