"""
Based on Nieman 2011
Model of bullfrog saccular hair cell
Hodgkin-Huxely type
"""

using Unitful
#using Plots
using Makie
using MakieLayout
using Printf
#gr()

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

###################################################################
#Construct GUI
###################################################################


scene, layout=layoutscene(nrows=1,ncols=3,resolution = (1400,700))

#set up axes for burn-in plot
layout[1,1]=BIplots_layout=GridLayout(2,1)
BIplots_layout[1,1]=BI_axis=LAxis(scene, title="Burn-In", ylabel="memb. pot. (mV)", xlabel="time (sec)")

#set up axes for pulse response plot
layout[1,3]=PRplots_layout=GridLayout(2,1)
PRplots_layout[1,1]=PR_axis=LAxis(scene, title="Pulse Response", ylabel="memb. pot. (mV)", xlabel="time (sec)")

#slider panel
layout[1,2]=slider_panel=GridLayout(13,1)
#EK slider
slider_panel[1,1]=EK_layout=GridLayout(1,3)
EK_layout[1,1]=LText(scene,"E_K", textsize=20)
EK_layout[1,2]=E_K_slider=LSlider(scene, range=-100.0:0.2:-70.0,startvalue=-95.0)
EK_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), E_K_slider.value))

#Eh slider
slider_panel[2,1]=Eh_layout=GridLayout(1,3)
Eh_layout[1,1]=LText(scene,"E_h", textsize=20)
Eh_layout[1,2]=E_h_slider=LSlider(scene, range=-50.0:0.1:-35.0,startvalue=-45.0)
Eh_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), E_h_slider.value))


#E_Ca slider
slider_panel[3,1]=ECa_layout=GridLayout(1,3)
ECa_layout[1,1]=LText(scene,"E_Ca", textsize=20)
ECa_layout[1,2]=E_Ca_slider=LSlider(scene, range=40.0:0.5:100.0,startvalue=42.5)
ECa_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), E_Ca_slider.value))

#EL slider
slider_panel[4,1]=EL_layout=GridLayout(1,3)
EL_layout[1,1]=LText(scene,"E_L", textsize=20)
EL_layout[1,2]=E_L_slider=LSlider(scene, range=-70.0:0.0,startvalue=-40.0)
EL_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), E_L_slider.value))

#gK1 slider
slider_panel[5,1]=gK1_layout=GridLayout(1,3)
gK1_layout[1,1]=LText(scene,"g_K1", textsize=20)
gK1_layout[1,2]=g_K1_slider=LSlider(scene, range=3.0:0.5:42.0,startvalue=20.0)
gK1_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), g_K1_slider.value))

#gh slider
slider_panel[6,1]=gh_layout=GridLayout(1,3)
gh_layout[1,1]=LText(scene,"g_h", textsize=20)
gh_layout[1,2]=g_h_slider=LSlider(scene, range=1.0:0.05:6.0,startvalue=2.2)
gh_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), g_h_slider.value))

#gCa slider
slider_panel[7,1]=gCa_layout=GridLayout(1,3)
gCa_layout[1,1]=LText(scene,"g_Ca", textsize=20)
gCa_layout[1,2]=g_Ca_slider=LSlider(scene, range=1.0:0.01:5.0,startvalue=1.2)
gCa_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), g_Ca_slider.value))

#gL slider
slider_panel[8,1]=gL_layout=GridLayout(1,3)
gL_layout[1,1]=LText(scene,"g_L", textsize=20)
gL_layout[1,2]=g_L_slider=LSlider(scene, range=0.1:0.01:0.9,startvalue=0.77)
gL_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), g_L_slider.value))

#b slider
slider_panel[9,1]=b_layout=GridLayout(1,3)
b_layout[1,1]=LText(scene,"b (BK)", textsize=20)
b_layout[1,2]=b_slider=LSlider(scene, range=0.01:0.01:1.0,startvalue=0.5)
b_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), b_slider.value))

#v0 slider
slider_panel[10,1]=v0_layout=GridLayout(1,3)
v0_layout[1,1]=LText(scene,"v0", textsize=20)
v0_layout[1,2]=v0_slider=LSlider(scene, range=-120.0:40.0,startvalue=-50.0)
v0_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), v0_slider.value))

#beta_c slider
slider_panel[11,1]=β_c_layout=GridLayout(1,3)
β_c_layout[1,1]=LText(scene,"β_c (BK)", textsize=20)
β_c_layout[1,2]=β_c_slider=LSlider(scene, range=1000:100:2500,startvalue=2500)
β_c_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), β_c_slider.value))

#pulse size slider
slider_panel[12,1]=p_layout=GridLayout(1,3)
p_layout[1,1]=LText(scene,"Pulse", textsize=20)
p_layout[1,2]=p_slider=LSlider(scene, range=0.005:0.005:0.1,startvalue=0.005)
p_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), p_slider.value))

#Redraw button
slider_panel[13,1]=button_grid=GridLayout(1,2)
button_grid[1,1]=go_button=LButton(scene,label="Redraw")


################################################################################
#Done constructing GUI
################################################################################

"""
    Contructor. Note: gK1 and b are the main control parameters of the model.
"""
function Hair_Cell(V::typeof(1.0u"mV"))#,g_K1=100.0u"nS",b=0.7) # g and b from (a)
    #Reversal Potentials
    E_K=E_K_slider.value[]*1.0u"mV" #-90.0u"mV" # (a,c)
    E_h=E_h_slider.value[]*1.0u"mV" #-40.0u"mV" #for cation h-current (e)
    E_Ca=E_Ca_slider.value[]*1.0u"mV" #42.5u"mV" # (a,f)
    E_L=E_L_slider.value[]*1.0u"mV" #-40.0u"mV" # (f)
    E_MET=0.0u"mV" # (a)
    #Maximal Conductances
    g_K1=g_K1_slider.value[]*1.0u"nS" #g_K1#g_K1 is one of the main control parameters of the model
    g_h=g_h_slider.value[]*1.0u"nS" #2.0u"nS" #for cation h-current (e)
    g_Ca=g_Ca_slider.value[]*1.0u"nS" #1.3u"nS" # (d,e)
    g_L=g_L_slider.value[]*1.0u"nS" #0.77u"nS" # (b)
    g_MET=0.65u"nS" # (a)
    #Max permeabilities of currents
    P_DRK=2.4e-14u"L/s" # (a)
    P_BKS=2e-13u"L/s" #calcium current (a,c)
    P_BKT=14e-13u"L/s" #calcium current (a,c)
    #b is another key control parameter (dimensionless)
    #b controls the strength of the above currents (BKS and BKT)
    b=b_slider.value[]
    #cell capacitance
    Cm=10.0u"pF" # (a,c)
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

"""
function taking a hair cell and a current, mutates hair cell based on effect of current
"""
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
β_c()=β_c_slider.value[]*1.0u"1/s"#1000.0u"1/s"#2500.0u"1/s" #s^-1 #changed to OG value *******************************

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
ΔC_2(V,C1,C2,O2)=uconvert(Unitful.NoUnits,(k2_Ca(V)*C1+α_c(V)*O2-(k_2+β_c())*C2)*dt)
O_2(V,C2,O2,O3)= O2 + ΔO_2(V,C2,O2,O3)
ΔO_2(V,C2,O2,O3)=uconvert(Unitful.NoUnits,(β_c()*C2+k_3*O3-(α_c(V)+k3_Ca(V))*O2)*dt)
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
        #display(plot(ustrip.(t),ustrip.(voltBI), title="Membrane Potential during Burn-in"))
    end
    if displayEachC
        #plotEach(ustrip.(t),ustrip.(currents_p))
    end
    if displayallC
        #plotTOGETHER(ustrip.(t),ustrip.(currents_p))
    end
    if printMP
        #println(string("RMP after burn-in: ",cell.x[1]))
    end
    return ustrip.(voltBI),ustrip.(t)
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
function pulse(time=0.5u"s", start=2.0e-1u"s", len=2.0e-1u"s", amplitude=p_slider.value[]*1.0u"nA")
    t=0.0u"s":dt:time
    u = Array{typeof(1.0u"nA")}(undef, length(t))
    u.=0.0u"nA"
    u[findall( t-> (t>=start) & ( t<start+len), t)] .= amplitude

    return u
end

"""
    Models response of cell to a pulse of given amplitude and length
"""
function pulseResponse(cell, time=0.5u"s", start=2.0e-1u"s", len=2.0e-1u"s", amplitude=p_slider.value[]*1.0u"nA", displayV=true,displayEachC=false,displayallC=true,printMP=false)
    t=0.0u"s":dt:time
    input = pulse(time,start,len,amplitude) #in nA
    voltages=Array{typeof(1.0u"mV")}(undef, length(t))
    currents_p=Array{Any,2}(undef,length(t),8)

    for i in 1:length(input)
        update(cell,input[i])
        voltages[i]=cell.x[1]
        currents_p[i,:].=cell.currents
    end
    if displayV
        #plot(ustrip.(t),ustrip.(voltages), title="Membrane Potential during Pulse",ylabel="Membrane Potential (mV)")
        #display(plot!(twinx(),ustrip.(t),ustrip.(input),ylims=(-1e-3,ustrip(amplitude)*3),ylabel="Pulse Amplitude (nA)",linecolor=:violet,label="pulse"))
    end
    if displayEachC
        #plotEach(ustrip.(t),ustrip.(currents_p))
    end
    if displayallC
        #plotTOGETHER(ustrip.(t),ustrip.(currents_p))
    end
    if printMP
        #println(string("RMP after pulse: ",cell.x[1]))
    end
    return voltages,t,cell
end

helga = Hair_Cell(v0_slider.value[]*1.0u"mV")

burn_arrays=burnin(helga,2u"s")
burn_time=burn_arrays[2]
burn_v=burn_arrays[1]
lines!(BI_axis,burn_time,burn_v, color=:grey)

pulse_arrays=pulseResponse(helga)
pulse_time=ustrip.(pulse_arrays[2])
pulse_v=ustrip.(pulse_arrays[1])
lines!(PR_axis,pulse_time,pulse_v, color=:grey)
#lines!(PR_axis,1:10,1:10)


on(go_button.clicks) do clicks
    println("should change")
    #lines!(BI_axis,burn_time1,burn_v1,color=:grey)
    romy = Hair_Cell(v0_slider.value[]*1.0u"mV")
    burn_arrays=burnin(romy,2u"s")
    burn_time=burn_arrays[2]
    burn_v=burn_arrays[1]
    #BIplots_layout[1,1]=BI_axis2=LAxis(scene, title="Burn-In", ylabel="memb. pot. (mV)", xlabel="time (sec)")
    lines!(BI_axis,burn_time,burn_v,color=:red)
    pulse_arrays=pulseResponse(romy)
    pulse_time=ustrip.(pulse_arrays[2])
    pulse_v=ustrip.(pulse_arrays[1])
    lines!(PR_axis,pulse_time,pulse_v, color=:red)
    update!(scene)
end



display(scene)
