
#生产能力kg/s
flow=30/24/3600*1000


```julia
#mole fraction
x_D=(0.94/46)/(0.94/46+0.06/18)
x_B=(0.001/46)/(0.001/46+0.999/18)
x_F=(0.25/46)/(0.25/46+0.75/18)
```


    0.11538461538461538


```julia
#minimun stages calculation
T_D=78.2+273.15
T_B=99.8+273.15
p_EtOH_D=10^(7.30243-1630.868/(T_D-43.569))
p_EtOH_B=10^(6.84806-1358.124/(T_B-71.034))
p_water_D=10^(7.11564-1687.537/(T_D-42.98))
p_water_B=10^(7.11564-1687.537/(T_B-42.98))
alpha_D=p_EtOH_D/p_water_D
alpha_B=p_EtOH_B/p_water_B
alpha=sqrt(alpha_D*alpha_B)
N_m=log(x_D/(1-x_D)*(1-x_B)/x_B)/log(alpha)
```


    11.836927128798163


```julia
#Underwood Eq
theta=1.97391645
R_min=(1-x_D)/(1-theta)+alpha*(x_D)/(alpha-theta)-1
R=1.8*R_min
```


    10.117415234239775


```julia
#mass Eq
#R_min from matlab
R_min=2.0645
R=1.8*R_min
D=(flow*0.94/46+flow*0.06/18)*1000;
solve=inv([1 -1 1;x_F -x_B 0;0 0 1])*[D;D*x_D;(R+1)*D]
F=solve[1]
B=solve[2]
V_0=solve[3]
```

```julia
[x_D x_B x_F]
```


    1×3 Array{Float64,2}:
     0.859756  0.000391543  0.115385


```julia
#unit mol/s
[F,V_0,D,B]
```

```julia
#V'
V_prime=V_0
L=R*D
L_prime=L+F
[L L_prime V_0 V_prime]

#=Aspen result
L=7.06
L_prime=21.7
V_0=8.99
V_prime=9.2=#
```



```julia
#idea stages from Grilliland figure
X=(R-R_min)/(R+1)
Y=1-exp((1+54.4*X)*(X-1)/(11+117.2*X)/sqrt(X))
N=(Y+N_m)/(1-Y)
```




    18.596843387766274




```julia
#stage efficiency
E_T=0.5
T=(106.861+78.8242)/2

0.17-0.616log(0.3)
```




    0.9116472474647767




```julia
#surface tension calculation
q=2
T_D=78.8255
T_B=106.815
T_F=89.9305
#surface tension,pure,unit:dyn/cm
sigma_e_D=17.3496144267588
sigma_e_B=14.6413477845619
sigma_e_F=16.2960387828613
sigma_w_D=62.3193571276298
sigma_w_B=56.9180940502544
sigma_w_F=60.1893524530583
#mole volume unit:ccm/mol
V_e_D=1/0.0159383960755079
V_e_B=1/0.0155345167312197
V_e_F=1/0.0156829574583475
V_w_D=1/0.053825020877694
V_w_B=1/0.0532409595776902
V_w_F=1/0.0534594853416629
phi_water_D=(1-x_D)*V_w_D/((1-x_D)*V_w_D+x_D*V_e_D)
phi_ethnaol_D=x_D*V_e_D/((1-x_D)*V_w_D+x_D*V_e_D)
phi_water_B=(1-x_B)*V_w_B/((1-x_B)*V_w_B+x_B*V_e_B)
phi_ethnaol_B=x_B*V_e_B/((1-x_B)*V_w_B+x_B*V_e_B)
phi_water_F=(1-x_F)*V_w_F/((1-x_F)*V_w_F+x_F*V_e_F)
phi_ethnaol_F=x_F*V_e_F/((1-x_F)*V_w_F+x_F*V_e_F)

Q_D=0.441*q/T_D*(sigma_e_D*V_e_D^(2/3)/q-sigma_w_D*V_w_D^(2/3))
Q_B=0.441*q/T_B*(sigma_e_B*V_e_B^(2/3)/q-sigma_w_B*V_w_B^(2/3))
Q_F=0.441*q/T_F*(sigma_e_F*V_e_F^(2/3)/q-sigma_w_F*V_w_F^(2/3))

B_D=log(phi_water_D^q/phi_ethnaol_D)
B_B=log(phi_water_B^q/phi_ethnaol_B)
B_F=log(phi_water_F^q/phi_ethnaol_F)

A_D=B_D+Q_D
A_B=B_B+Q_B
A_F=B_F+Q_F

phi_se_D=(2+10^A_D-sqrt((10^A_D+2)^2-4))/2
phi_sw_D=1-phi_se_D
phi_se_B=(2+10^A_B-sqrt((10^A_B+2)^2-4))/2
phi_sw_B=1-phi_se_B
phi_se_F=(2+10^A_F-sqrt((10^A_F+2)^2-4))/2
phi_sw_F=1-phi_se_F

sigma_m_D=(phi_sw_D*sigma_w_D^0.25+phi_se_D*sigma_e_D^0.25)^4
sigma_m_B=(phi_sw_B*sigma_w_B^0.25+phi_se_B*sigma_e_B^0.25)^4
sigma_m_F=(phi_sw_F*sigma_w_F^0.25+phi_se_F*sigma_e_F^0.25)^4
#compare with Aspen database
sigma_D=23.4028341490826
sigma_B=56.9272057495209
sigma_F=54.3729402955967

[sigma_m_D sigma_m_B sigma_m_F]
```




    1×3 Array{Float64,2}:
     17.3501  56.9145  17.8162




```julia
#properties
#Molar Density,kmol/m^3
MD_D_L=17.566120355233
MD_D_V=0.0355382795548142
MD_B_L=50.5309214121779
MD_B_V=0.0407411639790688
MD_F_L=40.6191292799249
MD_F_V=0.045852963826942

#pressure(bar)/temperature/vapor molewieght/liquid moleweight/rho_vapor/rho_liquid/sigma/V_sm/L_sm(m^3/s)
prop_D=[1.04 78.8255 42.3789 42.2578 1.5051 743.118 sigma_D 0 D/1000/MD_D_L]
prop_F=[1.2 89.9305 31.1044 21.2523 1.23644 868.019 sigma_F 0 F/1000/MD_F_L]
prop_B=[1.28679 106.815 18.0153 18.6153 0.735287 910.956 sigma_B 0 B/1000/MD_B_L]
prop_distill=(prop_D+prop_F)/2
prop_ex=(prop_F+prop_B)/2

#vapor flow(m^3/s) corrected by aspen simulation
prop_distill[8]=39.048/1000/((MD_D_V+MD_F_V)/2)
prop_ex[8]=39.2276/1000/((MD_B_V+MD_F_V)/2)
```

```julia
#精馏段塔径
#两相流动参数
F_LV=prop_distill[9]/prop_distill[8]*(prop_distill[6]/prop_distill[5])^0.5
#气体负荷因子
H_T=0.61
h_I=0.05
H=H_T-h_I
C_f=exp(-4.531+1.6562*H+5.5496*H^2-6.4695*H^3+(-0.474675+0.079*H-1.39*H^2+1.3212*H^3)*log(F_LV)+(-0.07291+0.088307*H-0.49123*H^2+0.43196*H^3)*log(F_LV)^2)
#液泛气速
u_f=C_f*(prop_distill[6]/prop_distill[5]-1)^0.5
u=0.7*u_f
#column diameter
D=sqrt(4*prop_distill[8]/pi/u)
```

```julia
#提馏段塔径
#两相流动参数
F_LV=prop_ex[9]/prop_ex[8]*(prop_ex[6]/prop_ex[5])^0.5
#气体负荷因子
H_T=0.61
h_I=0.05
H=H_T-h_I
C_f=exp(-4.531+1.6562*H+5.5496*H^2-6.4695*H^3+(-0.474675+0.079*H-1.39*H^2+1.3212*H^3)*log(F_LV)+(-0.07291+0.088307*H-0.49123*H^2+0.43196*H^3)*log(F_LV)^2)
#液泛气速
u_f=C_f*(prop_ex[6]/prop_ex[5]-1)^0.5
u=0.7*u_f
#column diameter
D=sqrt(4*prop_ex[8]/pi/u)
```

```julia
#塔高计算
H_b=4*60*prop_ex[9]/(pi*0.8^2/4)
n=38
n_f=2
n_p=4
H_t=0.61
H_d=0.72
H_b=0.8
H_f=0.61
H_p=0.5
H_all=(n-n_f-n_p-1)*H_t+H_d+H_b+H_f*n_f+H_p*n_p
```



```julia
#塔板结构尺寸
D=0.8
L=0.529
H=0.1

#堰上液层厚度
h_pw_dist=2.84/1000*1*(prop_distill[9]*3600/L)^(2/3)
#出口堰高度
h_w=0.04
#底隙高度
h_o=0.0263
#底隙处流速
u_o_prime_distill=prop_distill[9]/L/h_o
u_o_prime_ex=prop_ex[9]/L/h_o
#阀孔动能因子
u_o_distill=prop_distill[8]/(pi/4*0.8^2*0.109)
F_o_distill=u_o*sqrt(prop_distill[5]*0.7)
u_o_ex=prop_ex[8]/(pi/4*0.8^2*0.109)
F_o_ex=u_o*sqrt(prop_ex[5]*0.7)
#临界孔速
u_oc=(73.1/prop_distill[5])^(1/1.875)
#精馏段单板压降
h_c_distill=5.34*prop_distill[5]*u_o^2/2/prop_distill[6]/9.81
h_l_distill=0.5*(h_w+h_o)##dyn/cm=1mN/m
h_sigma_distill=2*prop_distill[7]/1000/0.0127/prop_distill[6]/9.81
drop_p_distill=prop_distill[6]*9.81*(h_c_distill+h_l_distill+h_sigma_distill)
#提馏段单板压降
h_c_ex=5.34*prop_ex[5]*u_o^2/2/prop_ex[6]/9.81
h_l_ex=0.5*(h_w+h_o)##dyn/cm=1mN/m
h_sigma_ex=2*prop_ex[7]/1000/0.0127/prop_ex[6]/9.81
drop_p_ex=prop_ex[6]*9.81*(h_c_ex+h_l_ex+h_sigma_ex)
#精馏段降液管液泛
H_d_distill=(h_c_distill+h_l_distill+h_sigma_distill)+(h_w+h_o)+0.153*u_o_prime_distill^2
phi_we_distill=H_d/(H_t+h_w)
#提馏段降液管液泛
H_d_ex=(h_c_ex+h_l_ex+h_sigma_ex)+(h_w+h_o)+0.153*u_o_prime_ex^2
phi_we_ex=H_d/(H_t+h_w)
#泛点率e
K=1
C_F_distill=0.11
C_F_ex=0.103
e_distill=prop_distill[8]*sqrt(prop_distill[5]/(prop_distill[6]-prop_distill[5]))/0.78/C_F_distill/(pi*0.8^2/4)
e_ex=prop_ex[8]*sqrt(prop_ex[5]/(prop_ex[6]-prop_ex[5]))/0.78/C_F_ex/(pi*0.8^2/4)
```


​    

```julia
#精馏段塔板负荷性能图
e_distill_max=0.8

#过量液沫夹带线m3/h
V_distill_max=0.78*C_F_distill*(pi*0.8^2/4)*e_distill_max/sqrt(prop_ex[5]/(prop_ex[6]-prop_ex[5]))*3600

#液相下限线m3/h
L_distill_min=0.006^1.5*L/(0.00284)^1.5

#严重漏液线m3/h
F_o_min=6
V_distill_min=(pi/4*0.8^2*0.109)*F_o_min/sqrt(prop_distill[5])*3600

#液相上限线m3/h
u_o_prime_distill_max=0.25
L_distill_max=u_o_prime_distill_max*L*h_o*3600

#降液管液泛线m3/h
H_d_distill_max=0.6*(0.6+h_w)
A_distill=5.34*prop_distill[5]/2/prop_distill[6]/9.81/(pi/4*0.8^2*0.109)^2/3600^2
B_distill=0.153/(L*h_o)^2/3600^2
C_distill=0.3840-0.03315-0.00075-0.0663

```



```julia
#提馏段塔板负荷性能图
e_ex_max=0.8

#过量液沫夹带线m3/h
V_ex_max=0.78*C_F_ex*(pi*0.8^2/4)*e_ex_max/sqrt(prop_ex[5]/(prop_ex[6]-prop_ex[5]))*3600

#液相下限线m3/h
L_ex_min=0.006^1.5*L/(0.00284)^1.5

#严重漏液线m3/h
F_o_min=6
V_ex_min=(pi/4*0.8^2*0.109)*F_o_min/sqrt(prop_ex[5])*3600

#液相上限线m3/h
u_o_prime_ex_max=0.25
L_ex_max=u_o_prime_ex_max*L*h_o*3600

#降液管液泛线m3/h
H_d_ex_max=0.6*(0.6+h_w)
A_ex=5.34*prop_ex[5]/2/prop_ex[6]/9.81/(pi/4*0.8^2*0.109)^2/3600^2
B_ex=0.153/(L*h_o)^2/3600^2
C_ex=0.3840-0.03315-0.00075-0.0663
```

```julia
#接管尺寸设计
#水力学核算后根据Aspen结果更新性质m3/s
prop_distill[8]=1.0922
prop_distill[9]=0.00222395
prop_ex[8]=0.99043
prop_ex[9]=0.00183092
prop_F[9]=0.0537359564733485
prop_B[9]=0.00183091582038441

#塔顶蒸汽出口管直径
u_distill_vapor=20
D_distill_vapor=sqrt(4*prop_distill[8]/pi/u_distill_vapor)

#回流液管管径
u_reflux=1
D_distill_reflux=sqrt(4*prop_distill[9]/pi/u_reflux)

#加料管径
u_feed=2
D_feed=sqrt(4*prop_F[9]/pi/u_feed)

#残液排出管径
u_bottom=0.8
D_bottom=sqrt(4*prop_B[9]/pi/u_feed)

#加热蒸汽管径
V_0_volumn=0.628904762274061
u_vapor=30
D_vapor=sqrt(4*V_0_volumn/pi/u_vapor)
```


