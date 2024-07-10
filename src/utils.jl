using Dates
import GridInterpolations as GI
import HDF5
using StaticArrays

function interpolate_value(ds1, ds2, field, t)
    df = dateformat"y-m-d_H:M:SZ"
    timestamp_ds1 = ds1["Times"]
    t_ds1 = DateTime(string(timestamp_ds1...),df)
    timestamp_ds2 = ds2["Times"]
    t_ds2 = DateTime(string(timestamp_ds2...),df)
    time_diff = (t_ds2 - t_ds1)*(0.001)/Millisecond(1)  #Gives time differewnce between the two dates in seconds
    interpolated_values = ( (ds2[field] - ds1[field])*(t/time_diff) ) + ds1[field]
    return interpolated_values
end

function interpolate_using_GI()

    #Demo code for interpolation. Poorly written. Modify it when dataset becomes clearer.
    # x_values = [i for i in 1:300]
    x_values = SVector{300,Int}(1:300)
    # y_values = [i for i in 1:300]
    y_values = SVector{300,Int}(1:300)
    # z_values = [i for i in 1:50]
    z_values = SVector{50,Int}(1:50)

    grid = GI.RectangleGrid(x_values,y_values,z_values)

    f1 = HDF5.h5open("./dataset/ENS_MEM_01/wrfout_d01_2020-04-24_21_00_00", "r")
    s = HDF5.read(f1)

    temp_values = s["T"]
    point = SVector(1.1,2.9,3.1)
    interpolated_temp_value = GI.interpolate(grid,temp_values,point)

end

function profile_interpolation(grid,value_array,point)
    interpolated_value = GI.interpolate(grid::GI.AbstractGrid,temp_values::DenseArray,point::SVector)
    return interpolated_value
end

function profile_interpolation2(grid,value_array,point)
    interpolated_value = GI.interpolate(grid,temp_values,point)
    return interpolated_value
end


x_ind, y_ind = 1,1
println( s["XLAT"][x_ind,y_ind,1],",", s["XLONG"][x_ind,y_ind,1] )

x_ind, y_ind = 300,1
println( s["XLAT"][x_ind,y_ind,1],",", s["XLONG"][x_ind,y_ind,1] )

x_ind, y_ind = 300,300
println( s["XLAT"][x_ind,y_ind,1],",", s["XLONG"][x_ind,y_ind,1] )

x_ind, y_ind = 1,300
println( s["XLAT"][x_ind,y_ind,1],",", s["XLONG"][x_ind,y_ind,1] )

h = ( s["PH"] .+ s["PHB"] )/9.81
z = (h[:,:,1:50,1] .+ h[:,:,2:51,1] )/2.0

#=
Keys in data

ACGRAUP
ACRUNOFF
ACSNOM
ACSNOW
ALBBCK
ALBEDO
BF
BH
BR
C1F
C1H
C2F
C2H
C3F
C3H
C4F
C4H
CANWAT
CF1
CF2
CF3
CFN
CFN1
CLAT
CLDFRA
COMPOSITE_REFL_10CM
COSALPHA
CUBOT
CUPPT
DIM0009
DIM0012
DN
DNW
DT100
DTBC
DTS
DTSEPS
DZS
DateStrLen
E
F
FCX
FLHC
FNDALBSI
FNDICEDEPTH
FNDSNOWH
FNDSNOWSI
FNDSOILW
FNM
FNP
GCX
GLW
GRAUPELNC
GRAUPELNCV
GRAUP_ACC_NC
GRDFLX
GRPL_MAX
GSW
HAILCAST_DIAM_MAX
HAILCAST_DIAM_MEAN
HAILCAST_DIAM_STD
HAILCAST_WDUR
HAILCAST_WUP_MASK
HAILDTACTTIME
HAIL_MAX2D
HAIL_MAXK1
HFX
HFX_FORCE
HFX_FORCE_TEND
HGT
H_DIABATIC
ISEEDARRAY_SPP_CONV
ISEEDARRAY_SPP_LSM
ISEEDARRAY_SPP_PBL
ISEEDARR_RAND_PERTURB
ISEEDARR_SKEBS
ISEEDARR_SPPT
ISLTYP
ITIMESTEP
IVGTYP
IWP
LAI
LAKEMASK
LAKE_DEPTH
LANDMASK
LANDUSEF
LAP_HGT
LH
LH_FORCE
LH_FORCE_TEND
LTG1_MAX
LTG2_MAX
LTG3_MAX
LU_INDEX
LWP
MAPFAC_M
MAPFAC_MX
MAPFAC_MY
MAPFAC_U
MAPFAC_UX
MAPFAC_UY
MAPFAC_V
MAPFAC_VX
MAPFAC_VY
MAVAIL
MAXCLDFRA
MAX_MSTFX
MAX_MSTFY
MF_VX_INV
MU
MU0
MUB
MUT
NCA_LTG
NCA_REFD
NCA_W
NCA_WQ
NCI_LTG
NCI_REFD
NCI_W
NCI_WQ
NEST_POS
OLR
P
P00
PB
PBLH
PC
PCB
PH
PHB
PRATEC
PREC_ACC_C
PREC_ACC_NC
PSFC
P_HYD
P_STRAT
P_TOP
Q2
QCG
QCLOUD
QFX
QGRAUP
QG_MAX_CI
QHAIL
QICE
QNCCN
QNDROP
QNGRAUPEL
QNHAIL
QNICE
QNRAIN
QNSNOW
QRAIN
QR_MAX_CI
QSNOW
QVAPOR
QVG
QVGRAUPEL
QVHAIL
QV_BASE
RAD_TTEN_DFI_1
RAD_TTEN_DFI_2
RAD_TTEN_DFI_3
RAD_TTEN_DFI_4
RAINC
RAINCV
RAINNC
RAINNCV
RDN
RDNW
RDX
RDY
REFD_MAX
REFL_10CM
REFL_10CM_1KM
REFL_10CM_4KM
REL_VORT
REL_VORT_MAX
RESM
RHOSNF
RLV
RLVN
RQVBLTEN
RQVFTEN
SAVE_TOPO_FROM_REAL
SEAICE
SFROFF
SH2O
SHDMAX
SHDMIN
SHEAR100
SINALPHA
SLOPECAT
SMOIS
SNOALB
SNOW
SNOWC
SNOWFALLAC
SNOWH
SNOWNC
SNOWNCV
SNOW_ACC_NC
SOILCAT
SOILCBOT
SOILCTOP
SOILHGT
SOILT1
SR
SST
SST_INPUT
STEP_NUMBER
SWDDIF
SWDDIFC
SWDDNI
SWDDNIC
SWDOWN
SWNORM
SWNORMMEAN
T
T00
T2
TH2
TISO
TKE
TLP
TLP_STRAT
TMN
TOPOSLPX
TOPOSLPY
TOPOSTDV
TSK
TSK_FORCE
TSK_FORCE_TEND
TSLB
TTEN_TIMES
T_BASE
T_INIT
Time
Times
U
U10
UDROFF
UH
UH16
UP_HELI_MAX
UP_HELI_MAX16
UST
U_BASE
U_FRAME
V
V10
VAR
VAR_SSO
VEGCAT
VEGFRA
VT_DBZ_WT
V_BASE
V_FRAME
W
WSPD10
WSPD10MAX
WSPD80
W_DN_MAX
W_MEAN
W_UP_MAX
XLAND
XLAT
XLAT_U
XLAT_V
XLONG
XLONG_U
XLONG_V
XTIME
ZETATOP
ZNT
ZNU
ZNW
ZOL
ZS
Z_BASE
bottom_top
bottom_top_stag
land_cat_stag
soil_cat_stag
soil_layers_stag
south_north
south_north_stag
west_east
west_east_stag




HGT - Terrain height
Pg 246 of this documentation - https://www2.mmm.ucar.edu/wrf/users/docs/user_guide_v4/v4.2/WRFUsersGuide_v42.pdf#page=24.00

=#