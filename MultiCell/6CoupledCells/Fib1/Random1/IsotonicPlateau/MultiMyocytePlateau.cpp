#include <iostream>
#include <boost/array.hpp>
#include <fstream>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

// dimension of problem
const int cells = 6;
const int segments = 50;
const double segments_double = 50.0;
const int StateNum = 57;
const int n = StateNum*segments*cells;
typedef boost::array< double , n > state_type;

// define structure
struct OpT
{
    std::array<double, cells*segments> SLsetAll;
    std::array<double, cells*segments> GapAll;
};

//output
double kPrint;

//internal shortening
double sumSL[cells] = {0.0};//, 0.0, 0.0};//, 0.0, 0.0};
const double pi = 3.14159265358979323846264338328;

// constants
const double xOpti[8] = {0.9572, 0.8075, 1.4267, 1.2535, 0.8565, 0.7009, 0.1260, 4.2595};
const int gto_num     = 5;
const int gca_num     = 3;
const int gkr_num     = 3;
const int gks_num     = 2;
const int gkix_num    = 1;
const int gnak_num    = 4;

const int f_inhib_num = 2;
const int K_o_num     = 2;
const int f_katp_num  = 1;
const int na_num      = 1;
const double naRange[2]  = {0.0, 0.3};

const double Bmax_Calsequestrin          = 0.14;         // millimolar (in Ca_buffer)
const double Bmax_SLB_SL                 = 0.0374;       // millimolar (in Ca_buffer)
const double Bmax_SLB_jct                = 0.0046;       // millimolar (in Ca_buffer)
const double Bmax_SLHigh_SL              = 0.0134;       // millimolar (in Ca_buffer)
const double Bmax_SLHigh_jct             = 0.00165;      // millimolar (in Ca_buffer)
const double koff_Calsequestrin          = 65.0;         // per_millisecond (in Ca_buffer)
const double koff_SLB                    = 1.3;          // per_millisecond (in Ca_buffer)
const double koff_SLHigh                 = 30.0e-3;      // per_millisecond (in Ca_buffer)
const double kon_Calsequestrin           = 100.0;        // per_millimolar_per_millisecond (in Ca_buffer)
const double kon_SL                      = 100.0;        // per_millimolar_per_millisecond (in Ca_buffer)

const double Fx_ICaL_SL                  = 0.1;          // dimensionless (in ICaL)
const double Fx_ICaL_jct                 = 0.9;          // dimensionless (in ICaL)
const double PCa                         = 5.4e-4;       // litre_per_farad_millisecond (in ICaL)
const double PK                          = 2.7e-7;       // litre_per_farad_millisecond (in ICaL)
const double PNa                         = 1.5e-8;       // litre_per_farad_millisecond (in ICaL)
const double Q10_CaL                     = 1.8;          // dimensionless (in ICaL)
const double gamma_Cai                   = 0.341;        // dimensionless (in ICaL)
const double gamma_Cao                   = 0.341;        // dimensionless (in ICaL)
const double gamma_Ki                    = 0.75;         // dimensionless (in ICaL)
const double gamma_Ko                    = 0.75;         // dimensionless (in ICaL)
const double gamma_Nai                   = 0.75;         // dimensionless (in ICaL)
const double gamma_Nao                   = 0.75;         // dimensionless (in ICaL)

const double Fx_CaBk_SL                  = 0.89;         // dimensionless (in ICab)
const double Fx_CaBk_jct                 = 0.11;         // dimensionless (in ICab)
const double G_CaBk                      = 0.0002513;    // milliS_per_microF (in ICab)

const double Fx_SLCaP_SL                 = 0.89;         // dimensionless (in ICap)
const double Fx_SLCaP_jct                = 0.11;         // dimensionless (in ICap)
const double H1                          = 1.6;          // dimensionless (H in ICap)
const double Km                          = 0.0005;       // millimolar (in ICap)
const double Q10_SLCaP                   = 2.35;         // dimensionless (in ICap)
const double V_maxAF                     = 0.0673;       // microA_per_microF (in ICap)

const double Fx_Cl_SL                    = 0.89;         // dimensionless (in ICl_Ca)
const double Fx_Cl_jct                   = 0.11;         // dimensionless (in ICl_Ca)
const double G_Cl                        = 0.109625;     // milliS_per_microF (in ICl_Ca)
const double Kd_ClCa                     = 0.1;          // millimolar (in ICl_Ca)

const double G_ClBk                      = 0.009;        // milliS_per_microF (in IClb)

const double g_Kp                        = 0.001;        // milliS_per_microF (in IKp)

const double Fx_Ks_SL                    = 0.89;         // dimensionless (in IKs)
const double Fx_Ks_jct                   = 0.11;         // dimensionless (in IKs)
const double pKNa                        = 0.01833;      // dimensionless (in IKs)

const double Fx_NCX_SL                   = 0.89;         // dimensionless (in INaCa)
const double Fx_NCX_jct                  = 0.11;         // dimensionless (in INaCa)
const double HNa                         = 3.0;          // dimensionless (in INaCa)
const double K_mCai                      = 0.00359;      // millimolar (in INaCa)
const double K_mCao                      = 1.3;          // millimolar (in INaCa)
const double K_mNai                      = 12.29;        // millimolar (in INaCa)
const double K_mNao                      = 87.5;         // millimolar (in INaCa)
const double Kd_act                      = 0.000256;     // millimolar (in INaCa)
const double Q10_NCX                     = 1.57;         // dimensionless (in INaCa)
const double V_max_1                     = 9.0;          // microA_per_microF (V_max in INaCa)
const double eta                         = 0.35;         // dimensionless (in INaCa)
const double ksat                        = 0.27;         // dimensionless (in INaCa)

const double Fx_Na_SL                    = 0.89;         // dimensionless (in INa)
const double Fx_Na_jct                   = 0.11;         // dimensionless (in INa)
const double G_INa                       = 16.0;         // milliS_per_microF (in INa)

const double Fx_NaK_SL                   = 0.89;         // dimensionless (in INaK)
const double Fx_NaK_jct                  = 0.11;         // dimensionless (in INaK)
const double H_NaK                       = 4.0;          // dimensionless (in INaK)
const double I_NaK_max                   = 1.90719;      // microA_per_microF (in INaK)
const double Km_Ko                       = 1.5;          // millimolar (in INaK)
const double Km_Nai                      = 11.0;         // millimolar (in INaK)

const double Fx_NaBk_SL                  = 0.89;         // dimensionless (in INab)
const double Fx_NaBk_jct                 = 0.11;         // dimensionless (in INab)
const double G_NaBk                      = 0.297e-3;     // milliS_per_microF (in INab)

const double G_tof                       = 0.02;         // milliS_per_microF (in Itof)

const double G_tos                       = 0.06;         // milliS_per_microF (in Itos)

const double KSRleak                     = 5.348e-6;     // per_millisecond (in Jleak_SR)

const double H2                          = 1.787;        // dimensionless (H in Jpump_SR)
const double Kmf                         = 0.000246;     // millimolar (in Jpump_SR)
const double Kmr                         = 1.7;          // millimolar (in Jpump_SR)
const double Q10_SRCaP                   = 2.6;          // dimensionless (in Jpump_SR)
const double V_max_2                     = 5.3114e-3;    // millimolar_per_millisecond (V_max in Jpump_SR)

const double EC50_SR                     = 0.45;         // millimolar (in Jrel_SR)
const double HSR                         = 2.5;          // dimensionless (in Jrel_SR)
const double Max_SR                      = 15.0;         // dimensionless (in Jrel_SR)
const double Min_SR                      = 1.0;          // dimensionless (in Jrel_SR)
const double kiCa                        = 0.5;          // per_millimolar_per_millisecond (in Jrel_SR)
const double koCa                        = 10.0;         // per_millimolar2_per_millisecond (in Jrel_SR)
const double kom                         = 0.06;         // per_millisecond (in Jrel_SR)
const double kim                         = 0.005;        // per_millisecond (in Jrel_SR)
const double ks                          = 25.0;         // per_millisecond (in Jrel_SR)

const double Bmax_SL                     = 1.65;         // millimolar (in Na_buffer)
const double Bmax_jct                    = 7.561;        // millimolar (in Na_buffer)
const double koff_EP                     = 1.0e-3;       // per_millisecond (in Na_buffer)
const double kon_EP                      = 0.0001;       // per_millimolar_per_millisecond (in Na_buffer)

const double offset                      = 10.1;		     // Offset of simulation from start
const double stim_amplitude              = 20.0*segments;// microA_per_microF (in cell)
const double stim_duration               = 3.0;          // millisecond (in cell)
const double stim_period                 = 400.0;        // millisecond (in cell)
const double stim_start                  = offset;       // millisecond (in cell)

const double Bmax_Calmodulin             = 0.024;        // millimolar (in cytosolic_Ca_buffer)
const double Bmax_Myosin_Ca              = 0.14;         // millimolar (in cytosolic_Ca_buffer)
const double Bmax_Myosin_Mg              = 0.14;         // millimolar (in cytosolic_Ca_buffer)
const double Bmax_SRB                    = 0.0171;       // millimolar (in cytosolic_Ca_buffer)
const double Bmax_TroponinC              = 0.07;         // millimolar (in cytosolic_Ca_buffer)
const double Bmax_TroponinC_Ca_Mg_Ca     = 0.14;         // millimolar (in cytosolic_Ca_buffer)
const double Bmax_TroponinC_Ca_Mg_Mg     = 0.14;         // millimolar (in cytosolic_Ca_buffer)
const double koff_Calmodulin             = 238.0e-3;     // per_millisecond (in cytosolic_Ca_buffer)
const double koff_Myosin_Ca              = 0.46e-3;      // per_millisecond (in cytosolic_Ca_buffer)
const double koff_Myosin_Mg              = 0.057e-3;     // per_millisecond (in cytosolic_Ca_buffer)
const double koff_SRB                    = 60.0e-3;      // per_millisecond (in cytosolic_Ca_buffer)
const double koff_TroponinC              = 19.6e-3;      // per_millisecond (in cytosolic_Ca_buffer)
const double koff_TroponinC_Ca_Mg_Ca     = 0.032e-3;     // per_millisecond (in cytosolic_Ca_buffer)
const double koff_TroponinC_Ca_Mg_Mg     = 3.33e-3;      // per_millisecond (in cytosolic_Ca_buffer)
const double kon_Calmodulin              = 34.0;         // per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
const double kon_Myosin_Ca               = 13.8;         // per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
const double kon_Myosin_Mg               = 15.7e-3;      // per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
const double kon_SRB                     = 100.0;        // per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
const double kon_TroponinC               = 32.7;         // per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
const double kon_TroponinC_Ca_Mg_Ca      = 2.37;         // per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
const double kon_TroponinC_Ca_Mg_Mg      = 3.0e-3;       // per_millimolar_per_millisecond (in cytosolic_Ca_buffer)

const double Cao                         = 4.2;          // millimolar (in model_parameters)
const double Ki                          = 135.0;        // millimolar (in model_parameters)
const double Cli                         = 15.0;         // millimolar (in model_parameters)
const double Clo                         = 150.0;        // millimolar (in model_parameters)
const double Cm                          = 1.381e-10/segments_double;    // farad (in model_parameters)
const double F                           = 96485.0;      // coulomb_per_mole (in model_parameters)
const double Mgi                         = 1.0;          // millimolar (in model_parameters)
const double Nao                         = 140.0;        // millimolar (in model_parameters)
const double R2                          = 8314.3;       // joule_per_kilomole_kelvin (R in model_parameters)
const double T                           = 310.0;        // kelvin (in model_parameters)


//// ---------------------------------------------------------------------------
//  Variable inputs
//  ---------------------------------------------------------------------------

const double  rangeGto[5]    = {0.7, 0.85, 1.0, 1.15, 1.3};
const double  rangeGCaL[5]   = {0.7, 0.85, 1.0, 1.15, 1.3};
const double  rangeGKr[5]    = {0.7, 0.85, 1.0, 1.15, 1.3};
const double  rangeGKs[5]    = {0.7, 0.85, 1.0, 1.15, 1.3};
const double  rangeGK1[5]    = {0.7, 0.85, 1.0, 1.15, 1.3};
const double  rangeGNaK[5]   = {0.7, 0.85, 1.0, 1.15, 1.3};

//// ---------------------------------------------------------------------------
//  Ischaemia variables
//  ---------------------------------------------------------------------------

const double  K_o_normal         = 5.4;
const double  g_katp             = 2.61;
const double  f_inhibRange[2]    = {0.75, 1.0};
const double  K_oRange[7]        = {4.0, 5.4, 7.0, 9.0, 12.0, 15.0, 17.0};
const double  f_katpRange[4]     = {0.0, 0.006, 0.008, 0.01};

//// ---------------------------------------------------------------------------
//   Input Nimrod variable values
//  ---------------------------------------------------------------------------
const double ScaleFactorGto  = rangeGto[gto_num-1];
const double ScaleFactorGCaL = rangeGCaL[gca_num-1];
const double ScaleFactorGKr  = rangeGKr[gkr_num-1];
const double ScaleFactorGKs  = rangeGKs[gks_num-1];
const double ScaleFactorGK1  = rangeGK1[gkix_num-1];
const double ScaleFactorGNaK = rangeGNaK[gnak_num-1]*(1-naRange[na_num-1]);

const double f_inhib         = f_inhibRange[f_inhib_num-1];
const double K_o             = K_oRange[K_o_num-1];            // millimolar [in model_parameters)
const double f_katp          = f_katpRange[f_katp_num-1];

//// ---------------------------------------------------------------------------
//  Computed variables
//  ---------------------------------------------------------------------------
//Cell geometry
const double cell_length     = 100.0/segments_double;        // micrometre (in model_parameters)
const double cell_radius     = 10.25;        // micrometre (in model_parameters)
const double Vol_Cell        = 3.141592654*pow(cell_radius/1000.0, 2.0)*cell_length/pow(1000.0, 3.0);
const double Vol_myo         = 0.65*Vol_Cell;
const double Vol_SR          = 0.035*Vol_Cell;
const double Vol_SL          = 0.02*Vol_Cell;
const double Vol_jct         = 0.0539*0.01*Vol_Cell;

// Temperature Dependence
const double Q_CaL       = pow(Q10_CaL, (T-310.0)/10.0);
const double Q_NCX       = pow(Q10_NCX, (T-310.0)/10.0);
const double Q_SLCaP     = pow(Q10_SLCaP, (T-310.0)/10.0);
const double Q_SRCaP     = pow(Q10_SRCaP, (T-310.0)/10.0);


// ---------------------------------------------------------------------------
//  Initialize state variables
// ---------------------------------------------------------------------------
double Vm, Ca_Calsequestrin, Ca_SL, Ca_SLB_SL, Ca_SLB_jct,
	Ca_SLHigh_SL, Ca_SLHigh_jct, Ca_SR, Ca_jct, Cai, d, fCaB_SL, fCaB_jct,
	fEP, Xr, Xs, h, j, m, X_tof, Y_tof, R_tof, X_tos, Y_tos, I, O, R1,
	Na_SL, Na_SL_buf, Na_jct, Na_jct_buf, Nai, Ca_Calmodulin, Ca_Myosin,
	Ca_SRB, Ca_TroponinC, Ca_TroponinC_Ca_Mg, Mg_Myosin, Mg_TroponinC_Ca_Mg,
  N_NoXB, P_NoXB, N, P, XBprer, XBpostr, SL, xXBpostr, xXBprer,
	TRPNCaL, TRPNCaH, force, intf,
  Vf, r_KV, s_KV, K_fi, Na_fi;

double dVm, dCalsequestrin, dCa_SL, dCa_SLB_SL, dCa_SLB_jct,
  dCa_SLHigh_SL, dCa_SLHigh_jct, dCa_SR, dCa_jct, dCai, dd, dfCaB_SL, dfCaB_jct,
  dfEP, dXr, dXs, dh, dj, dm, dX_tof, dY_tof, dR_tof, dX_tos, dY_tos, dI, dO, dR1,
  dNa_SL, dNa_SL_buf, dNa_jct, dNa_jct_buf, dNai, dCa_Calmodulin, dCa_Myosin,
  dCa_SRB, dCa_TroponinC, dCa_TroponinC_Ca_Mg, dMg_Myosin, dMg_TroponinC_Ca_Mg,
  dN_NoXB, dP_NoXB, dN, dP, dXBprer, dXBpostr, dSL, dxXBpostr, dxXBprer,
	dTRPNCaL, dTRPNCaH, dforce, dintf,
  dVf, dr_KV, ds_KV, dK_fi, dNa_fi;

double dSL1[segments*cells];

// ---------------------------------------------------------------------------
//  Initialize state variables
// ---------------------------------------------------------------------------
// Sarcomere Geometry
const double len_thin  = 1.2;      // (um) thin filament length
const double len_thick = 1.65;     // (um) thick filament length
const double len_hbare = 0.1;      // (um) length of bare portion of thick filament

// Temperature Dependence
const double Qkon = 1.5;
const double Qkoff = 1.3;
const double Qkn_p = 1.6;
const double Qkp_n = 1.6;
const double Qfapp = 6.25;
const double Qgapp = 2.5;
const double Qhf = 6.25;
const double Qhb = 6.25;
const double Qgxb = 6.25;

// Ca binding to troponin
const double kon     = 50.0e-3*xOpti[0];      // (1/[ms uM])
const double koffL   = 250.0e-3*xOpti[1];     // (1/ms)
const double koffH   = 25.0e-3;      // (1/ms)
const double perm50  = 0.45;       // perm variable that controls n to p transition
const double nperm   = 10.0;         //   in Hill-like fashion
const double kn_p    = 500.0e-3;     // (1/ms)
const double kp_n    = 50.0e-3;      // (1/ms)
const double koffmod = 0.9;        // mod to change species

// Thin filament regulation and crossbridge cycling
const double fapp    = 500.0e-3;     // (1/ms) XB on rate
const double gapp    = 70.0e-3;      // (1/ms) XB off rate
const double gSLmod  = 6.0;          // controls SL effect on gapp
const double hfXB    = 2000.0e-3;    // (1/ms) rate between pre-force and force states

//
const double hfmdc   = 5.0;          //
const double hb      = 400.0e-3;     // (1/ms) rate between pre-force and force states
const double hbmdc   = 0.0;          //
const double gxb     = 70.0e-3;      // (1/ms) ATP consuming transition rate
const double xbmodsp = 0.2;        // rabbit specific modification for XB cycling rates

// Mean strain of strongly-bound states
const double x_0     = 0.007;      // (um) strain induced by head rotation
const double xPsi    = 2.0;          // scaling factor balancing SL motion and XB cycling

// Normalized active and passive force
const double SLrest  = 1.89;       // (um) rest SL length for 0 passive force
const double PCon_t  = 0.002;      // (norm force) passive force due to titin
const double PExp_t  = 10.0;         //   these apply to trabeculae and single cells only
const double SL_c    = 2.25;       // (um) resting length for collagen
const double PCon_c  = 0.02;       // (norm force) passive force due to collagen
const double PExp_c  = 70.0;         //   these apply to trabeculae and single cells only

// Calculation of complete muscle response
const double massf   = 0.00025e6;  // ([norm force ms.^2]/um) muscle mass
const double visc    = 0.003e3;    // ([norm force ms]/um) muscle viscosity
const double KSE     = 1.0;          // (norm force/um) series elastic element
const double kxb     = 120.0;        // (mN/mm.^2) maximal force
const double Trop_conc = 70.0;       // (uM) troponin concentration

double SLset;

double CaiMech;
const double Temp      = 310.0;            // Temperature (K)

// Compute single overlap fractions
double sovr_ze, sovr_cle, len_sovr,
  SOVFThick, SOVFThin,
  Tropreg, permtot, inprmt,
  konT, hfmd, hbmd,
  kn_pT, kp_nT, fappT, gapSLmd, gappT,
  hfT, hbT, koffLT, koffHT,
  SSXBprer, SSXBpostr, Fnordv,
  active, ppforce, ppforce_t,
  ppforce_c, preload, afterload,
  lmbda, Cp, b_ff, passive,
  total_force, C2, GaussCurve,
  sigman, sigmap, gxbmd, gxbT,
  dutyprer, dutypostr,
  FrSBXB, dFrSBXB, TropTot,
  signxXBpostr, signxXBprer, signSL, signSLset,
  heavisidePostR, dsovr_ze, dsovr_cle,
  dlen_sovr, dSOVFThin, dSOVFThick,
  dTropTot, dactive, dppforce_t,
  dppforce_c, dppforce, dsfib;

// interaction between agjacent SR, cytosol and SL due to fluxes (based on Voigt et al.)
double J_diff_cyt, J_diff_SR, J_diff_SL;
const double tau_seg_cyt = 0.6*2.0;
const double tau_seg_SR = 15.0*10.0;
const double tau_seg_SL = 3.4*2/15.0;


// diffusion constants
const double Diff_Na_jct_SL = 1.8313e-14/segments_double;
const double Diff_Na_SL_myo = 1.6386e-12/segments_double;
const double Diff_Ca_jct_SL = 8.2413e-13/segments_double;
const double Diff_Ca_SL_myo = 3.7243e-12/segments_double;

double alpha_m, beta_m, alpha_h, beta_h,
alpha_j, beta_j, openProb, i_Na_jct,
i_Na_SL, i_Na, i_Nab_jct, i_Nab_SL,
i_Nab, f_NaK, i_NaK_jct, i_NaK_SL,
i_NaK,
Rr, i_Kr, Xr_infinity, tau_Xr,
Xs_infinity, tau_Xs, i_Ks_jct,
i_Ks_SL, i_Kp, i_Ks,
i_tof, i_tos, i_to, X_tof_infinity,
tau_X_tof, Y_tof_infinity,tau_Y_tof,
R_tof_infinity, tau_R_tof, X_tos_infinity,
tau_X_tos, Y_tos_infinity, tau_Y_tos,
alpha_K1, beta_K1, K1_infinity, i_K1,
i_Cl_Ca, i_Clb, d_infinity, tau_d,
f_infinity, tau_f, temp, fCa_jct,
fCa_SL, i_CaL_Ca_SL, i_CaL_Ca_jct,
i_CaL_Na_jct, i_CaL_Na_SL, i_CaL_K, i_CaL,
Ka_jct, temp_jct, i_NaCa_jct, Ka_SL,
temp_SL, i_NaCa_SL, i_NaCa,j_pump_SR,
i_Cap_jct, i_Cap_SL, i_Cap,
i_Cab_jct, i_Cab_SL, i_Cab,
kCaSR, koSRCa, kiSRCa, RI,
i_Ca_SL_tot, j_leak_SR, j_rel_SR,
J_Na_jct_SL, J_Na_SL_myo,
i_Ca_jct_tot, J_Ca_jct_SL, J_Ca_SL_myo, past, i_Stim,
a_katp, i_katp, J_diff_Vm, dCa_SL_tot_bound,
dCa_jct_tot_bound, dCa_cytosol_tot_bound;

// membrane potential
const double tau_Vm = .0001e1;


const double G_gap = 4.0e2;


// Nernst potential
double E_Na_SL, E_Na_jct;
const double E_Cl        = R2*T/F*log(Cli/Clo);
double E_Ca_jct;
const double E_K         = R2*T/F*log(K_o/Ki);
double E_Ks;
double E_Ca_SL;

// Conductances
double pCa_jct;
double pCa_SL;
const double G_K1        = 0.9*sqrt(K_o/5.4);
const double G_IKr       = 0.03*sqrt(K_o/5.4);
double G_Ks_jct;
double G_Ks_SL;

const double sigma       = (exp(Nao/67.3)-1.0)/7.0;
const double sigma2      = 20.0;
const double Gauss       = 60.0;
const double bcl         = 400.0;
const double GaussFactor = 3.8;


const double K_fo       = 5.3581;
const double Na_fo      = 130.0110;
const double Rf         = R2;
const double Tf         = T;
const double Ff         = F;

const double g_KV       = 0.25; // [milli_S./micro_F]
const double g_fK1      = 0.4822; //[milli_S./micro_F]
const double I_fNaK_inf = 1.644; //[micro_A./micro_F]
const double B_fNaK     = -200.0;  //[mV]
const double V_fNaK_rev = -150.0;  //[mV]
const double K_fmK      = 1.0;   //[mmol./L]
const double K_fmNa     = 11.0;  //[mmol./L]
const double G_fbNa     = 0.0095; //[milli_S./micro_F]
const double Cmf        = 1.0/6.3; //[F]
const double Vol_f      = 0.00137; //[microL]


double E_fK, E_fNa, tau_r, tau_s, r_inf, s_inf, I_KV,
    alpha_fK1, beta_fK1, I_fK1, I_fNaK, I_fbNa, I_f,
    E_NaGap, E_KGap, GapFib, I_gapNa, I_gapK, I_gap;
// ---------------------------------------------------------------------------
//  Run Algorithm
// ---------------------------------------------------------------------------

class Shannon
{
    struct OpT T1;

public:
    Shannon(struct OpT G) : T1(G) {}

  void operator()( const state_type &x , state_type &dxdt , double t)
{
    for (int jCount = 0; jCount<cells; jCount++){
      for (int i = 0; i < segments; i++){
          // input Mechanics model
          N_NoXB    = x[39 + StateNum*(i) + StateNum*(segments*jCount)];
          P_NoXB    = x[40 + StateNum*(i) + StateNum*(segments*jCount)];
          N         = x[41 + StateNum*(i) + StateNum*(segments*jCount)];
          P         = x[42 + StateNum*(i) + StateNum*(segments*jCount)];
          XBprer    = x[43 + StateNum*(i) + StateNum*(segments*jCount)];
          XBpostr   = x[44 + StateNum*(i) + StateNum*(segments*jCount)];
          SL        = x[45 + StateNum*(i) + StateNum*(segments*jCount)];
          xXBpostr  = x[46 + StateNum*(i) + StateNum*(segments*jCount)];
          xXBprer   = x[47 + StateNum*(i) + StateNum*(segments*jCount)];
          TRPNCaL   = x[48 + StateNum*(i) + StateNum*(segments*jCount)];
          TRPNCaH   = x[49 + StateNum*(i) + StateNum*(segments*jCount)];
          force     = x[50 + StateNum*(i) + StateNum*(segments*jCount)];
          intf      = x[51 + StateNum*(i) + StateNum*(segments*jCount)];

          // ---------------------------------------------------------------------------
        	//  Mechanics model
        	//  ---------------------------------------------------------------------------
          // set SLset
          SLset = T1.SLsetAll[i + segments*jCount];

        	// Time-Varying Parameters
        	CaiMech   = Cai*1000.0;        //cytosolic calcium concentration from EP model

          // Compute single overlap fractions
      		sovr_cle = max(SL/2.0-(SL-len_thin), len_hbare/2.0);

      		if (len_thick/2.0 < SL/2.0){
      			sovr_ze = len_thick/2.0; // z-line end
      		} else{
      			sovr_ze = SL/2.0;
      		}

      		len_sovr  = sovr_ze-sovr_cle;               // single overlap length
      		SOVFThick = len_sovr*2.0/(len_thick-len_hbare);  // thick filament overlap frac
      		SOVFThin  = len_sovr/len_thin;                 // thin filament overlap frac

      		// Compute combined Ca binding to high- (w/XB) and low- (no XB) sites
      		Tropreg = (1.0-SOVFThin)*TRPNCaL + SOVFThin*TRPNCaH;
      		permtot = sqrt(1.0/(1.0+pow(perm50/Tropreg, nperm)));
      		if (1.0/permtot < 100.0){
      			inprmt  = 1.0/permtot;
      		} else{
      			inprmt  = 100.0;
      		}


      		// Adjustments for Ca activation, temperature, SL, stress and strain
      		konT       = kon*pow(Qkon, ((Temp-310.0)/10.0));

      		if (xXBprer<0.0){
      			signxXBprer = -1.0;
      		}else{
      			signxXBprer = 1.0;
      		}
      		if ((xXBpostr-x_0)<0.0){
      			signxXBpostr = -1.0;
      		}else{
      			signxXBpostr = 1.0;
      		}

      		hfmd    = exp(-signxXBprer*hfmdc*pow((xXBprer/x_0), 2.0));
      		hbmd    = exp(signxXBpostr*hbmdc*pow((xXBpostr-x_0)/x_0, 2.0));

      		// Adjustments for Ca activation, temperature, SL, stress and strain
      		kn_pT   = kn_p*permtot*pow(Qkn_p, ((Temp-310.0)/10.0));
      		kp_nT   = kp_n*inprmt*pow(Qkp_n, ((Temp-310.0)/10.0));
      		fappT   = fapp*xbmodsp*pow(Qfapp, ((Temp-310.0)/10.0));
      		gapSLmd = 1.0 + (1.0-SOVFThick)*gSLmod;
      		gappT   = gapp*gapSLmd*xbmodsp*pow(Qgapp, ((Temp-310.0)/10.0));
      		hfT     = hfXB*hfmd*xbmodsp*pow(Qhf, ((Temp-310.0)/10.0));
      		hbT     = hb*hbmd*xbmodsp*pow(Qhb, ((Temp-310.0)/10.0));
      		koffLT  = koffL*pow(Qkoff, ((Temp-310.0)/10.0))*koffmod;
      		koffHT  = koffH*pow(Qkoff, ((Temp-310.0)/10.0))*koffmod;

      		// steady-state fractions in XBprer and XBpostr using King-Altman rule
      		SSXBprer = (hb*fapp+gxb*fapp)/
      		  (gxb*hfXB+fapp*hfXB+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);
      		SSXBpostr = fapp*hfXB/(gxb*hfXB+fapp*hfXB+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);

      		// normalization for scaling active and passive force (maximal force)
      		Fnordv = kxb*x_0*SSXBpostr;

      		// Calculate forces (active, passive, preload, afterload)
      		force = kxb*SOVFThick*(xXBpostr*XBpostr+xXBprer*XBprer);
      		active = force/Fnordv;

      		if ((SL-SLrest)<0.0){
      			signSL = -1.0;
      		}else{
      			signSL = 1.0;
      		}
      		ppforce_t = signSL*PCon_t*(exp(PExp_t*abs(SL-SLrest))-1.0);
      		if ((SL-SL_c) < 0.0){
      		    ppforce_c = 0.0*PCon_c*(exp(PExp_c*abs(SL-SL_c))-1.0);
      		} else if ((SL-SL_c) == 0.0){
      		    ppforce_c = 1.0/2*PCon_c*(exp(PExp_c*abs(SL-SL_c))-1.0);
      		} else {
      		    ppforce_c = 1.0*PCon_c*(exp(PExp_c*abs(SL-SL_c))-1.0);
      		}
      		ppforce = ppforce_t + ppforce_c;
      		if ((SLset-SLrest)<0.0){
      			int signSLset = -1.0;
      		}else{
      			int signSLset = 1.0;
      		}
      		preload = signSLset*PCon_t*(exp(PExp_t*abs(SLset-SLrest))-1.0);
      		afterload = 0.0;

      		// change in SL and forces
      		lmbda = SL/SLset;
          lmbda = (SL/SLset)*(lmbda>1.0);
          lmbda = lmbda+(lmbda==0.0);


      		Cp = 0.2;
      		b_ff = 50.0;
      		passive = Cp*b_ff*(-0.5 + 0.5*(lmbda*lmbda))*exp(b_ff*((-0.5 +0.5*(lmbda*lmbda))*(-0.5 + 0.5*(lmbda*lmbda))));
      		total_force = -preload + active + passive;
      		dintf = (-ppforce+preload-active+afterload);

          dSL = (-total_force + afterload)/visc;

          lmbda = (SL/SLset)*(lmbda>1.0);
          lmbda = lmbda+(lmbda==0.0);

      		dSL1[i+segments*jCount] = dSL;
        }
      }

      for (int jCount = 0; jCount<cells; jCount++){
        for (int i = 0; i < segments; i++){
        //input EP model
        Vm                  = x[ 0 + StateNum*i + StateNum*(segments*jCount)];
        Ca_Calsequestrin    = x[ 1 + StateNum*i + StateNum*(segments*jCount)];
        Ca_SL               = x[ 2 + StateNum*i + StateNum*(segments*jCount)];
        Ca_SLB_SL           = x[ 3 + StateNum*i + StateNum*(segments*jCount)];
        Ca_SLB_jct          = x[ 4 + StateNum*i + StateNum*(segments*jCount)];
        Ca_SLHigh_SL        = x[ 5 + StateNum*i + StateNum*(segments*jCount)];
        Ca_SLHigh_jct       = x[ 6 + StateNum*i + StateNum*(segments*jCount)];
        Ca_SR               = x[ 7 + StateNum*i + StateNum*(segments*jCount)];
        Ca_jct              = x[ 8 + StateNum*i + StateNum*(segments*jCount)];
        Cai                 = x[ 9 + StateNum*i + StateNum*(segments*jCount)];
        d                   = x[10 + StateNum*i + StateNum*(segments*jCount)];
        fCaB_SL             = x[11 + StateNum*i + StateNum*(segments*jCount)];
        fCaB_jct            = x[12 + StateNum*i + StateNum*(segments*jCount)];
        fEP                 = x[13 + StateNum*i + StateNum*(segments*jCount)];
        Xr                  = x[14 + StateNum*i + StateNum*(segments*jCount)];
        Xs                  = x[15 + StateNum*i + StateNum*(segments*jCount)];
        h                   = x[16 + StateNum*i + StateNum*(segments*jCount)];
        j                   = x[17 + StateNum*i + StateNum*(segments*jCount)];
        m                   = x[18 + StateNum*i + StateNum*(segments*jCount)];
        X_tof               = x[19 + StateNum*i + StateNum*(segments*jCount)];
        Y_tof               = x[20 + StateNum*i + StateNum*(segments*jCount)];
        R_tof               = x[21 + StateNum*i + StateNum*(segments*jCount)];
        X_tos               = x[22 + StateNum*i + StateNum*(segments*jCount)];
        Y_tos               = x[23 + StateNum*i + StateNum*(segments*jCount)];
        I                   = x[24 + StateNum*i + StateNum*(segments*jCount)];
        O                   = x[25 + StateNum*i + StateNum*(segments*jCount)];
        R1                  = x[26 + StateNum*i + StateNum*(segments*jCount)];
        Na_SL               = x[27 + StateNum*i + StateNum*(segments*jCount)];
        Na_SL_buf           = x[28 + StateNum*i + StateNum*(segments*jCount)];
        Na_jct              = x[29 + StateNum*i + StateNum*(segments*jCount)];
        Na_jct_buf          = x[30 + StateNum*i + StateNum*(segments*jCount)];
        Nai                 = x[31 + StateNum*i + StateNum*(segments*jCount)];
        Ca_Calmodulin       = x[32 + StateNum*i + StateNum*(segments*jCount)];
        Ca_Myosin           = x[33 + StateNum*i + StateNum*(segments*jCount)];
        Ca_SRB              = x[34 + StateNum*i + StateNum*(segments*jCount)];
        Ca_TroponinC        = x[35 + StateNum*i + StateNum*(segments*jCount)];
        Ca_TroponinC_Ca_Mg  = x[36 + StateNum*i + StateNum*(segments*jCount)];
        Mg_Myosin           = x[37 + StateNum*i + StateNum*(segments*jCount)];
        Mg_TroponinC_Ca_Mg  = x[38 + StateNum*i + StateNum*(segments*jCount)];

        // input Mechanics model
        N_NoXB              = x[39 + StateNum*i + StateNum*(segments*jCount)];
        P_NoXB              = x[40 + StateNum*i + StateNum*(segments*jCount)];
        N                   = x[41 + StateNum*i + StateNum*(segments*jCount)];
        P                   = x[42 + StateNum*i + StateNum*(segments*jCount)];
        XBprer              = x[43 + StateNum*i + StateNum*(segments*jCount)];
        XBpostr             = x[44 + StateNum*i + StateNum*(segments*jCount)];
        SL                  = x[45 + StateNum*i + StateNum*(segments*jCount)];
        xXBpostr            = x[46 + StateNum*i + StateNum*(segments*jCount)];
        xXBprer             = x[47 + StateNum*i + StateNum*(segments*jCount)];
        TRPNCaL             = x[48 + StateNum*i + StateNum*(segments*jCount)];
        TRPNCaH             = x[49 + StateNum*i + StateNum*(segments*jCount)];
        force               = x[50 + StateNum*i + StateNum*(segments*jCount)];
        intf                = x[51 + StateNum*i + StateNum*(segments*jCount)];

        // input fibroblast
        Vf                  = x[52 + StateNum*i + StateNum*(segments*jCount)];
        r_KV                = x[53 + StateNum*i + StateNum*(segments*jCount)];
        s_KV                = x[54 + StateNum*i + StateNum*(segments*jCount)];
        K_fi                = x[55 + StateNum*i + StateNum*(segments*jCount)];
        Na_fi               = x[56 + StateNum*i + StateNum*(segments*jCount)];


        // ---------------------------------------------------------------------------
      	//  Mechanics model
      	//  ---------------------------------------------------------------------------
        // set SLset
        SLset = T1.SLsetAll[i + segments*jCount];

      	// Time-Varying Parameters
      	CaiMech   = Cai*1000.0;        //cytosolic calcium concentration from EP model

        // Compute single overlap fractions
    		sovr_cle = max(SL/2.0-(SL-len_thin), len_hbare/2.0);

    		if (len_thick/2.0 < SL/2.0){
    			sovr_ze = len_thick/2.0; // z-line end
    		} else{
    			sovr_ze = SL/2.0;
    		}

    		len_sovr  = sovr_ze-sovr_cle;               // single overlap length
    		SOVFThick = len_sovr*2.0/(len_thick-len_hbare);  // thick filament overlap frac
    		SOVFThin  = len_sovr/len_thin;                 // thin filament overlap frac

    		// Compute combined Ca binding to high- (w/XB) and low- (no XB) sites
    		Tropreg = (1.0-SOVFThin)*TRPNCaL + SOVFThin*TRPNCaH;
    		permtot = sqrt(1.0/(1.0+pow(perm50/Tropreg, nperm)));
    		if (1.0/permtot < 100.0){
    			inprmt  = 1.0/permtot;
    		} else{
    			inprmt  = 100.0;
    		}


    		// Adjustments for Ca activation, temperature, SL, stress and strain
    		konT       = kon*pow(Qkon, ((Temp-310.0)/10.0));

    		if (xXBprer<0.0){
    			signxXBprer = -1.0;
    		}else{
    			signxXBprer = 1.0;
    		}
    		if ((xXBpostr-x_0)<0.0){
    			signxXBpostr = -1.0;
    		}else{
    			signxXBpostr = 1.0;
    		}

    		hfmd    = exp(-signxXBprer*hfmdc*pow((xXBprer/x_0), 2.0));
    		hbmd    = exp(signxXBpostr*hbmdc*pow((xXBpostr-x_0)/x_0, 2.0));

    		// Adjustments for Ca activation, temperature, SL, stress and strain
    		kn_pT   = kn_p*permtot*pow(Qkn_p, ((Temp-310.0)/10.0));
    		kp_nT   = kp_n*inprmt*pow(Qkp_n, ((Temp-310.0)/10.0));
    		fappT   = fapp*xbmodsp*pow(Qfapp, ((Temp-310.0)/10.0));
    		gapSLmd = 1.0 + (1.0-SOVFThick)*gSLmod;
    		gappT   = gapp*gapSLmd*xbmodsp*pow(Qgapp, ((Temp-310.0)/10.0));
    		hfT     = hfXB*hfmd*xbmodsp*pow(Qhf, ((Temp-310.0)/10.0));
    		hbT     = hb*hbmd*xbmodsp*pow(Qhb, ((Temp-310.0)/10.0));
    		koffLT  = koffL*pow(Qkoff, ((Temp-310.0)/10.0))*koffmod;
    		koffHT  = koffH*pow(Qkoff, ((Temp-310.0)/10.0))*koffmod;

    		// steady-state fractions in XBprer and XBpostr using King-Altman rule
    		SSXBprer = (hb*fapp+gxb*fapp)/
    		  (gxb*hfXB+fapp*hfXB+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);
    		SSXBpostr = fapp*hfXB/(gxb*hfXB+fapp*hfXB+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);

    		// normalization for scaling active and passive force (maximal force)
    		Fnordv = kxb*x_0*SSXBpostr;

    		// Calculate forces (active, passive, preload, afterload)
    		force = kxb*SOVFThick*(xXBpostr*XBpostr+xXBprer*XBprer);
    		active = force/Fnordv;

    		if ((SL-SLrest)<0.0){
    			signSL = -1.0;
    		}else{
    			signSL = 1.0;
    		}
    		ppforce_t = signSL*PCon_t*(exp(PExp_t*abs(SL-SLrest))-1.0);
    		if ((SL-SL_c) < 0.0){
    		    ppforce_c = 0.0*PCon_c*(exp(PExp_c*abs(SL-SL_c))-1.0);
    		} else if ((SL-SL_c) == 0.0){
    		    ppforce_c = 1.0/2*PCon_c*(exp(PExp_c*abs(SL-SL_c))-1.0);
    		} else {
    		    ppforce_c = 1.0*PCon_c*(exp(PExp_c*abs(SL-SL_c))-1.0);
    		}
    		ppforce = ppforce_t + ppforce_c;
    		if ((SLset-SLrest)<0.0){
    			int signSLset = -1.0;
    		}else{
    			int signSLset = 1.0;
    		}
    		preload = signSLset*PCon_t*(exp(PExp_t*abs(SLset-SLrest))-1.0);
    		afterload = 0.0;

    		// change in SL and forces
    		lmbda = SL/SLset;
        lmbda = (SL/SLset)*(lmbda>1.0);
        lmbda = lmbda+(lmbda==0.0);


    		Cp = 0.2;
    		b_ff = 50.0;
    		passive = Cp*b_ff*(-0.5 + 0.5*(lmbda*lmbda))*exp(b_ff*((-0.5 +0.5*(lmbda*lmbda))*(-0.5 + 0.5*(lmbda*lmbda))));
    		total_force = -preload + active + passive;
    		dintf = (-ppforce+preload-active+afterload);

        dSL = (-total_force + afterload)/visc;

        lmbda = (SL/SLset)*(lmbda>1.0);
        lmbda = lmbda+(lmbda==0.0);

        //Plateau
        if (t <= Gauss){
            GaussCurve = 1.0/(sqrt(GaussFactor*pi)) * (-1.0/2.0*(2.0*(t-(Gauss))/sigma2)*1.0/sigma2) * exp(-1.0/2.0*pow(((t-(Gauss))/sigma2),2.0));
        }
        else if (t>Gauss && t<400.0-Gauss){
            GaussCurve = 0.0;
        }
        else if (t>=400.0-Gauss){
            double Gauss2 = Gauss + (bcl-2.0*Gauss);
            GaussCurve = 1.0/(sqrt(GaussFactor*pi)) * (-1.0/2.0*(2.0*(t-(Gauss2))/sigma2)*1.0/sigma2) * exp(-1.0/2.0*pow(((t-(Gauss2))/sigma2),2.0));
        }
        if (t>400.0){
            GaussCurve = 0.0;
        }

        // internal shortening
        if (i != 0 && i != (segments-1)){
          dSL = -(dSL1[i-1+segments*jCount] + dSL1[i+1+segments*jCount] - 2.0*dSL1[i+segments*jCount]) + GaussCurve;
        } else if (i == 0){
          if (jCount != 0){
            dSL = -(dSL1[i+1+segments*jCount] + dSL1[(segments-1)+segments*(jCount-1)] - 2.0*dSL1[i+segments*jCount]) + GaussCurve;
          } else if (jCount == 0){
            dSL = -(dSL1[i+1+segments*jCount] - dSL1[i+segments*jCount]) + GaussCurve;
          }
        } else if(i == (segments-1)){
          if (jCount != (cells-1)){
            dSL = -(dSL1[i-1+segments*jCount] + dSL1[0+segments*(jCount+1)] - 2.0*dSL1[i+segments*jCount]) + GaussCurve;
          }else if (jCount == (cells-1)){
            dSL = -(dSL1[i-1+segments*jCount] - dSL1[i+segments*jCount]) + GaussCurve;
          }
        }


    		// change to Rice model
    		sigman = xOpti[6]*exp(xOpti[4]*(lmbda-1.0));
    		sigmap = xOpti[7]*exp(xOpti[5]*(lmbda-1.0));
    		if ((x_0-xXBpostr)<0.0){
    			heavisidePostR = 0.0;
    		}else{
    			heavisidePostR = 1.0;
    		}

    		gxbmd   = heavisidePostR*exp(sigmap*pow(((x_0-xXBpostr)/x_0), 2.0))+
    		  (1.0-heavisidePostR)*exp(sigman*pow(((xXBpostr-x_0)/x_0), 2.0));

    		gxbT    = gxb*gxbmd*xbmodsp*pow(Qgxb, ((Temp-310.0)/10.0));

    		//
    		// Regulation and corssbridge cycling state derivatives
    		dTRPNCaL  = konT*(exp(xOpti[2]*(lmbda-1.0)))*CaiMech*(1.0-TRPNCaL) - koffLT*pow(exp(xOpti[3]*(lmbda-1.0)), (-1.0))*TRPNCaL;
    		dTRPNCaH  = konT*(exp(xOpti[2]*(lmbda-1.0)))*CaiMech*(1.0-TRPNCaH) - koffHT*pow(exp(xOpti[3]*(lmbda-1.0)), (-1.0))*TRPNCaH;

    		dN_NoXB   = -kn_pT*N_NoXB+ kp_nT*P_NoXB;
    		dP_NoXB   = -kp_nT*P_NoXB + kn_pT*N_NoXB;

    		dN        = -kn_pT*N+ kp_nT*P;
    		// dP      = -kp_nT*P + kn_pT*N - fappT*P + gappT*XBprer + gxbT*XBpostr;
    		dXBprer   = fappT*P - gappT*XBprer - hfT*XBprer + hbT*XBpostr;
    		dXBpostr  = hfT*XBprer - hbT*XBpostr - gxbT*XBpostr;

    		dP        = -(dN+dXBprer+dXBpostr);

    		// Mean strain of strongly-bound states due to SL motion and XB cycling
    		dutyprer  = (hbT*fappT+gxbT*fappT)/
    			(fappT*hfT+gxbT*hfT+gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT); // duty fractions using the
    		dutypostr = fappT*hfT/
    		  (fappT*hfT+gxbT*hfT+gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT);// King-Alman Rule
    		dxXBprer = 0.5*dSL+xPsi/dutyprer*(-xXBprer*fappT+(xXBpostr-x_0-xXBprer)*hbT);
    		dxXBpostr = 0.5*dSL+ xPsi/dutypostr*(x_0+xXBprer-xXBpostr)*hfT;

    		// Ca buffering by low-affinity troponin C (LTRPNCa)
    		FrSBXB    = (XBpostr+XBprer)/(SSXBpostr + SSXBprer);
    		dFrSBXB   = (dXBpostr+dXBprer)/(SSXBpostr + SSXBprer);

    		if ((len_thick-SL) < 0.0){
    		   dsovr_ze  = -dSL/2.0*0.0;
    		} else if ((len_thick-SL) == 0.0){
    		   dsovr_ze  = -dSL/2.0*1.0/2.0;
    		} else {
    		   dsovr_ze  = -dSL/2.0*1.0;
    		}

    		if (((2.0*len_thin-SL)-len_hbare) < 0.0){
    		    dsovr_cle = -dSL/2.0*0.0;
    		} else if (((2.0*len_thin-SL)-len_hbare) == 0.0){
    		    dsovr_cle = -dSL/2.0*1.0/2.0;
    		} else{
    		    dsovr_cle = -dSL/2.0*1.0;
    		}

    		dlen_sovr = dsovr_ze-dsovr_cle;
    		dSOVFThin = dlen_sovr/len_thin;
    		dSOVFThick= 2.0*dlen_sovr/(len_thick-len_hbare);

    		TropTot = Trop_conc*((1.0-SOVFThin)*TRPNCaL +
    		  SOVFThin*(FrSBXB*TRPNCaH+(1.0-FrSBXB)*TRPNCaL));
    		dTropTot= (-dSOVFThin*TRPNCaL+(1.0-SOVFThin)*dTRPNCaL +
    		  dSOVFThin*(FrSBXB*TRPNCaH+(1.0-FrSBXB)*TRPNCaL) +
    		  SOVFThin*(dFrSBXB*TRPNCaH+FrSBXB*dTRPNCaH-dFrSBXB*TRPNCaL+
    		  (1.0-FrSBXB)*dTRPNCaL));


    		dforce = kxb*dSOVFThick*(xXBpostr*XBpostr+xXBprer*XBprer) +
    		  kxb*SOVFThick*(dxXBpostr*XBpostr+xXBpostr*dXBpostr +
    		  dxXBprer*XBprer+xXBprer*dXBprer);

    		dactive = 0.5*dforce/Fnordv;
    		dppforce_t = signSL*PCon_t*PExp_t*dSL*exp(PExp_t*abs(SL-SLrest));

    		if ((SL-SL_c) < 0.0){
    		    dppforce_c = 0.0*PCon_c*PExp_c*dSL*exp(PExp_c*abs(SL-SL_c));
    		} else if ((SL-SL_c) == 0.0){
    		    dppforce_c = 1.0/2.0*PCon_c*PExp_c*dSL*exp(PExp_c*abs(SL-SL_c));
    		} else{
    		    dppforce_c = 1.0*PCon_c*PExp_c*dSL*exp(PExp_c*abs(SL-SL_c));
    		}

    		dppforce = dppforce_t + dppforce_c;
    		dsfib = dppforce+dactive;


        // ---------------------------------------------------------------------------
        //  EP model
        //  ---------------------------------------------------------------------------

        // interaction between agjacent SR, cytosol and SL due to fluxes (based on Voigt et al.)
        if (i != 0 && i != (segments-1)){
          J_diff_cyt = (x[9 + StateNum*(i-1) + StateNum*(segments*jCount)]
                       + x[9 + StateNum*(i+1) + StateNum*(segments*jCount)]
                       - 2.0*x[9 + StateNum*i + StateNum*(segments*jCount)]) / tau_seg_cyt;
          J_diff_SR  = (x[7 + StateNum*(i-1) + StateNum*(segments*jCount)]
                       + x[7 + StateNum*(i+1) + StateNum*(segments*jCount)]
                       - 2.0*x[7 + StateNum*i + StateNum*(segments*jCount)]) / tau_seg_SR;
          J_diff_SL  = (x[2 + StateNum*(i-1) + StateNum*(segments*jCount)]
                       + x[2 + StateNum*(i+1) + StateNum*(segments*jCount)]
                       - 2.0*x[2 + StateNum*i + StateNum*(segments*jCount)]) / tau_seg_SL;
        } else if (i == 0){
          J_diff_cyt = (x[9 + StateNum*(i+1) + StateNum*(segments*jCount)]
                       - x[9 + StateNum*i + StateNum*(segments*jCount)]) / tau_seg_cyt;
          J_diff_SR  = (x[7 + StateNum*(i+1) + StateNum*(segments*jCount)]
                       - x[7 + StateNum*i + StateNum*(segments*jCount)]) / tau_seg_SR;
          J_diff_SL  = (x[2 + StateNum*(i+1) + StateNum*(segments*jCount)]
                       - x[2 + StateNum*i + StateNum*(segments*jCount)]) / tau_seg_SL;
        } else if(i == (segments-1)){
          J_diff_cyt = (x[9 + StateNum*(i-1) + StateNum*(segments*jCount)]
                       - x[9 + StateNum*i + StateNum*(segments*jCount)]) / tau_seg_cyt;
          J_diff_SR  = (x[7 + StateNum*(i-1) + StateNum*(segments*jCount)]
                       - x[7 + StateNum*i + StateNum*(segments*jCount)]) / tau_seg_SR;
          J_diff_SL  = (x[2 + StateNum*(i-1) + StateNum*(segments*jCount)]
                       - x[2 + StateNum*i + StateNum*(segments*jCount)]) / tau_seg_SL;
        }

        // Nernst potential
        E_Na_SL     = R2*T/F*log(Nao/Na_SL);
        E_Na_jct    = R2*T/F*log(Nao/Na_jct);
        E_Ca_jct    = R2*T/(2.0*F)*log(Cao/Ca_jct);
        E_Ks        = R2*T/F*log((K_o+pKNa*Nao)/(Ki+pKNa*Nai));
        E_Ca_SL     = R2*T/(2.0*F)*log(Cao/Ca_SL);

        // Conductances
        pCa_jct     = -log10(Ca_jct/1.0)+3.0;
        pCa_SL      = -log10(Ca_SL/1.0)+3.0;
        G_Ks_jct    = 0.07*(0.057+0.19/(1.0+exp((-7.2+pCa_jct)/0.6)));
        G_Ks_SL     = 0.07*(0.057+0.19/(1.0+exp((-7.2+pCa_SL)/0.6)));


        // Membrane CuRrents
        // I_Na: Fast Na CuRrent
        alpha_m = 0.32*(Vm+47.13)/1.0/(1.0-exp(-0.1*(Vm+47.13)));
        beta_m = 0.08*exp(-Vm/11.0);

        if (Vm < -40.0){
          alpha_h = 0.135*exp((80.0+Vm)/-6.8);
          beta_h = 3.56*exp(0.079*Vm)+3.1e5*exp(0.35*Vm);
          alpha_j = (-1.2714e5*exp(0.2444*Vm)-3.474e-5*exp(-0.04391*Vm))*(Vm+37.78)/1.0/(1.0+exp(0.311*(Vm+79.23)));
          beta_j = 0.1212*exp(-0.01052*Vm)/(1.0+exp(-0.1378*(Vm+40.14)));
        }else{
          alpha_h = 0.0;
          beta_h = 1.0/(0.13*(1.0+exp((Vm+10.66)/-11.1)));
          alpha_j = 0.0;
          beta_j = 0.3*exp(-2.535e-7*Vm)/(1.0+exp(-0.1*(Vm+32.0)));
        }

        dh = alpha_h*(1.0-h)-beta_h*h;
        dj = alpha_j*(1.0-j)-beta_j*j;
        dm = alpha_m*(1.0-m)-beta_m*m;


        openProb = pow(m, 3.0)*h*j;
        i_Na_jct = f_inhib*Fx_Na_jct*G_INa*openProb*(Vm-E_Na_jct);
        i_Na_SL = f_inhib*Fx_Na_SL*G_INa*openProb*(Vm-E_Na_SL);
        i_Na = i_Na_jct+i_Na_SL;

        // I_nabk: Na Background CuRrent
        i_Nab_jct = Fx_NaBk_jct*G_NaBk*(Vm-E_Na_jct);
        i_Nab_SL = Fx_NaBk_SL*G_NaBk*(Vm-E_Na_SL);
        i_Nab = i_Nab_jct+i_Nab_SL;

        // I_nak: Na/K Pump CuRrent
        f_NaK = 1.0/(1.0+0.1245*exp(-0.1*Vm*F/(R2*T))+0.0365*sigma*exp(-Vm*F/(R2*T)));
        i_NaK_jct = Fx_NaK_jct*I_NaK_max*f_NaK/(1.0+pow(Km_Nai/Na_jct, H_NaK))*K_o/(K_o+Km_Ko);
        i_NaK_SL = Fx_NaK_SL*I_NaK_max*f_NaK/(1.0+pow(Km_Nai/Na_SL, H_NaK))*K_o/(K_o+Km_Ko);
        i_NaK = ScaleFactorGNaK*(i_NaK_jct+i_NaK_SL);

        // i_Kr: Rapidly Activating K CuRrent
        Rr = 1.0/(1.0+exp((33.0+Vm)/22.4));
        i_Kr = ScaleFactorGKr*(G_IKr*Xr*Rr*(Vm-E_K));

        // I_ks: Slowly Activating K CuRrent
        Xr_infinity = 1.0/(1.0+exp(-(50.0+Vm)/7.5));
        tau_Xr = 1.0/(0.00138*(Vm+7.0)/(1.0-exp(-0.123*(Vm+7.0)))+0.00061*(Vm+10.0)/(exp(0.145*(Vm+10.0))-1.0));
        dXr = (Xr_infinity-Xr)/tau_Xr;

        Xs_infinity = 1.0/(1.0+exp(-(Vm-1.5)/16.7));
        tau_Xs = 1.0/(7.19e-5*(Vm+30.0)/(1.0-exp(-0.148*(Vm+30.0)))+1.31e-4*(Vm+30.0)/(-1.0+exp(0.0687*(Vm+30.0))));
        dXs = (Xs_infinity-Xs)/tau_Xs;

        i_Ks_jct = Fx_Ks_jct*G_Ks_jct*pow(Xs, 2.0)*(Vm-E_Ks);
        i_Ks_SL = Fx_Ks_SL*G_Ks_SL*pow(Xs, 2.0)*(Vm-E_Ks);
        i_Ks = ScaleFactorGKs*(i_Ks_jct+i_Ks_SL);

        // I_kp: Plateau K cuRrent
        i_Kp = g_Kp*(Vm-E_K)/(1.0+exp(7.488-Vm/5.98));

        // I_to: Transient Outward K CuRrent (slow and fast components)
        i_tof = G_tof*X_tof*Y_tof*(Vm-E_K);
        i_tos = G_tos*X_tos*(Y_tos+0.5*R_tof)*(Vm-E_K);
        i_to = ScaleFactorGto*(i_tof+i_tos);

        X_tof_infinity = 1.0/(1.0+exp(-(Vm+3.0)/15.0));
        tau_X_tof = 3.5*exp(-pow(Vm/30.0, 2.0))+1.5;
        dX_tof = (X_tof_infinity-X_tof)/tau_X_tof;
        Y_tof_infinity = 1.0/(1.0+exp((Vm+33.5)/10.0));
        tau_Y_tof = 20.0/(1.0+exp((Vm+33.5)/10.0))+20.0;
        dY_tof = (Y_tof_infinity-Y_tof)/tau_Y_tof;
        R_tof_infinity = 1.0/(1.0+exp((Vm+33.5)/10.0));
        tau_R_tof = 2.8e3/(1.0+exp((Vm+60.0)/10.0))+220.0;
        dR_tof = (R_tof_infinity-R_tof)/tau_R_tof;
        X_tos_infinity = 1.0/(1.0+exp(-(Vm+3.0)/15.0));
        tau_X_tos = 9.0/(1.0+exp((Vm+3.0)/15.0))+0.5;
        dX_tos = (X_tos_infinity-X_tos)/tau_X_tos;
        Y_tos_infinity = 1.0/(1.0+exp((Vm+33.5)/10.0));
        tau_Y_tos = 3000.0/(1.0+exp((Vm+60.0)/10.0))+30.0;
        dY_tos = (Y_tos_infinity-Y_tos)/tau_Y_tos;

        // I_ki: Time-Independent K CuRrent
        alpha_K1 = 1.02/(1.0+exp(0.2385*(Vm-(E_K+59.215))));
        beta_K1 = (0.49124*exp(0.08032*(Vm-E_K+5.476))+1.0*exp(0.06175*(Vm-(E_K+594.31))))/(1.0+exp(-0.5143*(Vm-E_K+4.753)));
        K1_infinity = alpha_K1/(alpha_K1+beta_K1);
        i_K1 = ScaleFactorGK1*(G_K1*K1_infinity*(Vm-E_K));

        // I_ClCa: Ca-activated Cl CuRrent).^(i_Clbk: background Cl CuRrent
        i_Cl_Ca = G_Cl*(Vm-E_Cl)*(Fx_Cl_jct/(1.0+Kd_ClCa/Ca_jct)+Fx_Cl_SL/(1.0+Kd_ClCa/Ca_SL));
        i_Clb = G_ClBk*(Vm-E_Cl);

        // I_Ca: L-type Calcium CuRrent
        d_infinity = 1.0/(1.0+exp(-(Vm+14.5)/6.0));
        tau_d = 1.0*d_infinity*(1.0-exp(-(Vm+14.5)/6.0))/(0.035*(Vm+14.5));
        dd = (d_infinity-d)/tau_d;
        dfCaB_SL = 1.7*Ca_SL*(1.0-fCaB_SL)-11.9e-3*fCaB_SL;
        dfCaB_jct = 1.7*Ca_jct*(1.0-fCaB_jct)-11.9e-3*fCaB_jct;
        f_infinity = 1.0/(1.0+exp((Vm+35.06)/3.6))+0.6/(1.0+exp((50.0-Vm)/20.0));
        tau_f = 1.0/(0.0197*exp(-pow(0.0337*(Vm+14.5), 2.0))+0.02);
        dfEP = (f_infinity-fEP)/tau_f;

        temp = 0.45*d*fEP*Q_CaL*Vm*pow(F, 2.0)/(R2*T);
        fCa_jct = 1.0-fCaB_jct;
        fCa_SL = 1.0-fCaB_SL;

        i_CaL_Ca_SL = temp*fCa_SL*Fx_ICaL_SL*PCa*4.0*(gamma_Cai*Ca_SL*exp(2.0*Vm*F/(R2*T))-gamma_Cao*Cao)/(exp(2.0*Vm*F/(R2*T))-1.0);
        i_CaL_Ca_jct = temp*fCa_jct*Fx_ICaL_jct*PCa*4.0*(gamma_Cai*Ca_jct*exp(2.0*Vm*F/(R2*T))-gamma_Cao*Cao)/(exp(2.0*Vm*F/(R2*T))-1.0);

        i_CaL_Na_jct = temp*fCa_jct*Fx_ICaL_jct*PNa*(gamma_Nai*Na_jct*exp(Vm*F/(R2*T))-gamma_Nao*Nao)/(exp(Vm*F/(R2*T))-1.0);
        i_CaL_Na_SL = temp*fCa_SL*Fx_ICaL_SL*PNa*(gamma_Nai*Na_SL*exp(Vm*F/(R2*T))-gamma_Nao*Nao)/(exp(Vm*F/(R2*T))-1.0);

        i_CaL_K = temp*(fCa_SL*Fx_ICaL_SL+fCa_jct*Fx_ICaL_jct)*PK*(gamma_Ki*Ki*exp(Vm*F/(R2*T))-gamma_Ko*K_o)/(exp(Vm*F/(R2*T))-1.0);

        i_CaL = f_inhib*ScaleFactorGCaL*(i_CaL_Ca_SL+i_CaL_Ca_jct+i_CaL_Na_SL+i_CaL_Na_jct+i_CaL_K);


        // I_ncx: Na/Ca Exchanger flux
        Ka_jct = 1.0/(1.0+pow(Kd_act/Ca_jct, 3.0));
        temp_jct = (exp(eta*Vm*F/(R2*T))*pow(Na_jct, HNa)*Cao-exp((eta-1.0)*Vm*F/(R2*T))*pow(Nao, HNa)*Ca_jct)/(1.0+ksat*exp((eta-1.0)*Vm*F/(R2*T)));
        i_NaCa_jct = Fx_NCX_jct*V_max_1*Ka_jct*Q_NCX*temp_jct/(K_mCai*pow(Nao, HNa)*(1.0+pow(Na_jct/K_mNai, HNa))
            +pow(K_mNao, HNa)*Ca_jct*(1.0+Ca_jct/K_mCai)+K_mCao*pow(Na_jct, HNa)+pow(Na_jct, HNa)*Cao+pow(Nao, HNa)*Ca_jct);

        Ka_SL = 1.0/(1.0+pow(Kd_act/Ca_SL, 3.0));
        temp_SL = (exp(eta*Vm*F/(R2*T))*pow(Na_SL, HNa)*Cao-exp((eta-1.0)*Vm*F/(R2*T))*pow(Nao, HNa)*Ca_SL)/(1.0+ksat*exp((eta-1.0)*Vm*F/(R2*T)));
        i_NaCa_SL = Fx_NCX_SL*V_max_1*Ka_SL*Q_NCX*temp_SL/(K_mCai*pow(Nao, HNa)*(1.0+pow(Na_SL/K_mNai, HNa))+
            pow(K_mNao, HNa)*Ca_SL*(1.0+Ca_SL/K_mCai)+K_mCao*pow(Na_SL, HNa)+pow(Na_SL, HNa)*Cao+pow(Nao, HNa)*Ca_SL);

        i_NaCa = i_NaCa_jct+i_NaCa_SL;

        // I_pca: Sarcolemmal Ca Pump CuRrent
        j_pump_SR = Q_SRCaP*V_max_2*(pow(Cai/Kmf, H2)-pow(Ca_SR/Kmr, H2))/(1.0+pow(Cai/Kmf, H2)+pow(Ca_SR/Kmr, H2));

        i_Cap_jct = Q_SLCaP*V_maxAF*Fx_SLCaP_jct/(1.0+pow(Km/Ca_jct, H1));
        i_Cap_SL = Q_SLCaP*V_maxAF*Fx_SLCaP_SL/(1.0+pow(Km/Ca_SL, H1));

        i_Cap = i_Cap_jct+i_Cap_SL;

        // I_cabk: Ca Background CuRrent
        i_Cab_jct = G_CaBk*Fx_CaBk_jct*(Vm-E_Ca_jct);
        i_Cab_SL = G_CaBk*Fx_CaBk_SL*(Vm-E_Ca_SL);

        i_Cab = i_Cab_SL+i_Cab_jct;

        // SR fluxes: Calcium Release).^(SR Ca pump).^(SR Ca leak
        kCaSR = Max_SR-(Max_SR-Min_SR)/(1.0+pow(EC50_SR/Ca_SR, HSR));
        koSRCa = koCa/kCaSR;
        kiSRCa = kiCa*kCaSR;
        RI = 1.0-R1-O-I;
        dR1 = kim*RI-kiSRCa*Ca_jct*R1-(koSRCa*pow(Ca_jct, 2.0)*R1-kom*O);
        dO = koSRCa*pow(Ca_jct, 2.0)*R1-kom*O-(kiSRCa*Ca_jct*O-kim*I);
        dI = kiSRCa*Ca_jct*O-kim*I-(kom*I-koSRCa*pow(Ca_jct, 2.0)*RI);

        i_Ca_SL_tot = i_CaL_Ca_SL-2.0*i_NaCa_SL+i_Cab_SL+i_Cap_SL;

        j_leak_SR = KSRleak*(Ca_SR-Ca_jct);
        j_rel_SR = ks*O*(Ca_SR-Ca_jct);

        // Sodium and Calcium Buffering
        dNa_jct_buf = kon_EP*Na_jct*(Bmax_jct-Na_jct_buf)-koff_EP*Na_jct_buf;
        dNa_SL_buf = kon_EP*Na_SL*(Bmax_SL-Na_SL_buf)-koff_EP*Na_SL_buf;

        // Cytosolic Ca Buffers
        //dCa_TroponinC = kon_TroponinC*Cai*(Bmax_TroponinC-Ca_TroponinC)-koff_TroponinC*Ca_TroponinC;
        dCa_TroponinC = Bmax_TroponinC*dTropTot;
        dCa_TroponinC_Ca_Mg = kon_TroponinC_Ca_Mg_Ca*Cai*(Bmax_TroponinC_Ca_Mg_Ca-(Ca_TroponinC_Ca_Mg+Mg_TroponinC_Ca_Mg))-koff_TroponinC_Ca_Mg_Ca*Ca_TroponinC_Ca_Mg;
        dMg_TroponinC_Ca_Mg = kon_TroponinC_Ca_Mg_Mg*Mgi*(Bmax_TroponinC_Ca_Mg_Mg-(Ca_TroponinC_Ca_Mg+Mg_TroponinC_Ca_Mg))-koff_TroponinC_Ca_Mg_Mg*Mg_TroponinC_Ca_Mg;
        dCa_Calmodulin = kon_Calmodulin*Cai*(Bmax_Calmodulin-Ca_Calmodulin)-koff_Calmodulin*Ca_Calmodulin;
        dCa_Myosin = kon_Myosin_Ca*Cai*(Bmax_Myosin_Ca-(Ca_Myosin+Mg_Myosin))-koff_Myosin_Ca*Ca_Myosin;
        dMg_Myosin = kon_Myosin_Mg*Mgi*(Bmax_Myosin_Mg-(Ca_Myosin+Mg_Myosin))-koff_Myosin_Mg*Mg_Myosin;
        dCa_SRB = kon_SRB*Cai*(Bmax_SRB-Ca_SRB)-koff_SRB*Ca_SRB;
        dCa_cytosol_tot_bound = dCa_TroponinC+dCa_TroponinC_Ca_Mg+dMg_TroponinC_Ca_Mg+dCa_Calmodulin+dCa_Myosin+dMg_Myosin+dCa_SRB;

        // Junctional and SL Ca Buffers
        dCa_SLB_SL = kon_SL*Ca_SL*(Bmax_SLB_SL*Vol_myo/Vol_SL-Ca_SLB_SL)-koff_SLB*Ca_SLB_SL;
        dCa_SLB_jct = kon_SL*Ca_jct*(Bmax_SLB_jct*0.1*Vol_myo/Vol_jct-Ca_SLB_jct)-koff_SLB*Ca_SLB_jct;
        dCa_SLHigh_SL = kon_SL*Ca_SL*(Bmax_SLHigh_SL*Vol_myo/Vol_SL-Ca_SLHigh_SL)-koff_SLHigh*Ca_SLHigh_SL;
        dCa_SLHigh_jct = kon_SL*Ca_jct*(Bmax_SLHigh_jct*0.1*Vol_myo/Vol_jct-Ca_SLHigh_jct)-koff_SLHigh*Ca_SLHigh_jct;
        dCa_jct_tot_bound = dCa_SLB_jct+dCa_SLHigh_jct;
        dCa_SL_tot_bound = dCa_SLB_SL+dCa_SLHigh_SL;

        // Ion concentrations
        // SR Ca Concentrations
        dCalsequestrin = kon_Calsequestrin*Ca_SR*(Bmax_Calsequestrin*Vol_myo/Vol_SR-Ca_Calsequestrin)-koff_Calsequestrin*Ca_Calsequestrin;

        // Sodium Concentrations
        J_Na_jct_SL = (Na_jct-Na_SL)*Diff_Na_jct_SL;
        dNa_jct = -Cm*(i_Na_jct+3.0*i_NaCa_jct+i_Nab_jct+3.0*i_NaK_jct+i_CaL_Na_jct)/(Vol_jct*F)-J_Na_jct_SL/Vol_jct-dNa_jct_buf;
        J_Na_SL_myo = (Na_SL-Nai)*Diff_Na_SL_myo;
        dNa_SL = -Cm*(i_Na_SL+3.0*i_NaCa_SL+i_Nab_SL+3.0*i_NaK_SL+i_CaL_Na_SL)/(Vol_SL*F)+(J_Na_jct_SL-J_Na_SL_myo)/Vol_SL-dNa_SL_buf;
        dNai = J_Na_SL_myo/Vol_myo;

        // Calcium concentrations
        i_Ca_jct_tot = i_CaL_Ca_jct-2.0*i_NaCa_jct+i_Cab_jct+i_Cap_jct;

        dCa_SR = j_pump_SR-(j_leak_SR*Vol_myo/Vol_SR+j_rel_SR)-dCalsequestrin + J_diff_SR;
        J_Ca_jct_SL = (Ca_jct-Ca_SL)*Diff_Ca_jct_SL;

        dCa_jct = -i_Ca_jct_tot*Cm/(Vol_jct*2.0*F)-J_Ca_jct_SL/Vol_jct+j_rel_SR*Vol_SR/Vol_jct+j_leak_SR*Vol_myo/Vol_jct-1.0*dCa_jct_tot_bound;
        J_Ca_SL_myo = (Ca_SL-Cai)*Diff_Ca_SL_myo;

        dCa_SL = -i_Ca_SL_tot*Cm/(Vol_SL*2.0*F)+(J_Ca_jct_SL-J_Ca_SL_myo)/Vol_SL-1.0*dCa_SL_tot_bound + J_diff_SL;
        dCai = -j_pump_SR*Vol_SR/Vol_myo+J_Ca_SL_myo/Vol_myo-1.0*dCa_cytosol_tot_bound + J_diff_cyt;


        // stimulus
        i_Stim = 0.0;
        if (t>10.1 & t<13.1){
          i_Stim = -stim_amplitude;
        }

        // K-ATP cuRrent
        a_katp = K_o/K_o_normal;
        i_katp = f_katp*g_katp*pow(a_katp, 0.24)*(R1-E_K);


        // fibroblasts
        E_fK = Rf*Tf/Ff * log(K_fo/K_fi);
        E_fNa = Rf*Tf/Ff * log(Na_fo/Na_fi);


        // Time- and voltage-dependent fibroblasts K+ current
        tau_r = 0.0203 + 0.138 * exp(pow((-(Vf + 20.0)/25.9),2.0));
        tau_s = 1.574 + 5.268 * exp(pow((-(Vf + 23.0)/22.7),2.0));
        r_inf  = 1.0/(1.0 + exp(-((Vf + 20.0)/11.0)));
        s_inf  = 1.0/(1.0 + exp((Vf + 23.0)/7.0));

        dr_KV= (r_inf - r_KV)/tau_r;
        ds_KV = (s_inf - s_KV)/tau_s;

        I_KV = g_KV*r_KV*s_KV*(Vm-E_fK);


        // Inward-rectfying K+ current in the fibroblasts
        alpha_fK1 = 0.1/(1.0 + exp(0.06*(Vf - E_fK - 200.0)));
        beta_fK1  = (3.0*exp(0.0002*(Vf - E_fK + 100.0)) + exp(0.1*(Vf - E_fK - 10.0)))/(1.0 + exp(-0.5*(Vf - E_fK)));

        I_fK1 = (g_fK1 * alpha_fK1 * (Vf - E_fK))/(alpha_fK1 + beta_fK1);


        // Na+ - K+ pump current in the fibroblast
        I_fNaK = I_fNaK_inf * (K_fo/(K_fo + K_fmK)) * (pow((Na_fi/(Na_fi + K_fmNa)),(3.0/2.0))) * ((Vf - V_fNaK_rev) / (Vf - B_fNaK));


        // Background Na+ current in the fibroblast
        I_fbNa = G_fbNa * (Vf - E_fNa);


        // intracellular ionic concentrations
        dK_fi = -(I_fK1 + I_KV - 2.0*I_fNaK)/(Vol_f*Ff);
        dNa_fi = -(I_fbNa + 3.0*I_fNaK)/(Vol_f*Ff);


        // total membrane current
        I_f = Cmf * (I_KV + I_fK1 + I_fNaK + I_fbNa);


        // Gap junction between myocyte and fibroblast
        E_NaGap = Rf*Tf/Ff * log(K_fi/Ki);
        E_KGap  = Rf*Tf/Ff * log(Na_fi/Nai);

        GapFib = T1.GapAll[i + segments*jCount];
        I_gapNa = GapFib * ((Vm - Vf)-E_NaGap);
        I_gapK = GapFib * ((Vm - Vf)-E_KGap);

        I_gap = Cmf * (I_gapNa + I_gapK)*GapFib;


        // membrane potential fibroblast
        dVf = -(I_f -I_gap)*GapFib;


        // membrane potential
        if (i != 0 && i != (segments-1)){
          J_diff_Vm = (x[0 + StateNum*(i-1) + StateNum*(segments*jCount)]
                      + x[0 + StateNum*(i+1) + StateNum*(segments*jCount)]
                      - 2.0*x[0 + StateNum*i + StateNum*(segments*jCount)]) / tau_Vm;
        } else if (i == 0){
          J_diff_Vm = (x[0 + StateNum*(i+1) + StateNum*(segments*jCount)]
                      - x[0 + StateNum*i + StateNum*(segments*jCount)]) / tau_Vm;
        } else if(i == (segments-1)){
          J_diff_Vm = (x[0 + StateNum*(i-1) + StateNum*(segments*jCount)]
                      - x[0 + StateNum*i + StateNum*(segments*jCount)]) / tau_Vm;
        }

        dVm = -(i_Na+i_Nab+i_NaK+i_Kr+i_Ks+i_to+i_K1+i_NaCa+
          i_Cl_Ca+i_Clb+i_CaL+i_Cab+i_Cap+i_Kp+i_katp+I_gap)+J_diff_Vm;

        if (i == 0 && jCount == 0){
          dVm = dVm - i_Stim;
        }

        // SolutionOutput
        dxdt[ 0+ StateNum*i + StateNum*(segments*jCount)]     = dVm;
        dxdt[ 1+ StateNum*i + StateNum*(segments*jCount)]     = dCalsequestrin;
        dxdt[ 2+ StateNum*i + StateNum*(segments*jCount)]     = dCa_SL;
        dxdt[ 3+ StateNum*i + StateNum*(segments*jCount)]     = dCa_SLB_SL;
        dxdt[ 4+ StateNum*i + StateNum*(segments*jCount)]     = dCa_SLB_jct;
        dxdt[ 5+ StateNum*i + StateNum*(segments*jCount)]     = dCa_SLHigh_SL;
        dxdt[ 6+ StateNum*i + StateNum*(segments*jCount)]     = dCa_SLHigh_jct;
        dxdt[ 7+ StateNum*i + StateNum*(segments*jCount)]     = dCa_SR;
        dxdt[ 8+ StateNum*i + StateNum*(segments*jCount)]     = dCa_jct;
        dxdt[ 9+ StateNum*i + StateNum*(segments*jCount)]     = dCai;
        dxdt[10+ StateNum*i + StateNum*(segments*jCount)]     = dd;
        dxdt[11+ StateNum*i + StateNum*(segments*jCount)]     = dfCaB_SL;
        dxdt[12+ StateNum*i + StateNum*(segments*jCount)]     = dfCaB_jct;
        dxdt[13+ StateNum*i + StateNum*(segments*jCount)]     = dfEP;
        dxdt[14+ StateNum*i + StateNum*(segments*jCount)]     = dXr;
        dxdt[15+ StateNum*i + StateNum*(segments*jCount)]     = dXs;
        dxdt[16+ StateNum*i + StateNum*(segments*jCount)]     = dh;
        dxdt[17+ StateNum*i + StateNum*(segments*jCount)]     = dj;
        dxdt[18+ StateNum*i + StateNum*(segments*jCount)]     = dm;
        dxdt[19+ StateNum*i + StateNum*(segments*jCount)]     = dX_tof;
        dxdt[20+ StateNum*i + StateNum*(segments*jCount)]     = dY_tof;
        dxdt[21+ StateNum*i + StateNum*(segments*jCount)]     = dR_tof;
        dxdt[22+ StateNum*i + StateNum*(segments*jCount)]     = dX_tos;
        dxdt[23+ StateNum*i + StateNum*(segments*jCount)]     = dY_tos;
        dxdt[24+ StateNum*i + StateNum*(segments*jCount)]     = dI;
        dxdt[25+ StateNum*i + StateNum*(segments*jCount)]     = dO;
        dxdt[26+ StateNum*i + StateNum*(segments*jCount)]     = dR1;
        dxdt[27+ StateNum*i + StateNum*(segments*jCount)]     = dNa_SL;
        dxdt[28+ StateNum*i + StateNum*(segments*jCount)]     = dNa_SL_buf;
        dxdt[29+ StateNum*i + StateNum*(segments*jCount)]     = dNa_jct;
        dxdt[30+ StateNum*i + StateNum*(segments*jCount)]     = dNa_jct_buf;
        dxdt[31+ StateNum*i + StateNum*(segments*jCount)]     = dNai;
        dxdt[32+ StateNum*i + StateNum*(segments*jCount)]     = dCa_Calmodulin;
        dxdt[33+ StateNum*i + StateNum*(segments*jCount)]     = dCa_Myosin;
        dxdt[34+ StateNum*i + StateNum*(segments*jCount)]     = dCa_SRB;
        dxdt[35+ StateNum*i + StateNum*(segments*jCount)]     = dCa_TroponinC;
        dxdt[36+ StateNum*i + StateNum*(segments*jCount)]     = dCa_TroponinC_Ca_Mg;
        dxdt[37+ StateNum*i + StateNum*(segments*jCount)]     = dMg_Myosin;
        dxdt[38+ StateNum*i + StateNum*(segments*jCount)]     = dMg_TroponinC_Ca_Mg;

        // output mechanics model
    		dxdt[39+ StateNum*i + StateNum*(segments*jCount)]     = dN_NoXB;
    		dxdt[40+ StateNum*i + StateNum*(segments*jCount)]     = dP_NoXB;
    		dxdt[41+ StateNum*i + StateNum*(segments*jCount)]     = dN;
    		dxdt[42+ StateNum*i + StateNum*(segments*jCount)]     = dP;
    		dxdt[43+ StateNum*i + StateNum*(segments*jCount)]     = dXBprer;
    		dxdt[44+ StateNum*i + StateNum*(segments*jCount)]     = dXBpostr;
    		dxdt[45+ StateNum*i + StateNum*(segments*jCount)]     = dSL;
    		dxdt[46+ StateNum*i + StateNum*(segments*jCount)]     = dxXBpostr;
    		dxdt[47+ StateNum*i + StateNum*(segments*jCount)]     = dxXBprer;
    		dxdt[48+ StateNum*i + StateNum*(segments*jCount)]     = dTRPNCaL;
    		dxdt[49+ StateNum*i + StateNum*(segments*jCount)]     = dTRPNCaH;
    		dxdt[50+ StateNum*i + StateNum*(segments*jCount)]     = dforce;
    		dxdt[51+ StateNum*i + StateNum*(segments*jCount)]     = dintf;

        // output fibroblast
    		dxdt[52+ StateNum*i + StateNum*(segments*jCount)]     = dVf;
    		dxdt[53+ StateNum*i + StateNum*(segments*jCount)]     = dr_KV;
    		dxdt[54+ StateNum*i + StateNum*(segments*jCount)]     = ds_KV;
    		dxdt[55+ StateNum*i + StateNum*(segments*jCount)]     = dK_fi;
    		dxdt[56+ StateNum*i + StateNum*(segments*jCount)]     = dNa_fi;
      }

      if (jCount>0){
        // Vm
        dxdt[0+ StateNum*(segments-1)+ StateNum*(segments*(jCount-1))] = dxdt[0+ StateNum*(segments-1)+ StateNum*(segments*(jCount-1))]
                                                              + G_gap*(x[0+ StateNum*0+ StateNum*(segments*jCount)]
                                                              - x[0+ StateNum*(segments-1)+ StateNum*(segments*(jCount-1))]);
        dxdt[0+ StateNum*0+ StateNum*(segments*jCount)] = dxdt[0+ StateNum*0+ StateNum*(segments*jCount)]
                                              - G_gap*(x[0+ StateNum*0+ StateNum*(segments*jCount)]
                                              - x[0+ StateNum*(segments-1)+ StateNum*(segments*(jCount-1))]);
      }
    }
  }
};

void write_Shannon( const state_type &x , const double t )
{
    if (t>=kPrint){
      cout << t << '\t';
      for (int l = 0; l < StateNum*cells*segments; l++){
        cout << x[l] << '\t';
      }
      cout << endl;
      kPrint = kPrint+0.2;
    }
}

void write_Initials2( const state_type &x , const double t )
{
    if (t>=kPrint){
      cout << t << '\t';
      cout << endl;
      kPrint = kPrint+0.1;
    }
}

void write_Initials( const state_type &x , const double t )
{
    if (t==399.9){
      cout << t << '\t';
      for (int l = 0; l < StateNum*cells*segments; l++){
        cout << x[l] << '\t';
      }
      cout << endl;
    }
}

int main(int argc, char **argv)
{
  struct OpT T2;

  int kSarcRef = 0;
  int kG_gap = 0;
  for (int iCount = 0; iCount < cells; iCount++){
      std::string str = "SarcRef" + std::to_string(iCount+1) + ".txt";
      ifstream file(str);
      if(file.is_open())
      {
        for(int jCount = 0; jCount < segments; jCount++){
            file >> T2.SLsetAll[kSarcRef];
            kSarcRef++;
        }
      }
      std::string strGap = "G_gap" + std::to_string(iCount+1) + ".txt";
      ifstream fileGap(strGap);
      if(fileGap.is_open())
      {
        for(int jCount = 0; jCount < segments; jCount++){
            fileGap >> T2.GapAll[kG_gap];
            kG_gap++;
        }
      }
  }

  // read single cell matrices
  double y[segments*cells*StateNum];
  state_type x = {};


    /*int kInit = 0;
    int kSLset = 0;
      for (int i = 0; i < cells; i++){
        std::string str = "Init" + std::to_string(i+1) + ".txt";
        ifstream file(str);
        if(file.is_open())
        {
        for(int jCount = 0; jCount < segments; ++jCount)
        {
          for(int iCount = 0; iCount < StateNum; ++iCount)
          {
              file >> y[kInit];
              x[kInit] = y[kInit];
              kInit++;
          }
        }
      }
    }*/

    int kInit = 0;
    int kSLset = 0;
      std::string str = "Init.txt";
      ifstream file(str);
      if(file.is_open())
      {
      for(int jCount = 0; jCount < segments*StateNum*cells; ++jCount)
      {
          file >> y[kInit];
          x[kInit] = y[kInit];
          kInit++;
      }
    }

  Shannon ODEShannon(T2);
  for (int i = 1; i < 2; i++){
    kPrint = 0.0;
    cout << i<< endl;
    integrate( ODEShannon , x , 0.0 , 3500.0 , 0.1 , write_Shannon);
  }
  /*for (int i = 1; i < 2; i++){
    cout << i<< endl;
    integrate( ODEShannon , x , 0.0 , 2500.0 , 0.1, write_Shannon);
  }*/
}
