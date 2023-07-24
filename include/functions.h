#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "includes.h"


//############################################
//##                                        ##
//##     Adapted FitzHugh-Nagumo (AFHN)     ##
//##                                        ##
//############################################
#if defined(AFHN)
double reactionV(double v, double w)
{
    return (1.0 / (Cm * chi)) * ((-G * v * (1.0 - (v / vth)) * (1.0 - (v / vp))) + (-eta1 * v * w));
}

double reactionW(double v, double w)
{
    return eta2 * ((v / vp) - (eta3 * w));
}
#endif // AFHN



//###########################################
//##                                       ##
//##     ten Tusscher 2006 model (TT2)     ##
//##                                       ##
//###########################################
#if defined(TT2)
#if defined(EPI) || defined(M) || defined(ENDO)
/*---------------------------------------------------------------------------------------------------------------------------------------------------
Functions for ten Tusscher model 2006 (https://journals.physiology.org/doi/full/10.1152/ajpheart.00109.2006)
from https://tbb.bio.uu.nl/khwjtuss/SourceCodes/HVM2/Source/Main.cc - ten Tusscher code
and https://github.com/rsachetto/MonoAlg3D_C/blob/master/src/models_library/ten_tusscher/ten_tusscher_2006_RS_CPU.c - Sachetto MonoAlg3D
-----------------------------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------
Currents functions
----------------------*/
// Reversal potentials for Na+, K+ and Ca++
double E_Na(double Na_i)
{
    return RTONF * log(Na_o / Na_i);
}
double E_K(double K_i)
{
    return RTONF * log(K_o / K_i);
}
double E_Ca(double Ca_i)
{
    return 0.5 * RTONF * log(Ca_o / Ca_i);
}

// Reversal potential for Ks
double E_Ks(double K_i, double Na_i)
{
    return RTONF * log((K_o + p_KNa * Na_o) / (K_i + p_KNa * Na_i));
}

// Fast sodium (Na+) current
double I_Na(double V, double m, double h, double j, double Na_i)
{
    return G_Na * (m*m*m) * h * j * (V - E_Na(Na_i));
}
double m_inf(double V)
{
    return 1.0 / ((1.0 + exp((-56.86 - V) / 9.03))*(1.0 + exp((-56.86 - V) / 9.03)));
}
double alpha_m(double V)
{
    return 1.0 / (1.0 + exp((-60.0 - V) / 5.0));
}
double beta_m(double V)
{
    return (0.1 / (1.0 + exp((V + 35.0) / 5.0))) + (0.1 / (1.0 + exp((V - 50.0) / 200.0)));
}
double tau_m(double V)
{
    return alpha_m(V) * beta_m(V);
}
double h_inf(double V)
{
    return 1.0 / ((1.0 + exp((V + 71.55) / 7.43))*(1.0 + exp((V + 71.55) / 7.43)));
}
double alpha_h(double V)
{
    if (V >= -40.0)
    {
        return 0.0;
    }
    else
    {
        return 0.057 * exp(-(80.0 + V) / 6.8);
    }
}
double beta_h(double V)
{
    if (V >= -40.0)
    {
        return 0.77 / (0.13 * (1.0 + exp((V + 10.66) / (-11.1))));
    }
    else
    {
        return 2.7 * exp(0.079 * V) + 3.1e5 * exp(0.3485 * V);
    }
}
double tau_h(double V)
{
    return 1.0 / (alpha_h(V) + beta_h(V));
}
double j_inf(double V)
{
    return 1.0 / ((1.0 + exp((V + 71.55) / 7.43))*(1.0 + exp((V + 71.55) / 7.43)));
}
double alpha_j(double V)
{
    if (V >= -40.0)
    {
        return 0.0;
    }
    else
    {
        return ((-25428.0 * exp(0.2444 * V) - (6.948e-6 * exp((-0.04391) * V))) * (V + 37.78)) / (1.0 + exp(0.311 * (V + 79.23)));
    }
}
double beta_j(double V)
{
    if (V >= -40.0)
    {
        return (0.6 * exp(0.057 * V)) / (1.0 + exp(-0.1 * (V + 32.0)));
    }
    else
    {
        return (0.02424 * exp(-0.01052 * V)) / (1.0 + exp(-0.1378 * (V + 40.14)));
    }
}
double tau_j(double V)
{
    return 1.0 / (alpha_j(V) + beta_j(V));
}

// L-type Ca2+ current
double I_CaL(double V, double d, double f, double f2, double fCass, double Ca_SS)   // !!!
{
    if (V < 15.0 - 1.0e-5)
    {
        return G_CaL * d * f * f2 * fCass * 4.0 * (V - 15.0) * (F*F) * (0.25 * Ca_SS * exp(2 * (V - 15.0) * FONRT) - Ca_o) / (R * T * (exp(2.0 * (V - 15.0) * FONRT) - 1.0));
    }
    else if (V > 15.0 + 1.0e-5)
    {
        return G_CaL * d * f * f2 * fCass * 2.0 * F * (0.25 * Ca_SS - Ca_o);
    }
}
double d_inf(double V)
{
    return 1.0 / (1.0 + exp((-8.0 - V) / 7.5));
}
double alpha_d(double V)
{
    return (1.4 / (1.0 + exp((-35.0 - V) / 13.0))) + 0.25;
}
double beta_d(double V)
{
    return 1.4 / (1.0 + exp((V + 5.0) / 5.0));
}
double gamma_d(double V)
{
    return 1.0 / (1.0 + exp((50.0 - V) / 20.0));
}
double tau_d(double V)
{
    return alpha_d(V) * beta_d(V) + gamma_d(V);
}
double f_inf(double V)
{
    return 1.0 / (1.0 + exp((V + 20.0) / 7.0));
}
double alpha_f(double V)
{
    return 1102.5 * exp(-(((V + 27.0)*(V + 27.0))) / 225.0);
}
double beta_f(double V)
{
    return 200.0 / (1.0 + exp((13.0 - V) / 10.0));
}
double gamma_f(double V)
{
    return 180.0 / (1.0 + exp((V + 30.0) / 10.0)) + 20.0;
}
double tau_f(double V)
{
    return alpha_f(V) + beta_f(V) + gamma_f(V);
}
double f2_inf(double V)
{
    return 0.67 / (1.0 + exp((V + 35.0) / 7.0)) + 0.33;
}
double alpha_f2(double V)   // !!!
{
    return 562.0 * exp(-(((V + 27.0)*(V + 27.0))) / 240.0);
}
double beta_f2(double V)
{
    return 31.0 / (1.0 + exp((25.0 - V) / 10.0));
}
double gamma_f2(double V)   // !!!
{
    return 80.0 / (1.0 + exp((V + 30.0) / 10.0));
}
double tau_f2(double V)
{
    return alpha_f2(V) + beta_f2(V) + gamma_f2(V);
}
double fCass_inf(double Ca_SS)
{
    return 0.6 / (1.0 + ((Ca_SS / 0.05)*(Ca_SS / 0.05))) + 0.4;
}
double tau_fCass(double Ca_SS)
{
    return 80.0 / (1.0 + ((Ca_SS / 0.05)*(Ca_SS / 0.05))) + 2.0;
}

// Transient outward current
double I_to(double V, double r, double s, double K_i)
{
    return G_to * r * s * (V - E_K(K_i));
}
double r_inf(double V)
{
    return 1.0 / (1.0 + exp((20.0 - V) / 6.0));
}
double tau_r(double V)
{
    return 9.5 * exp(-(((V + 40.0)*(V + 40.0))) / 1800.0) + 0.8;
}
#if defined(EPI) || defined(M)  // for epicardial and M cells
double s_inf(double V)
{
    return 1.0 / (1.0 + exp((V + 20.0) / 5.0));
}
double tau_s(double V)
{
    return 85.0 * exp(-(((V + 45.0)*(V + 45.0))) / 320.0) + 5.0 / (1.0 + exp((V - 20.0) / 5.0)) + 3.0;
}
#endif
#ifdef ENDO  // for endocardial cells
double s_inf(double V)
{
    return 1.0 / (1.0 + exp((V + 28.0) / 5.0));
}
double tau_s(double V)
{
    return 1000.0 * exp(-(((V + 67.0)*(V + 67.0))) / 1000.0) + 8.0;
}
#endif

// Slow delayed rectifier current
double I_Ks(double V, double X_s, double K_i, double Na_i)
{
    return G_Ks * (X_s*X_s) * (V - E_Ks(K_i, Na_i));
}
double x_s_inf(double V)
{
    return 1.0 / (1.0 + exp((-5.0 - V) / 14.0));
}
double alpha_x_s(double V)
{
    return 1400.0 / sqrt(1.0 + exp((5.0 - V) / 6.0));
}
double beta_x_s(double V)
{
    return 1.0 / (1.0 + exp((V - 35.0) / 15.0));
}
double tau_x_s(double V)
{
    return alpha_x_s(V) * beta_x_s(V) + 80.0;
}

// Rapid delayed rectifier current
double I_Kr(double V, double X_r1, double X_r2, double K_i)
{
    return G_Kr * sqrt(K_o / 5.4) * X_r1 * X_r2 * (V - E_K(K_i));
}
double x_r1_inf(double V)
{
    return 1.0 / (1.0 + exp((-26.0 - V) / 7.0));
}
double alpha_x_r1(double V)
{
    return 450.0 / (1.0 + exp((-45.0 - V) / 10.0));
}
double beta_x_r1(double V)
{
    return 6.0 / (1.0 + exp((V + 30.0) / 11.5));
}
double tau_x_r1(double V)
{
    return alpha_x_r1(V) * beta_x_r1(V);
}
double x_r2_inf(double V)
{
    return 1.0 / (1.0 + exp((V + 88.0) / 24.0));
}
double alpha_x_r2(double V)
{
    return 3.0 / (1.0 + exp((-60.0 - V) / 20.0));
}
double beta_x_r2(double V)
{
    return 1.12 / (1.0 + exp((V - 60.0) / 20.0));
}
double tau_x_r2(double V)
{
    return alpha_x_r2(V) * beta_x_r2(V);
}

// Inward rectifier K+ current
double alpha_K1(double V, double K_i)
{
    return 0.1 / (1.0 + exp(0.06 * (V - E_K(K_i) - 200.0)));
}
double beta_K1(double V, double K_i)
{
    return (3.0 * exp(0.0002 * (V - E_K(K_i) + 100.0)) + exp(0.1 * (V - E_K(K_i) - 10.0))) / (1.0 + exp(-0.5 * (V - E_K(K_i))));
}
double x_K1_inf(double V, double K_i)
{
    return alpha_K1(V, K_i) / (alpha_K1(V, K_i) + beta_K1(V, K_i));
}
double I_K1(double V, double K_i)
{
    return G_K1 * x_K1_inf(V, K_i) * (V - E_K(K_i));
}

// Na+/Ca++ exchanger current
double I_NaCa(double V, double Na_i, double Ca_i)   // !!!
{
    return (k_NaCa * ((exp((gamma_I_NaCa * V * FONRT)) * (Na_i*Na_i*Na_i) * Ca_o) - (exp(((gamma_I_NaCa - 1.0) * V * FONRT)) * (Na_o*Na_o*Na_o) * Ca_i * alpha))) / (((K_mNa_i*K_mNa_i*K_mNa_i) + (Na_o*Na_o*Na_o)) * (K_mCa + Ca_o) * (1.0 + (k_sat * exp(((gamma_I_NaCa) * V * FONRT)))));
}

// Na+/K+ pump current
double I_NaK(double V, double Na_i) // !!!
{
    return ((((p_KNa * K_o) / (K_o + K_mK)) * Na_i) / (Na_i + K_mNa)) / (1.0 + (0.1245 * exp(((-0.1) * V * FONRT))) + (0.0353 * exp(((-V) * FONRT))));
}

// I_pCa
double I_pCa(double V, double Ca_i)
{
    return (G_pCa * Ca_i) / (K_pCa + Ca_i);
}

// I_pK
double I_pK(double V, double K_i)
{
    return (G_pK * (V - E_K(K_i))) / (1.0 + exp((25.0 - V) / 5.98));
}

// Background currents
double I_bNa(double V, double Na_i)
{
    return G_bNa * (V - E_Na(Na_i));
}
double I_bCa(double V, double Ca_i)
{
    return G_bCa * (V - E_Ca(Ca_i));
}

// Calcium dynamics
double I_leak(double Ca_SR, double Ca_i)
{
    return V_leak * (Ca_SR - Ca_i);
}
double I_up(double Ca_i)
{
    return V_maxup / (1.0 + ((K_up*K_up) / (Ca_i*Ca_i)));
}
double k_casr(double Ca_SR)
{
    return max_SR - ((max_SR - min_SR) / (1.0 + ((EC / Ca_SR)*(EC / Ca_SR))));
}
double k1(double Ca_SR)
{
    return k1_prime / k_casr(Ca_SR);
}
double O(double Ca_SR, double Ca_SS, double R_prime)
{
    return (k1(Ca_SR) * (Ca_SS*Ca_SS) * R_prime) / (k3 + (k1(Ca_SR) * (Ca_SS*Ca_SS)));
}
double I_rel(double Ca_SR, double Ca_SS, double R_prime)
{
    return V_rel * O(Ca_SR, Ca_SS, R_prime) * (Ca_SR - Ca_SS);
}
double I_xfer(double Ca_SS, double Ca_i)
{
    return V_xfer * (Ca_SS - Ca_i);
}
double k2(double Ca_SR)
{
    return k2_prime * k_casr(Ca_SR);
}
double Ca_ibufc(double Ca_i)    // !!!
{
    return 1.0 / (1.0 + ((Buf_C * K_bufc) / ((Ca_i + K_bufc)*(Ca_i + K_bufc))));
}
double Ca_srbufsr(double Ca_SR) // !!!
{
    return 1.0 / (1.0 + ((Buf_SR * K_bufsr) / ((Ca_SR + K_bufsr)*(Ca_SR + K_bufsr))));
}
double Ca_ssbufss(double Ca_SS) // !!!
{
    return 1.0 / (1.0 + ((Buf_SS * K_bufss) / ((Ca_SS + K_bufss)*(Ca_SS + K_bufss))));
}


/*-----------------------------------------------------
Differential equations for each variable
-----------------------------------------------------*/
double Itotal(double I_stim, double V, double m, double h, double j, double Na_i, double K_i, double r, double s, double X_r1, double X_r2, double X_s, double d, double f, double f2, double fCass, double Ca_SS, double Ca_i)
{
    double VmENa = V - E_Na(Na_i);
    double VmEK = V - E_K(K_i);

    double INa = G_Na * (m*m*m) * h * j * VmENa;
    double IbNa = G_bNa * VmENa;
    double IK1 = G_K1 * x_K1_inf(V, K_i) * VmEK;
    double Ito = G_to * r * s * VmEK;
    double IKr = G_Kr * sqrt(K_o / 5.4) * X_r1 * X_r2 * VmEK;
    double IKs = I_Ks(V, X_s, K_i, Na_i);
    double ICaL = I_CaL(V, d, f, f2, fCass, Ca_SS);
    double INaK = I_NaK(V, Na_i);
    double INaCa = I_NaCa(V, Na_i, Ca_i);
    double IpCa = I_pCa(V, Ca_i);
    double IpK = (G_pK * VmEK) / (1.0 + exp((25.0 - V) / 5.98));
    double IbCa = I_bCa(V, Ca_i);

    return I_stim + INa + IbNa + IK1 + Ito + IKr + IKs + ICaL + INaK + INaCa + IpCa + IpK + IbCa;
}

double dRprimedt(double Ca_SS, double R_prime)
{
    return ((-k2(Ca_SS)) * Ca_SS * R_prime) + (k4 * (1.0 - R_prime));
}

double dCaidt(double Ca_i, double Ca_SR, double Ca_SS, double V, double Na_i)
{
    return Ca_ibufc(Ca_i) * (((((I_leak(Ca_SR, Ca_i) - I_up(Ca_i)) * V_SR) / V_C) + I_xfer(Ca_SS, Ca_i)) - ((((I_bCa(V, Ca_i) + I_pCa(V, Ca_i)) - (2.0 * I_NaCa(V, Na_i, Ca_i))) * Cm) / (2.0 * V_C * F)));
}

double dCaSRdt(double Ca_SR, double Ca_i, double Ca_SS, double R_prime)
{
    return Ca_srbufsr(Ca_SR) * (I_up(Ca_i) - (I_rel(Ca_SR, Ca_SS, R_prime) + I_leak(Ca_SR, Ca_i)));
}

double dCaSSdt(double Ca_SS, double V, double d, double f, double f2, double fCass, double Ca_SR, double R_prime, double Ca_i)
{
    return Ca_ssbufss(Ca_SS) * (((((-I_CaL(V, d, f, f2, fCass, Ca_SS)) * Cm) / (2.0 * V_SS * F)) + ((I_rel(Ca_SR, Ca_SS, R_prime) * V_SR) / V_SS)) - ((I_xfer(Ca_SS, Ca_i) * V_C) / V_SS));
}

double dNaidt(double V, double m, double h, double j, double Na_i, double Ca_i)
{
    return ((-(I_Na(V, m, h, j, Na_i) + I_bNa(V, Na_i) + (3.0 * I_NaK(V, Na_i)) + (3.0 * I_NaCa(V, Na_i, Ca_i)))) / (V_C * F)) * Cm;
}

double dKidt(double I_stim, double V, double K_i, double r, double s, double X_r1, double X_r2, double X_s, double Na_i)
{
    return ((-((I_stim + I_K1(V, K_i) + I_to(V, r, s, K_i) + I_Kr(V, X_r1, X_r2, K_i) + I_Ks(V, X_s, K_i, Na_i) + I_pK(V, K_i)) - (2.0 * I_NaK(V, Na_i)))) / (V_C * F)) * Cm;
}


/*-----------------------------------------------------
Differential equations for each variable
-----------------------------------------------------*/
double updateXr1(double X_r1, double V, double dt)
{
    double xr1inf = x_r1_inf(V);
    return xr1inf - (xr1inf - X_r1) * exp(-dt / tau_x_r1(V));
}

double updateXr2(double X_r2, double V, double dt)
{
    double xr2inf = x_r2_inf(V);
    return xr2inf - (xr2inf - X_r2) * exp(-dt / tau_x_r2(V));
}

double updateXs(double X_s, double V, double dt)
{
    double xsinf = x_s_inf(V);
    return xsinf - (xsinf - X_s) * exp(-dt / tau_x_s(V));
}

double updater(double r, double V, double dt)
{
    double rinf = r_inf(V);
    return rinf - (rinf - r) * exp(-dt / tau_r(V));
}

double updates(double s, double V, double dt)
{
    double sinf = s_inf(V);
    return sinf - (sinf - s) * exp(-dt / tau_s(V));
}

double updatem(double m, double V, double dt)
{
    double minf = m_inf(V);
    return minf - (minf - m) * exp(-dt / tau_m(V));
}

double updateh(double h, double V, double dt)
{
    double hinf = h_inf(V);
    return hinf - (hinf - h) * exp(-dt / tau_h(V));
}

double updatej(double j, double V, double dt)
{
    double jinf = j_inf(V);
    return jinf - (jinf - j) * exp(-dt / tau_j(V));
}

double updated(double d, double V, double dt)
{
    double dinf = d_inf(V);
    return dinf - (dinf - d) * exp(-dt / tau_d(V));
}

double updatef(double f, double V, double dt)
{
    double finf = f_inf(V);
    return finf - (finf - f) * exp(-dt / tau_f(V));
}

double updatef2(double f2, double V, double dt)
{
    double f2inf = f2_inf(V);
    return f2inf - (f2inf - f2) * exp(-dt / tau_f2(V));
}

double updatefCass(double fCass, double V, double dt)
{
    double fCassinf = fCass_inf(V);
    return fCassinf - (fCassinf - fCass) * exp(-dt / tau_fCass(V));
}
#endif  // EPI || M || ENDO
#endif  // TT2


#endif // FUNCTIONS_H