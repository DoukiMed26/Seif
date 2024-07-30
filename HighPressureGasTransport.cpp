/**
 *  @file HighPressureGasTransport.cpp
 *  Implementation file for class HighPressureGasTransport
 **/

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

// #include "cantera/thermo/PengRobinson.h"

#include "cantera/transport/HighPressureGasTransport.h"
#include "cantera/numerics/ctlapack.h"
#include "cantera/base/utilities.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/transport/MultiTransport.h"

using namespace std;

namespace Cantera
{

double HighPressureGasTransport::thermalConductivity()
{
    // Retrieve necessary properties from the state class
    double Tc = m_thermo->critTemperature();
    // cout << "Tc_lambda = " << Tc << "\n";
    double T = m_thermo->temperature();
    // cout << "T_lambda = " << T << "\n";
    double Vc = m_thermo->critVolume() * 1000;                  // Convert from [m^3/kmol] to [cm^3/mol] 
    // cout << "Vc_lambda = " << Vc << "\n";
    double M = m_thermo->meanMolecularWeight();                 // Convert from [kg/kmol] to [g/mol]    
    // cout << "M_lambda = " << M << "\n";
    double rho = ( m_thermo->density() / M ) * 0.001;           // Convert from kg/m^3 to mol/cm^3
    // cout << "rho_lambda = " << rho << "\n";
    double acentric = 0.0372;                                   // NIST for Nitrogen N2
    // double cv = m_thermo->cv_mass() * 0.2388 * M * 0.001;    // Convert heat capacity at constant volume unit from [J/kg.K] to [cal/mol.K]
                                                                // Restituisce cv usando la PR EoS. 
    size_t nsp = m_thermo->nSpecies();
    vector<double> cp_0_R(nsp);
    m_thermo->getCp_R(&cp_0_R[0]);
    double cp_0 = 0.;
    for (size_t i = 0; i < m_nsp; i++) {
        cp_0 = cp_0_R[i] * GasConstant;                         // cp_0 in [J/kmol/K] (prende l'unità di misura di R)
    }
    cout << "cp_0 = " << cp_0 << "\n";
    // double cp_0 = cp_0_r * GasConstant;                       // restituisce il Cp_0 in [J/kmol/K]
    double cv_0 = cp_0 - GasConstant;                            // restituisce il cv_0 in [J/kmol/K]
    double cv = cv_0 * 0.2388 * 0.001;                           // restituisce il cv_0 in [cal/mol.K]

    cout << "cv_0 = " << cv_0 << "\n";
    
    // Universal Gas Constantin Cantera: GasConstant in [J/kmol/K]
    
    double R = GasConstant/(4.186*1e3);                         // universal gas constat in [cal/mol.K]
    // cout << "R_lambda = " << R << "\n";
    double dipole = 0.;
    double kappa = 0.;

    const double A = 1.16145, B = 0.14874, C = 0.52487, D = 0.77320, E = 2.16178;
    const double F = 2.43787, G = -6.435e-4, H = 7.27371, S = 18.0323, W = -0.76830;

    // Calculate other required parameters
    double epsilon_over_k = Tc / 1.2593; // [K]
    double Tstar = T / epsilon_over_k;
    // cout << "Tstar_lambda = " << Tstar << "\n";

    double mu_r = 131.3 * dipole / sqrt(Vc * Tc);
    // cout << "mu_r_lambda = " << mu_r << "\n";
    double F_c = 1 - 0.2756 * acentric + 0.059035 * pow(mu_r, 4) + kappa; // [-]
    // cout << "F_c_lambda = " << F_c << "\n";

    // Calculate the reduced collision integral Omega* (=Omegast) in eq (2)  
    double Omegast = A * pow(Tstar, -B) + C * exp(-D * Tstar) + E * exp(-F * Tstar)
                           +G * pow(Tstar, B) * sin(S * pow(Tstar, W) - H); // [-]
    // cout << "Omegast_lambda = " << Omegast << "\n";

    // Calculate eta_0 (=eta0) in eq (6) 
    double eta0 = 4.0785e-5 * sqrt(M * T) / (pow(Vc, 2.0 / 3.0) * Omegast) * F_c; // [P]
    // cout << "eta0_lambda = " << eta0 << "\n";
    
    // Calculate Y term in eq (10) 
    double Y = (rho * Vc) / 6.0;
    // cout << "Y_lambda = " << Y << "\n";
   
    // Calculate G1 term in eq (10) 
    double G1 = (1.0 - 0.5 * Y) / pow(1 - Y, 3);
    // cout << "G1_lambda = " << G1 << "\n";
    
    // Constants for Chung viscosity model
    double b0[] = {2.41657, -0.50924, 6.61069, 14.54250, 0.79274, -5.8634, 81.171};
    double b1[] = {0.74824, -1.50936, 5.62073, -8.91387, 0.82019, 12.8005, 114.1580};
    double b2[] = {-0.91858, -49.99120, 64.7599, -5.63794, -0.69369, 9.58926, -60.841};
    double b3[] = {121.721, 69.9834, 27.0389, 74.3435, 6.31734, -65.5292, 466.7750};
    double MB[7];

    // Calculate the coefficients Bi in eq (13) 
    for (int j=0; j<7; j++ ){
       MB[j] = b0[j] + b1[j]*acentric + b2[j]*pow(mu_r, 4.0) + b3[j]*kappa; 
    } 

    // Calculate alpha (=alpha) in eq (9) 
    double alpha = (cv/R) - 3.0/2.0;

    // Calculate beta (=beta) in eq (9) 
    double beta = 0.7862 - 0.7109*acentric + 1.1368*pow(acentric, 2.0); 

    // Calculate Tr (=Tr) in eq (9)
    double Tr = T/Tc; 

    // Calculate Z (=Z) in eq (9)
    double Z = 2.0 + 10.5*pow(Tr, 2.0); 

    // Calculate Psi (=Psi) in eq (9) 
    double Psi = 1.0 + alpha*(0.215 + 0.28288*alpha - 1.061*beta + 0.26665*Z)/ 
                          (0.6366 + beta*Z + 1.061*alpha*beta); 

    // Calculate the thermal conductivity for dilute gases: lamda0 (=lamda0) in eq (9) 
    double lamda0 = 7.452*(eta0/M)*Psi; 
    cout << "lamda0 = " << lamda0 << " [W/(m.K)] " << "\n";

    // Calculate H2 (=H2) in eq (12)      
    double H2 = (MB[0]*(1.0-exp(-MB[3]*Y))/Y + MB[1]*G1*exp(MB[4]*Y)+MB[2]*G1)/
                      (MB[0]*MB[3] + MB[1] + MB[2]); 
    // cout << "H2 = " << H2 << "\n";

    // Calculate lamda_k(=lamdak) in eq (12)  
    double lamdak = lamda0*(1.0/H2 + MB[5]*Y);      
    // cout << "lamdak = " << lamdak << "\n";

    // Calculate lamda_p(=lamdap) in eq (12)
    double lamdap = (3.039e-4*pow(Tc/M, 1.0/2.0)/pow(Vc, 2.0/3.0))*MB[6]*pow(Y, 2.0)*H2*pow(Tr, 1.0/2.0);             
    // cout << "lamdap = " << lamdap << "\n";

    // Return the final thermal conductivity and change unit from [cal/cm.s.K] to [w/m.K] in eq (12) 
    return 418.6798 * (lamdak + lamdap);
    cout << "lambda = " << 418.6798 * (lamdak + lamdap) << " [W/(m.K)] " << "\n"; 
}

void HighPressureGasTransport::getThermalDiffCoeffs(double* const dt)
{
    // Method for MultiTransport class:
    // solveLMatrixEquation();
    // const double c = 1.6/GasConstant;
    // for (size_t k = 0; k < m_nsp; k++) {
    // dt[k] = c * m_mw[k] * m_molefracs[k] * m_a[k];
    // }
    throw NotImplementedError("HighPressureGasTransport::getThermalDiffCoeffs");
}

void HighPressureGasTransport::getBinaryDiffCoeffs(const size_t ld, double* const d)
{
    vector<double> PcP(5);
    size_t nsp = m_thermo->nSpecies();
    vector<double> molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    update_T();
    // Evaluate the binary diffusion coefficients from the polynomial fits.
    // This should perhaps be preceded by a check to see whether any of T, P, or
    //   C have changed.
    //if (!m_bindiff_ok) {
    updateDiff_T();
    //}
    if (ld < nsp) {
        throw CanteraError("HighPressureGasTransport::getBinaryDiffCoeffs",
                           "ld is too small");
    }
    double rp = 1.0/m_thermo->pressure();
    for (size_t i = 0; i < nsp; i++) {
        for (size_t j = 0; j < nsp; j++) {
            // Add an offset to avoid a condition where x_i and x_j both equal
            // zero (this would lead to Pr_ij = Inf):
            double x_i = std::max(Tiny, molefracs[i]);
            double x_j = std::max(Tiny, molefracs[j]);

            // Weight mole fractions of i and j so that X_i + X_j = 1.0:
            x_i = x_i/(x_i + x_j);
            x_j = x_j/(x_i + x_j);

            //Calculate Tr and Pr based on mole-fraction-weighted crit constants:
            double Tr_ij = m_temp/(x_i*Tcrit_i(i) + x_j*Tcrit_i(j));
            double Pr_ij = m_thermo->pressure()/(x_i*Pcrit_i(i) + x_j*Pcrit_i(j));

            double P_corr_ij;
            if (Pr_ij < 0.1) {
                // If pressure is low enough, no correction is needed:
                P_corr_ij = 1;
            }else {
                // Otherwise, calculate the parameters for Takahashi correlation
                // by interpolating on Pr_ij:
                P_corr_ij = setPcorr(Pr_ij, Tr_ij);

                // If the reduced temperature is too low, the correction factor
                // P_corr_ij will be < 0:
                if (P_corr_ij<0) {
                    P_corr_ij = Tiny;
                }
            }

            // Multiply the standard low-pressure binary diffusion coefficient
            // (m_bdiff) by the Takahashi correction factor P_corr_ij:
            d[ld*j + i] = P_corr_ij*rp * m_bdiff(i,j);
        }
    }
}

void HighPressureGasTransport::getMultiDiffCoeffs(const size_t ld, double* const d)
{
    // Not currently implemented.  m_Lmatrix inversion returns NaN.  Needs to be
    //   fixed.  --SCD - 2-28-2014
    throw NotImplementedError("HighPressureGasTransport:getMultiDiffCoeffs");
    // Calculate the multi-component Stefan-Maxwell diffusion coefficients,
    // based on the Takahashi-correlation-corrected binary diffusion coefficients.

    // update the mole fractions
    update_C();

    // update the binary diffusion coefficients
    update_T();
    updateThermal_T();

    // Correct the binary diffusion coefficients for high-pressure effects; this
    // is basically the same routine used in 'getBinaryDiffCoeffs,' above:
    size_t nsp = m_thermo->nSpecies();
    vector<double> molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);
    update_T();
    // Evaluate the binary diffusion coefficients from the polynomial fits -
    // this should perhaps be preceded by a check for changes in T, P, or C.
    updateDiff_T();

    if (ld < m_nsp) {
        throw CanteraError("HighPressureGasTransport::getMultiDiffCoeffs",
                           "ld is too small");
    }
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            // Add an offset to avoid a condition where x_i and x_j both equal
            //   zero (this would lead to Pr_ij = Inf):
            double x_i = std::max(Tiny, molefracs[i]);
            double x_j = std::max(Tiny, molefracs[j]);
            x_i = x_i/(x_i+x_j);
            x_j = x_j/(x_i+x_j);
            double Tr_ij = m_temp/(x_i*Tcrit_i(i) + x_j*Tcrit_i(j));
            double Pr_ij = m_thermo->pressure()/(x_i*Pcrit_i(i) + x_j*Pcrit_i(j));

            double P_corr_ij;
            if (Pr_ij < 0.1) {
                P_corr_ij = 1;
            }else {
                P_corr_ij = setPcorr(Pr_ij, Tr_ij);
                if (P_corr_ij<0) {
                    P_corr_ij = Tiny;
                }
            }

            m_bdiff(i,j) *= P_corr_ij;
        }
    }
    m_bindiff_ok = false; // m_bdiff is overwritten by the above routine.

    // Having corrected m_bdiff for pressure and concentration effects, the
    //    routine now proceeds the same as in the low-pressure case:

    // evaluate L0000 if the temperature or concentrations have
    // changed since it was last evaluated.
    if (!m_l0000_ok) {
        eval_L0000(molefracs.data());
    }

    // invert L00,00
    int ierr = invert(m_Lmatrix, m_nsp);
    if (ierr != 0) {
        throw CanteraError("HighPressureGasTransport::getMultiDiffCoeffs",
                           "invert returned ierr = {}", ierr);
    }
    m_l0000_ok = false; // matrix is overwritten by inverse
    m_lmatrix_soln_ok = false;

    double prefactor = 16.0 * m_temp
        *m_thermo->meanMolecularWeight()/(25.0*m_thermo->pressure());

    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            double c = prefactor/m_mw[j];
            d[ld*j + i] = c*molefracs[i]*(m_Lmatrix(i,j) - m_Lmatrix(i,i));
        }
    }
}

double HighPressureGasTransport::viscosity() 
{
    // Retrieve necessary properties from the state class
    double Tc_mix = m_thermo->critTemperature();
    // cout << "Tc_visc = " << Tc_mix << "\n";
    double acentric = 0.0372;   // NIST
    double T = m_thermo->temperature();
    // cout << "T_visc = " << T << "\n";
    double M_gmol = m_thermo->meanMolecularWeight(); // Convert from [kg/kmol] to [g/mol]
    // cout << "M_visc = " << M_gmol << "\n";
    double rho_molcm3 = ( m_thermo->density() / M_gmol ) * 0.001; // Convert from kg/m^3 to mol/cm^3
    // cout << "rho_visc" << rho_molcm3 << "\n";     
    double kappa = 0.;
    // double mu_r = 0.; 
    double Vc_cm3mol = m_thermo->critVolume() * 1000; // Convert from [m^3/kmol] to [cm^3/mol] 
    // cout << "Vc_visc" << Vc_cm3mol << "\n";     
    double dipole = 0.;

    // Constants for Chung viscosity model
    double a0[] = {0, 6.32402, 0.12102e-2, 5.28346, 6.62263, 19.74540, -1.89992, 24.27450, 0.79716, -0.23816, 0.68629e-1};
    double a1[] = {0, 50.41190, -0.11536e-2, 254.20900, 38.09570, 7.63034, -12.53670, 3.44945, 1.11764, 0.67695e-1, 0.34793};
    double a2[] = {0, -51.68010, -0.62571e-2, -168.48100, -8.46414, -14.35440, 4.98529, -11.29130, 0.12348e-1, -0.81630, 0.59256};
    double a3[] = {0, 1189.02000, 0.37283e-1, 3898.27000, 31.41780, 31.52670, -18.15070, 69.34660, -4.11661, 4.02528, -0.72663};
    double A[11];

    double mu_r = 131.3 * dipole / sqrt(Vc_cm3mol * Tc_mix);  // !!!!!!!!!!!!!!!!
    // cout << "mu_r_visc = " << mu_r << "\n";
    
    double F_c = 1 - 0.2756 * acentric + 0.059035 * pow(mu_r, 4) + kappa; // [-]
    // cout << "F_c_visc = " << F_c << "\n";

    // Calculate coefficients A
    for (int i = 1; i <= 10; ++i) {
        A[i] = a0[i] + a1[i] * acentric + a2[i] * pow(mu_r, 4) + a3[i] * kappa; // Assuming kappa is 0
    }

    // Calculate other required parameters
    double epsilon_over_k = Tc_mix / 1.2593; // [K]
    double Tstar = T / epsilon_over_k;
    // cout << "Tstar_visc = " << Tstar << "\n";
    
    double Y = rho_molcm3 * Vc_cm3mol / 6.0;
    // cout << "Y_visc = " << Y << "\n";
    
    double Omega_2_2 = 1.16145 * pow(Tstar, -0.14874) + 0.52487 * exp(-0.77320 * Tstar) + 2.16178 * exp(-2.43787 * Tstar)
                           - 6.435e-4 * pow(Tstar, 0.14874) * sin(18.0323 * pow(Tstar, -0.76830) - 7.27371); // [-]
    // cout << "Omega_lambda = " << Omega_2_2 << "\n";

    double eta0_P = 4.0785e-5 * sqrt(M_gmol * T) / (pow(Vc_cm3mol, 2.0 / 3.0) * Omega_2_2) * F_c; // [P]
    cout << "viscosity_0 = " << eta0_P/10 << " [Pa.s] " << "\n";

    // Calculate viscosity using Chung model
    double G_1 = (1.0 - 0.5 * Y) / pow(1 - Y, 3);
    // cout << "G1_visc = " << G_1 << "\n";

    double G_2 = (A[1] * (1 - exp(-A[4] * Y)) / Y + A[2] * G_1 * exp(A[5] * Y) + A[3] * G_1) / (A[1] * A[4] + A[2] + A[3]);
    double eta_k_P = eta0_P * (1 / G_2 + A[6] * Y); // [P]

    double eta_p_P = (36.344e-6 * sqrt(M_gmol * Tc_mix) / pow(Vc_cm3mol, 2.0 / 3.0)) * A[7] * pow(Y, 2) * G_2
                         * exp(A[8] + A[9] / Tstar + A[10] / pow(Tstar, 2)); // [P]

    return (eta_k_P + eta_p_P)/10; // [Pa*s]
    cout << "viscosity = " << (eta_k_P + eta_p_P)/10 << " [Pa.s] " << "\n";
}

// Pure species critical properties - Tc, Pc, Vc, Zc:
double HighPressureGasTransport::Tcrit_i(size_t i)
{
    // Store current molefracs and set temp molefrac of species i to 1.0:
    vector<double> molefracs = store(i, m_thermo->nSpecies());

    double tc = m_thermo->critTemperature();
    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return tc;
}

double HighPressureGasTransport::Pcrit_i(size_t i)
{
    // Store current molefracs and set temp molefrac of species i to 1.0:
    vector<double> molefracs = store(i, m_thermo->nSpecies());

    double pc = m_thermo->critPressure();
    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return pc;
}

double HighPressureGasTransport::Vcrit_i(size_t i)
{
    // Store current molefracs and set temp molefrac of species i to 1.0:
    vector<double> molefracs = store(i, m_thermo->nSpecies());

    double vc = m_thermo->critVolume();
    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return vc;
}

double HighPressureGasTransport::Zcrit_i(size_t i)
{
    // Store current molefracs and set temp molefrac of species i to 1.0:
    vector<double> molefracs = store(i, m_thermo->nSpecies());

    double zc = m_thermo->critCompressibility();
    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return zc;
}

vector<double> HighPressureGasTransport::store(size_t i, size_t nsp)
{
    vector<double> molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);
    vector<double> mf_temp(nsp, 0.0);
    mf_temp[i] = 1;
    m_thermo->setMoleFractions(&mf_temp[0]);
    return molefracs;
}

// Calculates quantum correction term for a species based on Tr and MW, used in
//   viscosity calculation:
double HighPressureGasTransport::FQ_i(double Q, double Tr, double MW)
{
    return 1.22*pow(Q,0.15)*(1 + 0.00385*pow(pow(Tr - 12.,2.),1./MW)
                             *fabs(Tr-12)/(Tr-12));
}

// Set value of parameter values for Takahashi correlation, by interpolating
//   table of constants vs. Pr:
double HighPressureGasTransport::setPcorr(double Pr, double Tr)
{
    const static double Pr_lookup[17] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0,
        1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0};
    const static double DP_Rt_lookup[17] = {1.01, 1.01, 1.01, 1.01, 1.01, 1.01,
        1.01, 1.02, 1.02, 1.02, 1.02, 1.03, 1.03, 1.04, 1.05, 1.06, 1.07};
    const static double A_ij_lookup[17] = {0.038042, 0.067433, 0.098317,
        0.137610, 0.175081, 0.216376, 0.314051, 0.385736, 0.514553, 0.599184,
        0.557725, 0.593007, 0.696001, 0.790770, 0.502100, 0.837452, 0.890390};
    const static double B_ij_lookup[17] = {1.52267, 2.16794, 2.42910, 2.77605,
        2.98256, 3.11384, 3.50264, 3.07773, 3.54744, 3.61216, 3.41882, 3.18415,
        3.37660, 3.27984, 3.39031, 3.23513, 3.13001};
    const static double C_ij_lookup[17] = {0., 0., 0., 0., 0., 0., 0., 0.141211,
        0.278407, 0.372683, 0.504894, 0.678469, 0.665702, 0., 0.602907, 0., 0.};
    const static double E_ij_lookup[17] = {1., 1., 1., 1., 1., 1., 1., 13.45454,
        14., 10.00900, 8.57519, 10.37483, 11.21674, 1., 6.19043, 1., 1.};

    // Interpolate Pr vs. those used in Takahashi table:
    int Pr_i = 0;
    double frac = 0.;

    if (Pr < 0.1) {
        frac = (Pr - Pr_lookup[0])/(Pr_lookup[1] - Pr_lookup[0]);
    } else {
        for (int j = 1; j < 17; j++) {
            if (Pr_lookup[j] > Pr) {
                frac = (Pr - Pr_lookup[j-1])/(Pr_lookup[j] - Pr_lookup[j-1]);
                break;
            }
            Pr_i++;
        }
    }
    // If Pr is greater than the greatest value used by Takahashi (5.0), use the
    //   final table value.  Should eventually add in an extrapolation:
    if (Pr_i == 17) {
        frac = 1.0;
    }

    double P_corr_1 = DP_Rt_lookup[Pr_i]*(1.0 - A_ij_lookup[Pr_i]
        *pow(Tr,-B_ij_lookup[Pr_i]))*(1-C_ij_lookup[Pr_i]
        *pow(Tr,-E_ij_lookup[Pr_i]));
    double P_corr_2 = DP_Rt_lookup[Pr_i+1]*(1.0 - A_ij_lookup[Pr_i+1]
        *pow(Tr,-B_ij_lookup[Pr_i+1]))*(1-C_ij_lookup[Pr_i+1]
        *pow(Tr,-E_ij_lookup[Pr_i+1]));
    return P_corr_1*(1.0-frac) + P_corr_2*frac;
}

}
