// Steady 1D Nozzle using SIMPLE and SIMPLER, upwind differencing
// See Versteeg, An Introduction to Computational Fluid Dynamics, Example 6.2
// Travis Burrows

#include <cmath>
#include <cstdlib>
#include <math.h>
#include <omp.h>
#include <string.h>

#include <array>
#include <iostream>
#include <vector>

//#define SR

double tempy;

//  Parameters
const int N = 10;            //  Number of points
const double speedoflight = 3e8;            //  Number of points
const double THRESH = 1e-3;  //  Outer Iteration L2 norm threshold
// const double THRESH = 1E-2;  //  Outer Iteration L2 norm threshold
const double OMEGAu = 0.05;  //  Momentum relaxation coefficient
const double OMEGAp = 0.01;  //  Pressure relaxation coefficient
// const double OMEGAu = 0.02;  //  Momentum relaxation coefficient
// const double OMEGAp = 0.02;  //  Pressure relaxation coefficient
const double MAXITER = 1E4;  //  Maximum iterations
const int DEBUG = 0;         //  Print extra information
const int printint = 10;     //  Interval to print convergence information

//  Global Constants
const double L = 2.0;    //  Length of domain (1D)
const double Ai = 100;   //  Initial area in meters squared
const double Af = 1;   //  Final area in meters squared
const double rho = 1.0;  //  Fluid density
const double Pi = 1e17;    //  Inlet pressure
const double Po = 0;     //  Outlet pressure
// const double Mi = 1000000.0;   //  initial mass flow rate
const double Mi = pow(2 * rho * (Pi - Po) * Ai * Af / abs(pow(Af, 2) - pow(Ai, 2)), 0.5);   //  initial mass flow rate
const double OmOMEGAu = 1.0 - OMEGAu;
const double OmOMEGAp = 1.0 - OMEGAp;
const int Nm1 = N - 1;
const int Nm2 = N - 2;
const double dx = L / Nm1;  //  delta x
const double THRESHinner = THRESH;
    //THRESH / 150.0;               //  Inner iteration residual threshold
const double Mexact = sqrt(0.2);  //  Exact solution mass flow rate
const double De = 0.0, Dw = 0.0;  //  Diffusion coefficients

const double adiabatic_constant = 4.0/3;
const double enthalpy_coeff = adiabatic_constant/(adiabatic_constant-1);

//  Data types
typedef std::array<double, N> Pvec;    //  1D vector of doubles for pressure
typedef std::array<double, Nm1> Uvec;  //  1D vector of doubles for velocity

//  Function Declarations
template <std::size_t SIZE>
void printMatrix(const std::string &name, const std::array<double, SIZE> &vec);
double P2(const double &value);
double Pexact(const double &area);
double Uexact(const double &area);
void buildMomentumCoeffs(const Pvec &p, const Pvec &areaP, const Uvec &uPrev,
                         const Uvec &areaU, Uvec &Su, Uvec &aPu, Uvec &aWu,
                         Uvec &aEu, Uvec &aPuinv, Uvec &d, const double rho);
void buildPressureCoeffs(Pvec &aWp, Pvec &aEp, Pvec &aPp, Pvec &aPpinv,
                         Pvec &bPrimep, const Uvec &d, const Uvec &usource,
                         const Uvec &areaU, const double rho);
void solnError(double &perror, double &uerror, const Pvec &p,
               const Pvec &pexact, const Uvec &u, const Uvec &uexact);
void throwError(const std::string &message);
void checkparams();

int main(int argc, char** argv) {

    //  Check for valid parameters
    checkparams();

    //  Declare variables
    Pvec xP{}, areaP{}, p{}, pPrime{}, aWp{}, aEp{}, aPp{}, aPpinv{}, bPrimep{},
        pPrev{}, Pexsoln{};
    Uvec xU{}, areaU{}, u{}, uStar{}, d{}, mfr{}, Su{}, uPrev{}, Uexsoln{},
        uHat{}, aPu{}, aEu{}, aWu{}, aPuinv{};
    double dif{}, UP{}, UE{}, UW{}, temp{}, PP{}, totaldif{}, totaldif0{},
        start{}, stop{}, uError{}, pError{};
    int Mcount{}, Pcount{}, Pprimecount{};

    //  Declare AMGCL-specific variables

    double rhs{};

    //  Initialize OpenMP with maximum threads
    const int maxthreads = omp_get_max_threads() - 1;
    omp_set_num_threads(maxthreads);

    //  Print some parameters
    printf("Threads: %d\n", maxthreads);
    printf("dx: %.2e\n", dx);
    printf("N: %d\n", N);
    std::cout << std::endl;

    //  Initialize Variables
    const double Aslope = (Af - Ai) / L;
    const double Pslope = (Po - Pi) / L;
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        xP[i] = i * dx;
        //  area at P nodes
        areaP[i] = Ai + Aslope * i * dx;
        //  exact pressure solution
        Pexsoln[i] = Pexact(areaP[i]);
        //  initial guess of pressure
        p[i] = Pi + Pslope * i * dx;
        if (i < Nm1) {
            xU[i] = 0.5 * dx + i * dx;
            //  area at U nodes
            areaU[i] = Ai + Aslope * 0.5 * dx + Aslope * i * dx;
            //  exact u solution
            Uexsoln[i] = Uexact(areaU[i]);
            //  initial guess of velocity
            u[i] = Mi / (rho * areaU[i]);
        }
    }

    if (DEBUG) {
        printMatrix("xP", xP);
        printMatrix("areaP", areaP);
        printMatrix("areaU", areaU);
        printMatrix("u", u);
        printMatrix("p", p);
    }

    printf(
        "Iterations\tU* iters\tP' iters\tL2 norm\t\tP error\t\tU error\n");

    //  Outer iterations
    start = omp_get_wtime();
    for (int outer = 0; outer < MAXITER; outer++) {
        //  Save previous iteration values
        pPrev = p;
        uPrev = u;

        //  Solve momentum equation for u*
        Mcount = 0;
        uStar = u;
        //  build momentum coefficients
        buildMomentumCoeffs(p, areaP, uPrev, areaU, Su, aPu, aWu, aEu, aPuinv,
                            d, rho);

        //  Iterate until uStar converges
        dif = 1.0;
        while (dif > THRESHinner) {
            dif = 0.0;
#pragma omp parallel for private(UP, UE, UW, temp, rhs) reduction(+ : dif)
            for (int i = 0; i < Nm1; i++) {
                UP = uStar[i];
                switch (i) {
                    //  Edge cases
                    case 0: {
                        rhs = Su[i];
                    } break;
                    case Nm2: {
                        UW = uStar[i - 1];
                        rhs = aWu[i] * UW + Su[i];
                    } break;

                    //  Interior
                    default: {
                        UE = uStar[i + 1];
                        UW = uStar[i - 1];
                        rhs = aWu[i] * UW + aEu[i] * UE + Su[i];
                    }
                }
                //  calculate residual
                temp = UP;
                UP = rhs * aPuinv[i];
                //  Under-relaxation
                UP = OmOMEGAu * temp + OMEGAu * UP;
                uStar[i] = UP;
                dif += P2((UP - temp) / (UP + temp));///P2(rhs);
            }
            dif = sqrt(dif / Nm1);
            // std::cout<<dif<<" "<<THRESHinner<<std::endl;
            Mcount++;

            if (DEBUG) std::cout << dif << std::endl;
        }
        //std::cout<<"sad\n";

        if (DEBUG) printMatrix("uStar", uStar);

        //  Iteratively solve for pressure correction
        pPrime[0] = 0.0;
        pPrime[Nm1] = 0.0;
        dif = 1.0;
        Pprimecount = 0;

        //  build vector of coefficients
        buildPressureCoeffs(aWp, aEp, aPp, aPpinv, bPrimep, d, uStar, areaU, rho);


        //  Iterate until P' converges
        // int hi;
        while (dif > THRESHinner) {
            dif = 0.0;
            #pragma omp parallel for private(temp, PP, rhs) reduction(+ : dif)
            // hi = 0;
            for (int i = 1; i < Nm1; i++) {
                temp = pPrime[i];
                // std::cout << "AAADYOTTTTT" << temp << std::endl;
                rhs = aWp[i] * pPrime[i - 1] + aEp[i] * pPrime[i + 1] +
                      bPrimep[i];
                //  calculate residual
                PP = rhs * aPpinv[i];
                //  under-relaxation
                PP = OmOMEGAp * temp + OMEGAp * PP;
                pPrime[i] = PP;
                dif +=P2((PP - temp) / (PP + temp));
            }
            // std::cout << "PATSALLLL" << hi << std::endl;
            // if (hi == Nm2) break;
            dif = sqrt(dif / Nm2);
            Pprimecount++;
            if (DEBUG) std::cout << dif << std::endl;
        }

        if (DEBUG) printMatrix("pPrime", pPrime);

        //  Correct pressure and velocity
        totaldif = 0.0;
#pragma omp parallel for private(PP, UP, temp) reduction(+ : totaldif)
        for (int i = 0; i < Nm1; i++) {
            //  Correct all velocities
            UP = u[i];
            temp = UP;
            UP = uStar[i] + d[i] * (pPrime[i] - pPrime[i + 1]);
            UP = OMEGAu * UP + OmOMEGAu * temp;
            u[i] = UP;
            // totaldif += P2(uPrev[i] - UP);///P2(UP);
            totaldif += P2((uPrev[i] - UP) / (uPrev[i] + UP));///P2(UP);

            //  Correct pressure, except on edges
            PP = p[i];
            if (i == 0) {
                PP = Pi - 0.5 * P2(u[0]) * P2(areaU[0] / areaP[0]);
            } else {
                //  under-relaxation
                PP += OMEGAp * pPrime[i];
            }
            p[i] = PP;
            //std::cout<<"oho"<<totaldif<<" "<<p[i]<<std::endl;
            // totaldif += P2(pPrev[i] - p[i]);//P2(p[i]);
            totaldif += P2((pPrev[i] - p[i]) / (pPrev[i] + p[i]));//P2(p[i]);
        }

        //  Calculate L2 norm of pressure and velocity difference
        if (outer == 0) totaldif0 = 1.0 / sqrt(totaldif / (2.0 * Nm1));
        //std::cout<<"aaaaaa"<<totaldif<<" "<<totaldif0<<std::endl;
        totaldif = sqrt(totaldif / ((2.0 * Nm1))) * totaldif0;

        //  Print information
        if (outer % printint == 0) {
            //  Calculate error from exact solution
            solnError(pError, uError, p, Pexsoln, u, Uexsoln);
            printf("%.2e\t%.2e\t%.2e\t%.2e\t%.3e\t%.3e\t\n",
                       static_cast<double>(outer), static_cast<double>(Mcount),
                       static_cast<double>(Pprimecount), totaldif, pError,
                       uError);

            //  If converged, print final information
            if (totaldif < THRESH) {
                //  Calculate mass flow rate
#pragma omp parallel for
                for (int i = 0; i < N; i++) mfr[i] = rho * areaU[i] * u[i];
                #ifdef SR
                for (int i = 0; i < N; i++) u[i] = u[i] / pow(1 + pow(u[i] / speedoflight, 2), 0.5);
                #endif

                stop = omp_get_wtime();
                printMatrix("Mass Flow Rate", mfr);
                printMatrix("u", u);
                printMatrix("p", p);
                printf("\nTime: %.3e s", stop - start);
                break;
            }
        }
        if (isnan(totaldif) || isinf(totaldif) || !isnormal(dif)){
            std::cout<<"cry"<<isnan(totaldif) << isinf(totaldif) << !isnormal(dif) << std::endl;
            throwError("dif error\n");
        }
    }
    return 0;
}

//  power of 2
double P2(const double &value) { return abs(value) ;}//* value; }

//  Prints a matrix
template <std::size_t SIZE>
void printMatrix(const std::string &name, const std::array<double, SIZE> &vec) {
    printf("\n%s:\n", name.c_str());
    for (size_t i = 0; i < SIZE; i++) printf("%.4f\t", vec[i]);
    printf("\n");
}

//  Calculates solution error
void solnError(double &perror, double &uerror, const Pvec &p,
               const Pvec &pexact, const Uvec &u, const Uvec &uexact) {
    perror = 0.0;
    uerror = 0.0;
#pragma omp parallel for reduction(+ : perror, uerror)
    for (int i = 0; i < N; i++) {
        perror += P2(p[i] - pexact[i]);
        if (i < Nm1) uerror += P2(u[i] - uexact[i]);
    }
    uerror = sqrt(uerror / Nm1);
    perror = sqrt(perror / N);
}

//  Exact solution for pressure
double Pexact(const double &area) { return 10.0 - 0.5 * P2(Mexact) / P2(area); }

//  Exact solution for velocity
double Uexact(const double &area) { return Mexact / area; }


//  Build Momentum coefficient vectors
void buildMomentumCoeffs(const Pvec &p, const Pvec &areaP, const Uvec &uPrev,
                         const Uvec &areaU, Uvec &Su, Uvec &aPu, Uvec &aWu,
                         Uvec &aEu, Uvec &aPuinv, Uvec &d, const double rho) {
    double UP{}, UE{}, UW{}, Aw{}, Ae{}, Fw{}, Fe{}, Pe{}, Pw{}, Vol{}, dPdx{};


#pragma omp parallel for private(UP, UE, UW, Aw, Ae, Fw, Fe, Pe, Pw, Vol, dPdx)
    for (int i = 0; i < Nm1; i++) {

        #ifdef SR

        double P_face = 0.5 * (p[i] + p[i + 1]);  // interior faces
        // std::cout << p[0] << "AAAAAA" << std::endl;
        // if (i == 0){
        //     P_face = Pi;
        // }
        double rho_new = (rho + enthalpy_coeff * P_face / pow(speedoflight, 2));
        #else
        double rho_new = rho;
        #endif

        // std::cout << rho_new << "AAAAAA" << std::endl;
        // exit(0);


        //  Cell-face pressures
        if (i == 0)
            Pw = Pi;
        else
            Pw = p[i];
        Pe = p[i + 1];

        //  Cell Volume
        Vol = areaU[i] * dx;

        //  dp/dx
        dPdx = (Pw - Pe) / dx;

        //  Source term
        UP = uPrev[i];
        #ifdef SR
        Su[i] = dPdx * (1 + enthalpy_coeff * UP * UP / pow(speedoflight, 2)) * Vol;
        #else
        Su[i] = dPdx * Vol;
        #endif

        // if (i == 0) Su[i] += P2(UP * areaU[i]) / areaP[i];
        if (i==0) Su[0] += P2(UP * areaU[i]) / areaP[i];


        //  Neighbor velocities
        if (i != Nm2) UE = uPrev[i + 1];
        if (i != 0) UW = uPrev[i - 1];

        //  Cell-face areas
        Aw = areaP[i];
        Ae = areaP[i + 1];

        //  Flux terms (old code assumed rho == 1)
        if (i == 0){
            // Fw = areaU[i] * UP;
            Fw = rho_new * areaU[i] * UP;
        }else{
            // Fw = 0.5 * Aw * (UP + UW);
            Fw = rho_new * 0.5 * Aw * (UP + UW);
        }

        if (i == Nm2){
            // Fe = areaU[i] * UP;
            Fe = rho_new * areaU[i] * UP;
        }else{
            // Fe = 0.5 * Ae * (UP + UE);
            Fe = rho_new * 0.5 * Ae * (UP + UE);
        }
        //  Upwind difference coefficients
        if (i == 0)
            aWu[i] = 0.0;
        else
            aWu[i] = Dw + std::max(Fw, 0.0);

        if (i == Nm2)
            aEu[i] = 0.0;
        else
            aEu[i] = De + std::max(0.0, -Fe);

        if (i == 0)
            aPu[i] = Fe + Fw * 0.5 * P2(areaU[i] / areaP[i]);
        else
            aPu[i] = aWu[i] + aEu[i] + (Fe - Fw);

        aPuinv[i] = 1.0 / aPu[i];

        d[i] = areaU[i] * aPuinv[i];
    }
}

//  build vector of pressure coefficients
void buildPressureCoeffs(Pvec &aWp, Pvec &aEp, Pvec &aPp, Pvec &aPpinv,
                         Pvec &bPrimep, const Uvec &d, const Uvec &usource,
                         const Uvec &areaU, const double rho) {
    double Fw{}, Fe{};
#pragma omp parallel for private(Fw, Fe)
    for (int i = 1; i < Nm1; i++) {  //old code assumed rho = 1
        aWp[i] = d[i - 1] * areaU[i - 1] * rho;
        aEp[i] = d[i] * areaU[i] * rho;
        Fw = usource[i - 1] * areaU[i - 1] * rho;
        Fe = usource[i] * areaU[i] * rho;

        aPp[i] = aWp[i] + aEp[i];
        aPpinv[i] = 1.0 / aPp[i];
        bPrimep[i] = Fw - Fe;
    }
}

//  Check for valid parameter values
void checkparams() {
    if (N < 3) throwError("Error: invalid grid size.  Pick N > 2\n");

    if (THRESH < 0 || THRESH > 1)
        throwError("Error: invalid convergence threshold\n");

    if (THRESH > 1E-5)
        printf("Warning: threshold is recommended to be at or below 1E-5\n");

    if (OMEGAu > 1 || OMEGAu < 0 || OMEGAp > 1 || OMEGAp < 0)
        throwError("Error: relaxation should be between 0 and 1\n");

    if (DEBUG != 0 && DEBUG != 1)
        throwError("Error: invalid debug value.  Pick 0 or 1\n");

    if (printint < 1)
        throwError("Invalid print interval.  Must be greater than 0.");

}

//  Throw error message and exit
void throwError(const std::string &message) {
    std::cout << message << std::endl;
    exit(1);
}
