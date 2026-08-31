/* ================================================================
   ternary_fdtd_fixed_scattering.c

   Ternary (-1,0,+1) FDTD compared against a floating-point Yee
   reference on the same 80x80 PEC cavity.

   Key differences from the original ternary scheme:

     1. The state accumulator is a TRUE integrator.  It is never
        drained by the quantiser.  (In the original code the
        quantiser subtracted THRESHOLD from the accumulator, which
        destroyed the time integration of the Yee update.)

     2. Quantisation is done by a separate first-order noise-shaped
        ENCODER with its own error register.  Sub-threshold cells
        emit 0, not the previous symbol.

     3. CFL appears explicitly, as an exact integer ratio in
        fixed point (MC = M * CFL).

   Fixed-point convention:

     symbol s in {-1,0,+1}   represents the physical field  s * AMP
     state  a (int32)        represents the physical field  a * AMP / M

   ================================================================ */

#include <stdio.h>
#include <stdint.h>
#include <math.h>

#define NX      80
#define NY      80
#define STEPS   8192

#define CFL     0.1f

/* fixed-point units per ternary symbol */
#define M       100000
/* M * CFL, exact integer */
#define MC      10000

/* physical field value represented by symbol = 1 (full scale) */
#define AMP     20.0

#define SOURCE_AMPLITUDE  4.0
#define SOURCE_FREQUENCY  0.01

/* ----------------------------------------------------------------
   PEC scatterer in the top-right quadrant.

   SHAPE 0 = none, 1 = square, 2 = circle.

   The wavelength here is CFL/SOURCE_FREQUENCY = 10 cells, so the
   grid cannot resolve curvature: a diameter-9 circle digitises to
   an octagon whose boundary deviates from the true circle by up to
   half a cell (lambda/20), and whose area is 8.5 % too large.  A
   grid-aligned square is exactly representable, so the geometry
   contributes no error at all and any ternary/FP discrepancy is
   attributable to the encoder rather than to two different
   staircase realisations.  Side 9 = 0.9 lambda gives strong
   scattering.  Use SHAPE 2 with SCAT_HALF 6.5 if curvature is
   wanted; the area error is then 3.2 %.
   ---------------------------------------------------------------- */
#define SCAT_SHAPE   2
#define SCAT_CI      60.0        /* centre row    */
#define SCAT_CJ      60.0        /* centre column */
#define SCAT_HALF    4.0         /* half-side, or radius if circle */

/* ---------------- ternary symbols ---------------- */
static int8_t  Ez[NX][NY], Hx[NX][NY], Hy[NX][NY];

/* ---------------- state integrators (never drained) -------- */
static int32_t aEz[NX][NY], aHx[NX][NY], aHy[NX][NY];

/* ---------------- encoder error registers ------------------ */
static int32_t eEz[NX][NY], eHx[NX][NY], eHy[NX][NY];
static int32_t fEz[NX][NY], fHx[NX][NY], fHy[NX][NY];

/* ---------------- floating point reference ----------------- */
static float fp_Ez[NX][NY], fp_Hx[NX][NY], fp_Hy[NX][NY];

/* ---------------- PEC scatterer mask ----------------------- */
static int8_t pec[NX][NY];
static int    pec_cells = 0;

/* saturation counter (diagnostic) */
static long sat_count = 0;

/* ----------------------------------------------------------------
   First-order noise-shaped ternary encoder.

     v = x + err          (x = state, err = previous rounding error)
     s = clamp(round(v/M), -1, +1)
     err = v - s*M

   Local mean of s equals x/M ; the rounding error is high-pass
   shaped, i.e. pushed towards the Nyquist frequency.
   ---------------------------------------------------------------- */
static inline int8_t encode(int32_t x, int32_t *e1, int32_t *e2)
{
    int32_t v = x + (3*(*e1) - (*e2))/2;   /* NTF = (1-z^-1)(1-0.5z^-1) */
    int32_t s;

    if (v >= 0) s =  (v + M/2) / M;
    else        s = -((-v + M/2) / M);

    if (s >  1) { s =  1; sat_count++; }
    if (s < -1) { s = -1; sat_count++; }

    *e2 = *e1;
    *e1 = v - s * M;

    if (*e1 >  2*M) *e1 =  2*M;      /* anti-windup */
    if (*e1 < -2*M) *e1 = -2*M;

    return (int8_t)s;
}


/* ======== Dump Ez field (ternary) ======== */
void dump_snapshot(int step)
{
    int i,j;
    char filename[64];
    sprintf(filename,"snapshot_%04d.txt",step);
    FILE *fp=fopen(filename,"w");
    if(fp==NULL){ perror("file"); return; }
    fprintf(fp,"STEP %d\n\n",step);
    for(i=0;i<NX;i++){
        for(j=0;j<NY;j++) fprintf(fp,"%6d ",Ez[i][j]);
        fprintf(fp,"\n");
    }
    fclose(fp);
}

/* ======== Dump Ez field (floating point) ======== */
void dump_fp_snapshot(int step)
{
    int i,j;
    char filename[64];
    sprintf(filename,"fp_snapshot_%04d.txt",step);
    FILE *fp=fopen(filename,"w");
    if(fp==NULL){ perror("file"); return; }
    fprintf(fp,"STEP %d\n\n",step);
    for(i=0;i<NX;i++){
        for(j=0;j<NY;j++) fprintf(fp,"%12.6f ",fp_Ez[i][j]);
        fprintf(fp,"\n");
    }
    fclose(fp);
}

/* ----------------------------------------------------------------
   Build the scatterer mask.
   ---------------------------------------------------------------- */
static void build_scatterer(void)
{
    int i, j;

    pec_cells = 0;

    for (i = 0; i < NX; i++)
        for (j = 0; j < NY; j++)
        {
            double di = i - SCAT_CI;
            double dj = j - SCAT_CJ;
            int in = 0;

#if   SCAT_SHAPE == 1
            in = (fabs(di) <= SCAT_HALF && fabs(dj) <= SCAT_HALF);
#elif SCAT_SHAPE == 2
            in = (di*di + dj*dj <= SCAT_HALF*SCAT_HALF);
#endif
            pec[i][j] = (int8_t)in;
            pec_cells += in;
        }
}

/* ----------------------------------------------------------------
   Enforce Ez = 0 inside the PEC.

   For the ternary field it is not enough to zero the emitted
   symbol: the state accumulator and both encoder error registers
   must be cleared too.  Otherwise the accumulator keeps integrating
   the curl inside the conductor and the error registers wind up
   against their anti-windup limits, which is not the discrete
   analogue of fp_Ez = 0.
   ---------------------------------------------------------------- */
static void apply_pec(void)
{
    int i, j;

    for (i = 0; i < NX; i++)
        for (j = 0; j < NY; j++)
            if (pec[i][j])
            {
                Ez[i][j]  = 0;
                aEz[i][j] = 0;
                eEz[i][j] = 0;
                fEz[i][j] = 0;

                fp_Ez[i][j] = 0.0f;
            }
}

int main(void)
{
    int n, i, j;

    const int sx = NX/2, sy = NY/2;
    const int px = 32,   py = 32;

    double srem = 0.0;              /* source rounding residual */

    FILE *probe = fopen("probe_fixed.txt", "w");
    if (!probe) { perror("probe_fixed.txt"); return 1; }

    fprintf(probe, "# step ternary_symbol fp_Ez state_Ez\n");

    for (i = 0; i < NX; i++)
        for (j = 0; j < NY; j++) {
            Ez[i][j]=Hx[i][j]=Hy[i][j]=0;
            aEz[i][j]=aHx[i][j]=aHy[i][j]=0;
            eEz[i][j]=eHx[i][j]=eHy[i][j]=0;
            fEz[i][j]=fHx[i][j]=fHy[i][j]=0;
            fp_Ez[i][j]=fp_Hx[i][j]=fp_Hy[i][j]=0.0f;
        }

    build_scatterer();

    fprintf(stderr, "PEC scatterer: shape %d, %d cells\n",
            SCAT_SHAPE, pec_cells);

    for (n = 0; n < STEPS; n++)
    {
        /* ======== ternary: H update ======== */
        for (i = 0; i < NX-1; i++)
            for (j = 0; j < NY-1; j++)
            {
                aHx[i][j] += MC * (Ez[i][j] - Ez[i][j+1]);
                aHy[i][j] += MC * (Ez[i+1][j] - Ez[i][j]);

                Hx[i][j] = encode(aHx[i][j], &eHx[i][j], &fHx[i][j]);
                Hy[i][j] = encode(aHy[i][j], &eHy[i][j], &fHy[i][j]);
            }

        /* ======== ternary: E update ======== */
        for (i = 1; i < NX-1; i++)
            for (j = 1; j < NY-1; j++)
            {
                aEz[i][j] += MC * ((Hy[i][j] - Hy[i-1][j])
                                 - (Hx[i][j] - Hx[i][j-1]));

                Ez[i][j] = encode(aEz[i][j], &eEz[i][j], &fEz[i][j]);
            }

        /* ======== ternary: source ======== */
        {
            double v = (M/AMP) * SOURCE_AMPLITUDE *
                       sin(2.0*M_PI*SOURCE_FREQUENCY*n) + srem;
            int32_t inc = (int32_t)lrint(v);
            srem = v - (double)inc;

            aEz[sx][sy] += inc;
            Ez[sx][sy]   = encode(aEz[sx][sy], &eEz[sx][sy], &fEz[sx][sy]);
        }

        /* ======== ternary: PEC boundaries ======== */
        for (i = 0; i < NX; i++) { Ez[i][0]=0; Ez[i][NY-1]=0; }
        for (j = 0; j < NY; j++) { Ez[0][j]=0; Ez[NX-1][j]=0; }

        /* ======== ternary + FP: PEC scatterer ======== */
        apply_pec();

        /* ======== floating point reference ======== */
        for (i = 0; i < NX-1; i++)
            for (j = 0; j < NY-1; j++)
            {
                fp_Hx[i][j] += CFL * (fp_Ez[i][j] - fp_Ez[i][j+1]);
                fp_Hy[i][j] += CFL * (fp_Ez[i+1][j] - fp_Ez[i][j]);
            }

        for (i = 1; i < NX-1; i++)
            for (j = 1; j < NY-1; j++)
                fp_Ez[i][j] += CFL * ((fp_Hy[i][j] - fp_Hy[i-1][j])
                                    - (fp_Hx[i][j] - fp_Hx[i][j-1]));

        fp_Ez[sx][sy] += (float)(SOURCE_AMPLITUDE *
                          sin(2.0*M_PI*SOURCE_FREQUENCY*n));

        for (i = 0; i < NX; i++) { fp_Ez[i][0]=0.0f; fp_Ez[i][NY-1]=0.0f; }
        for (j = 0; j < NY; j++) { fp_Ez[0][j]=0.0f; fp_Ez[NX-1][j]=0.0f; }

        apply_pec();

        fprintf(probe, "%d %d %.9f %.9f\n", n, (int)Ez[px][py], fp_Ez[px][py], aEz[px][py]*AMP/(double)M);

        /* ======== Snapshots ======== */
        if ((n % 20) == 0)
        {
            dump_snapshot(n);
            dump_fp_snapshot(n);
        }
    }

    fclose(probe);

    fprintf(stderr, "encoder saturation events: %ld  (%.4f %% of encodes)\n",
            sat_count, 100.0*sat_count/((double)STEPS*3.0*NX*NY));

    return 0;
}
