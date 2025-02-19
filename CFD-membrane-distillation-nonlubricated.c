#include "udf.h"

#define r 1.1e-7      // Average pore diameter of the membrane
#define tor 1.12      // Membrane tortuosity factor
#define porosity 0.75 // Membrane porosity
#define height 1e-5   // First layer grid height
#define kc 0.0546     // Membrane thermal conductivity
#define b 1.2e-4      // Membrane thickness

DEFINE_INIT(init_udm, domain)
{
    Thread *thread_membrane = Lookup_Thread(domain, 19); // Membrane boundary
    Thread *thread_membraneshadow = Lookup_Thread(domain, 2); // Membraneshadow boundary
    Thread *thread_feed = Lookup_Thread(domain, 13); // Feed region
    Thread *thread_permeat = Lookup_Thread(domain, 12); // Permeate region

    face_t face;
    cell_t cell;

    begin_c_loop(cell, thread_feed)
    {
        C_UDMI(cell, thread_feed, 0) = -1;
    }
    end_c_loop(cell, thread_feed);

    begin_c_loop(cell, thread_permeat)
    {
        C_UDMI(cell, thread_permeat, 0) = -1;
    }
    end_c_loop(cell, thread_permeat);

    begin_f_loop(face, thread_membrane)
    {
        C_UDMI(F_C0(face, thread_membrane), THREAD_T0(thread_membrane), 0) =
            F_C0(face, thread_membraneshadow);
        C_UDMI(F_C0(face, thread_membraneshadow), THREAD_T0(thread_membraneshadow), 0) =
            F_C0(face, thread_membrane);
    }
    end_f_loop(face, thread_membrane);
}

DEFINE_SOURCE(feed_water_source, cell_feed, thread_feed, dS, eqn)
{
    real source, temp0, temp1, xm, xa, rw, p0, pa, pc, mean_t, k1, k2, km;
    Domain *domain = Get_Domain(1);
    Thread *thread_permeat = Lookup_Thread(domain, 12);
    cell_t cell_permeat;

    if (C_UDMI(cell_feed, thread_feed, 0) == -1)
    {
        source = 0.0;
    }
    else
    {
        cell_permeat = C_UDMI(cell_feed, thread_feed, 0);
        temp0 = C_T(cell_permeat, thread_permeat); // Temperature on the permeate (cold) side
        temp1 = C_T(cell_feed, thread_feed);         // Temperature on the feed (hot) side
        mean_t = (temp0 + temp1) / 2;
        xm = C_YI(cell_feed, thread_feed, 1);          // Mass fraction of NaCl (hot side)
        xa = 18 * xm / (58.5 - 40.5 * xm);              // Mole fraction of NaCl (hot side)
        rw = 1 - 0.5 * xa - 10 * xa * xa;               // Activity coefficient of water
        p0 = exp(23.238 - 3841 / (temp1 - 45));         // Saturation vapor pressure (hot side)
        pa = p0 * (1 - xa) * rw;                        // Partial vapor pressure of water (hot side)
        pc = exp(23.238 - 3841 / (temp0 - 45));         // Saturation vapor pressure (cold side)
        k1 = 1.064 * r * porosity / tor / b * pow(0.018 / 8.314 / mean_t, 0.5); // Knudsen diffusion coefficient
        k2 = porosity / tor / b * (0.018 / 8.315 / mean_t) * 1.895e-5 * pow(mean_t, 2.072) / 1.01e5; // Molecular diffusion coefficient
        km = 1 / (1 / k1 + 1 / k2);                     // Overall mass transfer coefficient
        source = -km * (pa - pc) / height;              // Mass transfer source term (kg/m3/s)
    }
    dS[eqn] = 0.0;
    return source;
}

DEFINE_SOURCE(permeat_water_source, cell_permeat, thread_permeat, dS, eqn)
{
    real source, temp0, temp1, xm, xa, rw, p0, pa, pc, mean_t, k1, k2, km;
    Domain *domain = Get_Domain(1);
    Thread *thread_feed = Lookup_Thread(domain, 13);
    cell_t cell_feed;

    if (C_UDMI(cell_permeat, thread_permeat, 0) == -1)
    {
        source = 0.0;
    }
    else
    {
        cell_feed = C_UDMI(cell_permeat, thread_permeat, 0);
        temp0 = C_T(cell_permeat, thread_permeat); // Temperature on the permeate side
        temp1 = C_T(cell_feed, thread_feed);         // Temperature on the feed side
        mean_t = (temp0 + temp1) / 2;
        xm = C_YI(cell_feed, thread_feed, 1);          // Mass fraction of NaCl (hot side)
        xa = 18 * xm / (58.5 - 40.5 * xm);              // Mole fraction of NaCl (hot side)
        rw = 1 - 0.5 * xa - 10 * xa * xa;               // Activity coefficient of water
        p0 = exp(23.238 - 3841 / (temp1 - 45));         // Saturation vapor pressure (hot side)
        pa = p0 * (1 - xa) * rw;                        // Partial vapor pressure of water (hot side)
        pc = exp(23.238 - 3841 / (temp0 - 45));         // Saturation vapor pressure (cold side)
        k1 = 1.064 * r * porosity / tor / b * pow(0.018 / 8.314 / mean_t, 0.5); // Knudsen diffusion coefficient
        k2 = porosity / tor / b * (0.018 / 8.315 / mean_t) * 1.895e-5 * pow(mean_t, 2.072) / 1.01e5; // Molecular diffusion coefficient
        km = 1 / (1 / k1 + 1 / k2);                     // Overall mass transfer coefficient
        source = km * (pa - pc) / height;              // Mass transfer source term (kg/m3/s)
    }
    dS[eqn] = 0.0;
    return source;
}

DEFINE_PROFILE(memshadow_heat_flux, thread, i)
{
    face_t f;
    Thread *t0;
    Thread *t1;
    real tem0, tem1, qr, pa, xm, xa, rw, p0, pc, Qc, Jv, mean_t, k1, k2, km;
    begin_f_loop(f, thread)
    {
        cell_t c0, c1;
        c0 = F_C0(f, thread); // Hot side membrane cell
        t0 = THREAD_T0(thread);
        c1 = F_C1(f, thread); // Cold side membrane cell
        t1 = THREAD_T1(thread);
        tem0 = C_T(c0, t0); // Temperature at the hot side
        tem1 = C_T(c1, t1); // Temperature at the cold side
        mean_t = (tem0 + tem1) / 2;
        xm = C_YI(c0, t0, 1); // Mass fraction of NaCl (hot side)
        xa = 18 * xm / (58.5 - 40.5 * xm); // Mole fraction of NaCl (hot side)
        rw = 1 - 0.5 * xa - 10 * xa * xa;    // Activity coefficient of water
        p0 = exp(23.238 - 3841 / (tem0 - 45)); // Saturation vapor pressure (hot side)
        pa = p0 * (1 - xa) * rw;              // Partial vapor pressure of water (hot side)
        pc = exp(23.238 - 3841 / (tem1 - 45)); // Saturation vapor pressure (cold side)
        k1 = 1.064 * r * porosity / tor / b * pow(0.018 / 8.314 / mean_t, 0.5);
        k2 = porosity / tor / b * (0.018 / 8.315 / mean_t) * 1.895e-5 * pow(mean_t, 2.072) / 1.01e5;
        km = 1 / (1 / k1 + 1 / k2);
        Jv = km * (pa - pc);                // Mass transfer flux
        qr = (-0.001351 * tem0 * tem0 - 1.4461 * tem0 + 2986.5) * 1e+3; // Latent heat of vaporization
        Qc = (kc / b) * (tem0 - tem1);         // Membrane heat conduction
        F_PROFILE(f, thread, i) = -(Jv * qr + Qc); // Total heat flux
    }
    end_f_loop(f, thread)
}

DEFINE_PROFILE(mem_heat_flux, thread, i)
{
    face_t f;
    Thread *t0;
    Thread *t1;
    real tem0, tem1, qr, pa, xa, rw, p0, pc, Qc, Jv, xm, mean_t, k1, k2, km;
    begin_f_loop(f, thread)
    {
        cell_t c0, c1;
        c0 = F_C0(f, thread); // Hot side membrane cell
        t0 = THREAD_T0(thread);
        c1 = F_C1(f, thread); // Cold side membrane cell
        t1 = THREAD_T1(thread);
        tem0 = C_T(c0, t0); // Temperature at the cold side
        tem1 = C_T(c1, t1); // Temperature at the hot side
        mean_t = (tem0 + tem1) / 2;
        xm = C_YI(c1, t1, 1); // Mass fraction of NaCl (hot side)
        xa = 18 * xm / (58.5 - 40.5 * xm); // Mole fraction of NaCl (hot side)
        rw = 1 - 0.5 * xa - 10 * xa * xa;    // Activity coefficient of water
        p0 = exp(23.238 - 3841 / (tem1 - 45)); // Saturation vapor pressure (hot side)
        pa = p0 * (1 - xa) * rw;              // Partial vapor pressure of water (hot side)
        pc = exp(23.238 - 3841 / (tem0 - 45)); // Saturation vapor pressure (cold side)
        k1 = 1.064 * r * porosity / tor / b * pow(0.018 / 8.314 / mean_t, 0.5);
        k2 = porosity / tor / b * (0.018 / 8.315 / mean_t) * 1.895e-5 * pow(mean_t, 2.072) / 1.01e5;
        km = 1 / (1 / k1 + 1 / k2);
        Jv = km * (pa - pc);                // Mass transfer flux
        qr = (-0.001351 * tem1 * tem1 - 1.4461 * tem1 + 2986.5) * 1e+3; // Latent heat of vaporization
        Qc = (kc / b) * (tem1 - tem0);         // Membrane heat conduction
        F_PROFILE(f, thread, i) = (Jv * qr + Qc);  // Total heat flux
    }
    end_f_loop(f, thread)
}
