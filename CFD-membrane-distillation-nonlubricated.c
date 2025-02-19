#include "udf.h" 

#define r 1.1e-7  /* Membrane average pore size */ 
#define tor 1.12  /* Membrane tortuosity factor */ 
#define porosity 0.75 /* Membrane porosity */ 
#define height 1e-5  /* First layer grid height */ 
#define kc 0.0546  /* Membrane thermal conductivity */ 
#define b 1.2e-4  /* Membrane thickness */ 

DEFINE_INIT(init_udm, domain) /* Initialize UDM values and assign cell IDs */ 
{ 
    Thread *thread_membrane = Lookup_Thread(domain, 19); /* Membrane boundary thread */ 
    Thread *thread_membraneshadow = Lookup_Thread(domain, 2); /* Shadow membrane boundary thread */ 
    Thread *thread_feed = Lookup_Thread(domain, 13); /* Feed region thread */ 
    Thread *thread_permeat = Lookup_Thread(domain, 12); /* Permeate region thread */ 
    face_t face; 
    cell_t cell; 

    begin_c_loop(cell, thread_feed) /* Loop through feed region and set initial UDM values */ 
    { 
        C_UDMI(cell, thread_feed, 0) = -1; 
    } 
    end_c_loop(cell, thread_feed); 

    begin_c_loop(cell, thread_permeat) /* Loop through permeate region and set initial UDM values */ 
    { 
        C_UDMI(cell, thread_permeat, 0) = -1; 
    } 
    end_c_loop(cell, thread_permeat); 

    begin_f_loop(face, thread_membrane) /* Loop through membrane boundary and mark both sides */ 
    { 
        C_UDMI(F_C0(face, thread_membrane), THREAD_T0(thread_membrane), 0) = 
            F_C0(face, thread_membraneshadow); 

        C_UDMI(F_C0(face, thread_membraneshadow), 
            THREAD_T0(thread_membraneshadow), 0) = F_C0(face, thread_membrane); 
    } 
    end_f_loop(face, thread_membrane); 
} 

DEFINE_SOURCE(feed_water_source, cell_feed, thread_feed, dS, eqn) /* Mass source term for feed side */ 
{ 
    real source, temp0, temp1, xm, xa, rw, p0, pa, pc, mean_t, k1, k2, km; 
    Domain *domain = Get_Domain(1); /* Domain pointer */ 
    Thread *thread_permeat = Lookup_Thread(domain, 12); /* Permeate region pointer */ 
    cell_t cell_permeat; 

    if (C_UDMI(cell_feed, thread_feed, 0) == -1) /* Check if not a membrane cell */ 
    { 
        source = 0.; /* Set source term to zero */ 
    } 
    else /* If it is a membrane cell */ 
    { 
        cell_permeat = C_UDMI(cell_feed, thread_feed, 0); /* Get corresponding cell ID on permeate side */ 

        /* Compute mass transfer source term */ 
        temp0 = C_T(cell_permeat, thread_permeat); /* Permeate side temperature */ 
        temp1 = C_T(cell_feed, thread_feed); /* Feed side temperature */ 
        mean_t = (temp0 + temp1) / 2; 

        xm = C_YI(cell_feed, thread_feed, 1); /* NaCl mass fraction at membrane surface (feed side) */ 
        xa = 18 * xm / (58.5 - 40.5 * xm); /* Mole fraction of NaCl at membrane surface */ 
        rw = 1 - 0.5 * xa - 10 * xa * xa; /* Water activity coefficient */ 

        p0 = exp(23.238 - 3841 / (temp1 - 45)); /* Saturation vapor pressure on feed side */ 
        pa = p0 * (1 - xa) * rw; /* Water vapor partial pressure on feed side */ 
        pc = exp(23.238 - 3841 / (temp0 - 45)); /* Saturation vapor pressure on permeate side */ 

        k1 = 1.064 * r * porosity / tor / b * pow(0.018 / 8.314 / mean_t, 0.5); /* Knudsen diffusion coefficient */ 
        k2 = porosity / tor / b * (0.018 / 8.315 / mean_t) * 1.895e-5 * pow(mean_t, 2.072) / 1.01e5; /* Molecular diffusion coefficient */ 

        km = 1 / (1 / k1 + 1 / k2); /* Total mass transfer coefficient */ 
        source = -km * (pa - pc) / height; /* Mass source term, unit: kg/m³/s */ 
    } 

    dS[eqn] = 0.0; 
    return source; 
} 

DEFINE_SOURCE(permeat_water_source, cell_permeat, thread_permeat, dS, eqn) /* Mass source term for permeate side */ 
{ 
    real source, temp0, temp1, xm, xa, rw, p0, pa, pc, mean_t, k1, k2, km; 
    Domain *domain = Get_Domain(1); /* Domain pointer */ 
    Thread *thread_feed = Lookup_Thread(domain, 13); /* Feed region pointer */ 
    cell_t cell_feed; 

    if (C_UDMI(cell_permeat, thread_permeat, 0) == -1) /* Check if not a membrane cell */ 
    { 
        source = 0.; /* Set source term to zero */ 
    } 
    else /* If it is a membrane cell */ 
    { 
        cell_feed = C_UDMI(cell_permeat, thread_permeat, 0); /* Get corresponding feed side cell */ 

        /* Compute mass transfer source term */ 
        temp0 = C_T(cell_permeat, thread_permeat); /* Permeate side temperature */ 
        temp1 = C_T(cell_feed, thread_feed); /* Feed side temperature */ 
        mean_t = (temp0 + temp1) / 2; 

        xm = C_YI(cell_feed, thread_feed, 1); /* NaCl mass fraction at membrane surface (feed side) */ 
        xa = 18 * xm / (58.5 - 40.5 * xm); /* Mole fraction of NaCl at membrane surface */ 
        rw = 1 - 0.5 * xa - 10 * xa * xa; /* Water activity coefficient */ 

        p0 = exp(23.238 - 3841 / (temp1 - 45)); /* Saturation vapor pressure on feed side */ 
        pa = p0 * (1 - xa) * rw; /* Water vapor partial pressure on feed side */ 
        pc = exp(23.238 - 3841 / (temp0 - 45)); /* Saturation vapor pressure on permeate side */ 

        k1 = 1.064 * r * porosity / tor / b * pow(0.018 / 8.314 / mean_t, 0.5); /* Knudsen diffusion coefficient */ 
        k2 = porosity / tor / b * (0.018 / 8.315 / mean_t) * 1.895e-5 * pow(mean_t, 2.072) / 1.01e5; /* Molecular diffusion coefficient */ 

        km = 1 / (1 / k1 + 1 / k2); /* Total mass transfer coefficient */ 
        source = km * (pa - pc) / height; /* Mass source term, unit: kg/m³/s */ 
    } 

    dS[eqn] = 0.0; 
    return source; 
} 
