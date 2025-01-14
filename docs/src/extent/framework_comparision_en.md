# Comparing Programming Styles for Custom Models in superflexpy, MARRMoT, and HydroModels.jl

This content corresponds to Section 4.1 of the paper, where we discuss the differences in programming styles for custom model implementation across these three software libraries. Using the GR4J model implementation as an example, we compare the programming styles for customizing model computation formulas.

## [Superflexpy](https://github.com/dalmo1991/superflexPy/tree/master)

Superflexpy is a Python implementation based on the SUPERFLEX concept. For implementing hydrological model computations, this framework typically uses Python syntax to express model calculation formulas, as shown below.

```python
@staticmethod
def _flux_function_python(S, S0, ind, P, x1, alpha, beta, ni, PET, dt):
    if ind is None:
        return (
            [
                P * (1 - (S / x1) ** alpha), # Ps
                -PET * (2 * (S / x1) - (S / x1) ** alpha), # Evaporation
                -((x1 ** (1 - beta)) / ((beta - 1) * dt)) * (ni ** (beta - 1)) * (S**beta), # Perc
            ],
            0.0,
            S0 + P * (1 - (S / x1) ** alpha) * dt,
        )
    else:
        return (
            [
                P[ind] * (1 - (S / x1[ind]) ** alpha[ind]), # Ps
                -PET[ind] * (2 * (S / x1[ind]) - (S / x1[ind]) ** alpha[ind]), # Evaporation
                -((x1[ind] ** (1 - beta[ind])) / ((beta[ind] - 1) * dt[ind]))
                * (ni[ind] ** (beta[ind] - 1))
                * (S ** beta[ind]), # Perc
            ],
            0.0,
            S0 + P[ind] * (1 - (S / x1[ind]) ** alpha[ind]) * dt[ind],
            [
                -(P[ind] * alpha[ind] / x1[ind]) * ((S / x1[ind]) ** (alpha[ind] - 1)),
                -(PET[ind] / x1[ind]) * (2 - alpha[ind] * ((S / x1[ind]) ** (alpha[ind] - 1))),
                -beta[ind]
                * ((x1[ind] ** (1 - beta[ind])) / ((beta[ind] - 1) * dt[ind]))
                * (ni[ind] ** (beta[ind] - 1))
                * (S ** (beta[ind] - 1)),
            ],
        )

@staticmethod
@nb.jit(
    "Tuple((UniTuple(f8, 3), f8, f8, UniTuple(f8, 3)))"
    "(optional(f8), f8, i4, f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:])",
    nopython=True,
)
def _flux_function_numba(S, S0, ind, P, x1, alpha, beta, ni, PET, dt):
    return (
        (
            P[ind] * (1 - (S / x1[ind]) ** alpha[ind]),  # Ps
            -PET[ind] * (2 * (S / x1[ind]) - (S / x1[ind]) ** alpha[ind]),  # Evaporation
            -((x1[ind] ** (1 - beta[ind])) / ((beta[ind] - 1) * dt[ind]))
            * (ni[ind] ** (beta[ind] - 1))
            * (S ** beta[ind]),  # Perc
        ),
        0.0,
        S0 + P[ind] * (1 - (S / x1[ind]) ** alpha[ind]) * dt[ind],
        (
            -(P[ind] * alpha[ind] / x1[ind]) * ((S / x1[ind]) ** (alpha[ind] - 1)),
            -(PET[ind] / x1[ind]) * (2 - alpha[ind] * ((S / x1[ind]) ** (alpha[ind] - 1))),
            -beta[ind]
            * ((x1[ind] ** (1 - beta[ind])) / ((beta[ind] - 1) * dt[ind]))
            * (ni[ind] ** (beta[ind] - 1))
            * (S ** (beta[ind] - 1)),
        ),
    )
```

As shown, when implementing the GR4J production module calculations [origin code](https://github.com/dalmo1991/superflexPy/blob/master/superflexpy/implementation/elements/gr4j.py), the superflexpy framework requires separate construction of numpy and numba computation code. The function must specify hydrological fluxes `Ps`, `Evaporation`, `Perc`, and state variable balance equations. This approach couples multiple computational process implementations in a single code block, resulting in reduced code reusability (only reusable at the module level) and increased code volume and redundancy.

## [MARRMoT](https://github.com/wknoben/MARRMoT)

MARRMoT is a hydrological modeling framework developed in MATLAB. This framework constructs model types separately based on different hydrological models, storing the calculation formulas and state change equations involved. The simplified implementation code for the GR4J model in this framework is shown below.


```matlab
function [dS, fluxes] = model_fun(obj, S)
    % definitions
    ...
    % fluxes functions
    flux_pn   = max(P-Ep,0);
    flux_en   = max(Ep-P,0);
    flux_ef   = P - flux_pn;
    flux_ps   = saturation_4(S1,x1,flux_pn);
    flux_es   = evap_11(S1,x1,flux_en);
    flux_perc = percolation_3(S1,x1);
    flux_q9   = route(.9.*(flux_pn - flux_ps + flux_perc), uh_q9);
    flux_q1   = route(.1.*(flux_pn - flux_ps + flux_perc), uh_q1);
    flux_fr   = recharge_2(3.5,S2,x3,x2);
    flux_fq   = flux_fr;
    flux_qr   = baseflow_3(S2,x3);
    flux_qt   = flux_qr + max(flux_q1 + flux_fq,0);
    % this flux is not included in original MARRMoT,
    % but it is useful to calculate the water balance
    flux_ex = flux_fr + max(flux_q1 + flux_fq,0) - flux_q1;      
    % stores ODEs
    dS1 = flux_ps - flux_es - flux_perc;
    dS2 = flux_q9 + flux_fr - flux_qr;
    % outputs
    dS = [dS1 dS2];
    fluxes = [flux_pn,   flux_en, flux_ef, flux_ps, flux_es,...
                flux_perc, flux_q9, flux_q1, flux_fr, flux_fq,...
                flux_qr,   flux_qt, flux_ex];
end
```

As evident, the MARRMoT framework expresses each Flux calculation formula through functions, allowing users to directly call functions for calculations (according to provided documentation). This approach significantly simplifies code writing requirements compared to superflexpy. However, MARRMoT has poor support for module reusability and cannot be reused like superflexpy. Additionally, combining variable assignments, flux calculations, and state change equations in one function results in poor decoupling of computational content.

## HydroModels.jl

```julia
#* define the production store
prod_funcs = [
    HydroFlux([prcp, ep] => [pn, en], exprs=[prcp - min(prcp, ep), ep - min(prcp, ep)]),
    HydroFlux([pn, soilwater] => [ps], [x1], exprs=[max(0.0, pn * (1 - (soilwater / x1)^2))]),
    HydroFlux([en, soilwater] => [es], [x1], exprs=[en * (2 * soilwater / x1 - (soilwater / x1)^2)]),
    HydroFlux([soilwater] => [perc], [x1], exprs=[((x1)^(-4)) / 4 * ((4 / 9)^(4)) * (soilwater^5)]),
    HydroFlux([pn, ps, perc] => [pr], exprs=[pn - ps + perc]),
    HydroFlux([pr] => [slowflow, fastflow], exprs=[0.9 * pr, 0.1 * pr]),
]
prod_dfuncs = [StateFlux([ps] => [es, perc], soilwater)]
prod_ele = HydroBucket(name=:gr4j_prod, funcs=prod_funcs, dfuncs=prod_dfuncs)

#* uh function
uh_flux_1 = HydroModels.UnitHydrograph([slowflow] => [slowflow_routed], x4, uhfunc=HydroModels.UHFunction(:UH_1_HALF), solvetype=:SPARSE)
uh_flux_2 = HydroModels.UnitHydrograph([fastflow] => [fastflow_routed], x4, uhfunc=HydroModels.UHFunction(:UH_2_FULL), solvetype=:SPARSE)

#* define the routing store
rst_funcs = [
    HydroFlux([routingstore] => [exch], [x2, x3], exprs=[x2 * abs(routingstore / x3)^3.5]),
    HydroFlux([routingstore, slowflow, exch] => [routedflow], [x3], exprs=[x3^(-4) / 4 * (routingstore + slowflow + exch)^5]),
    HydroFlux([routedflow, fastflow_routed, exch] => [flow], exprs=[routedflow + max(fastflow_routed + exch, 0.0)]),
]
rst_dfuncs = [StateFlux([exch, slowflow] => [routedflow], routingstore)]
rst_ele = HydroBucket(name=:gr4j_rst, funcs=rst_funcs, dfuncs=rst_dfuncs)

gr4j_model = HydroModel(name=:gr4j, components=[prod_ele, uh_flux_1, uh_flux_2, rst_ele])
```

Returning to the GR4J model implementation in HydroModels.jl framework, we sequentially define the production store, unit hydrograph, and routing store. In this process, the model can achieve reusability while maintaining variable correspondence. The construction process includes instantiation of classes such as HydroFlux, StateFlux, UnitHydrograph, HydroModel, and HydroBucket. Unlike the previous two frameworks that use functional construction, these class instantiations express calculation formulas through symbolic programming. After model construction, calculations are performed through input data and parameters, demonstrating a design philosophy that completely decouples model structure and formulas from parameters and data, achieving high reusability.