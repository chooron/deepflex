@reexport module SIMHYD
# https://github.com/hydrogo/LHMP/blob/master/models/simhyd_cemaneige.py
using ..HydroModels
using ..HydroModels: @variables,@parameters
using ..HydroModels.Symbolics: Num
"""
elements in SIMHYD
"""
function SIMHYD_ELE(; name::Symbol, mtk::Bool=true)
    @variables prcp = 0.0
    @variables pet = 0.0

    @variables U = 0.0
    @variables IMAX = 0.0
    @variables INT = 0.0
    @variables INR = 0.0
    @variables RMO = 0.0
    @variables IRUN = 0.0
    @variables EVAP = 0.0
    @variables SRUN = 0.0
    @variables REC = 0.0
    @variables SMF = 0.0
    @variables POT = 0.0
    @variables BAS = 0.0
    @variables GWF = 0.0
    @variables SMS = 0.0
    @variables GW = 0.0

    @parameters INSC = 0.0
    @parameters COEFF = 0.0
    @parameters SQ = 0.0
    
    @parameters SMSC = 0.0
    @parameters SUB = 0.0
    @parameters CRAK = 0.0
    @parameters K = 0.0
    @parameters etmul = 0.0
    @parameters DELAY = 0.0
    @parameters X_m = 0.0

    funcs = [
        SimpleFlux([pet] => [IMAX], [INSC], exprs=@.[min(INSC, pet)]),
        SimpleFlux([prcp, IMAX] => [INT], Num[], exprs=@.[min(IMAX, prcp)]),
        SimpleFlux([prcp, INT] => [INR], Num[], exprs=@.[prcp - INT]),
        SimpleFlux([pet, INT] => [POT], Num[], exprs=@.[pet - INT]),
        SimpleFlux([POT, SMS] => [EVAP], [SMSC], exprs=@.[min(10 * (SMS / SMSC), POT)]),
        SimpleFlux([INR, SMS] => [RMO], [COEFF, SQ, SMSC], exprs=@.[min(COEFF * exp(-SQ * (SMS / SMSC)), INR)]),
        SimpleFlux([INR, RMO] => [IRUN], Num[], exprs=@.[INR - RMO]),
        SimpleFlux([RMO, SMS] => [SRUN], [SUB, SMSC], exprs=@.[SUB * (SMS / SMSC) * RMO]),
        SimpleFlux([RMO, SRUN, SMS] => [REC], [CRAK, SMSC], exprs=@.[CRAK * (SMS / SMSC) * (RMO - SRUN)]),
        SimpleFlux([RMO, SRUN, REC] => [SMF], Num[], exprs=@.[RMO - SRUN - REC]),
        SimpleFlux([SMS, SMF] => [GWF], [SMSC], exprs=@.[ifelse(SMS == SMSC, SMF, 0.0)]),
        SimpleFlux([GW] => [BAS], [K], exprs=@.[K * GW]),
        SimpleFlux([IRUN, SRUN, BAS] => [U], [K], exprs=@.[IRUN + SRUN + BAS]),
    ]

    dfuncs = [
        StateFlux([SMF] => [EVAP, GWF], SMS),
        StateFlux([REC, GWF] => [BAS], GW)
    ]

    HydroBucket(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk,
    )
end

function Model(; name::Symbol, mtk::Bool=true)
    elements = [
        SIMHYD_ELE(name=name, mtk=mtk),
    ]

    HydroModel(
        name,
        components=elements,
    )
end
end