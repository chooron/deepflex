# 导入模块
using LumpedHydro

model = LumpedHydro.ExpHydro.Node(name=:exphydro, mtk=true, step=false)

@info LumpedHydro.get_input_names(model)
@info LumpedHydro.get_output_names(model)
@info LumpedHydro.get_param_names(model)
@info LumpedHydro.get_state_names(model)

unit = model.units[:exphydro];
@info LumpedHydro.get_input_names(unit)
@info LumpedHydro.get_output_names(unit)
@info LumpedHydro.get_param_names(unit)
@info LumpedHydro.get_state_names(unit)

soil = unit.soil;
@info LumpedHydro.get_input_names(soil)
@info LumpedHydro.get_output_names(soil)
@info LumpedHydro.get_param_names(soil)
@info LumpedHydro.get_state_names(soil)