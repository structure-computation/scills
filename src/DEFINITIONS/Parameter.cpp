#include "Parameter.h"
#include "../../../LMTpp/include/containers/mat.h"

ParameterController &Parameter::controller = ParameterController();
Parameter Parameter::t("t",Parameter::Time);
Parameter Parameter::x("x",Parameter::Space);
Parameter Parameter::y("y",Parameter::Space);
Parameter Parameter::z("z",Parameter::Space);

Parameter(Sc2String _symbol, Parameter::Type type):
__symbol(_symbol),
__type(type),
__ex(symbol(_symbol)),
__has_changed(true)
{
    switch(__type){
        case Multiresolution:
            controller.multiresolution_parameters.push_back(*this);
            controller.multiresolution_symbols.push_back(__ex);
            controller.multiresolution_values[__ex] = 0.0;
            break;
        case Time:
            controller.time_parameters.push_back(*this);
            controller.time_symbols.push_back(__ex);
            controller.time_values[__ex] = 0.0;
            break;
        case Space:
            controller.space_parameters.push_back(*this);
            controller.space_symbols.push_back(__ex);
            controller.space_values[__ex] = 0.0;
            break;
        default:
            std::cerr << "Type de Parameter invalid pour " << __symbol << std::endl;
            assert(0);
    }
}

~Parameter(){}

void Parameter::set_expressions(const Vec<Sc2String> &expr,std::vector<Ex>  &symbols){
    const unsigned n = expr.size();
    for(unsigned i = 0; i < n; i++){
        __expressions.push_back(read_ex(expr[i],symbols));
    }
}

void Parameter::update_value(){
    if (__expressions.size() == 0) return;
    unsigned i_step;
    TYPEREEL new_value;
    switch(__type){
        case Multiresolution:
            new_value = (TYPEREEL) __expressions[i_step].subs_numerical(controller.multiresolution_values);
            break;
        case Time:
            new_value = (TYPEREEL) __expressions[i_step].subs_numerical(controller.time_values);
            break;
        case Space:
            new_value = (TYPEREEL) __expressions[i_step].subs_numerical(controller.space_values);
            break;
        default:
            std::cerr << "Type de Parameter invalid pour " << __symbol << std::endl;
            assert(0);
    }
    if(new_value != __value){
        __has_changed = true;
        __value = new_value;
    } else {
        __has_changed = false;
    }
}

void Parameter::has_changed() const{
    return __has_changed;
}

Parameter::operator TYPEREEL() const{
    return __value;
}

TYPEREEL Parameter::operator=(TYPEREEL value){
    if(value != __value){
        __has_changed = true;
        __value = value;
    } else {
        __has_changed = false;
    }
}
/*
void Parameter::update_multiresolution_parameters(const Process &process){
    
}

void Parameter::update_time_parameters(const Process &process){
    
}

void Parameter::update_space_parameters(const Process &process){
    
}
*/