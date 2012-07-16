/*
#include "Parameter.h"


/// Constructeur. 'symbol_' est la chaine de caractere associe au parametre et 'owner_' son groupe (aucun i.e. 0 par defaut)
MainParameter(const Sc2String& symbol_ = "",ParameterGroup* owner_ = 0){
    
}
/// Constructeur de copie (si 'parameter' appartient a un ParameterGroup, il en sera retirer et remplacer par sa copie)
MainParameter(const MainParameter& parameter){}
/// Destructeur
~MainParameter(){}
/// Fonction de desallocation memoire
void free(){}
/// Reparente le parametre au ParameterGroup passe en argument
void setOwner(ParameterGroup* owner_){}
/// Pour modifier la valeur stockee
TYPEREEL operator=(TYPEREEL val){}
/// Permet de simuler un comportement de nombre reel (Ex: TYPEREEL a = 2*E;)
operator TYPEREEL() const{}

/// Fonction de debug
void affiche() const{}
*/
/****************************************************************************/
/*                             PARAMETER                                    */
/****************************************************************************/
/*
/// Constructeur
Parameter::Parameter(const Sc2String& symbol_,ParameterGroup* owner_):
symbol(symbol_),
self_ex(symbol_.c_str()),
main_parameter(main_parameter_)
{
    setOwner(owner_);
    
    //if(symbol_ != "PI"){
        //PI.setOwner(&ParameterGroup::global_constants);
        //PI.value = M_PI;
    //}Je ne savais pas où le mettre... 
}

/// Constructeur de copie
Parameter::Parameter(const Parameter& parameter):
symbol(parameter.symbol),
self_ex(parameter.self_ex),
main_parameter(parameter.main_parameter)
{
    setOwner(parameter.owner);
}

/// Destructeur
Parameter::~Parameter(){
    free();
}

void Parameter::free()
{
    for(unsigned i = 0; i < expressions.size(); i++)
        expressions[i].clear();
    expressions.free();
    for(unsigned i = 0; i < expr.size(); i++)
        if(expr[i].op != 0)
            delete expr[i].op;
    expr.free();
}

/// Reparente le parametre au ParameterGroup passe en argument
void Parameter::setOwner(ParameterGroup* owner_){
    if(owner != 0 and owner!=owner_){
        
    }
    if(owner!=owner_){
        owner = owner_;
        if(owner != 0){
            owner->addParameter(this);
        }
    }
}

/// Ajoute une expression au parametre
void Parameter::add_expression(const Sc2String& expr_){
    expressions.push_back(expr_);
}

/// Transforme les chaines de caracteres en expressions "pretes-a-evaluer"
void Parameter::prepare_expressions(std::vector<Codegen::Ex> &symbols_){
    const int n = expressions.size();
    expr.resize(n);
    for(int i = 0; i < n; i++){
        expr[i] = Codegen::read_ex(expressions[i],symbols_);
    }
}

/// Actualise la valeur du parametre a partir de l'expression numéro 'step', connaissant la carte de valeurs 'values'
TYPEREEL Parameter::update_value(int step_,Codegen::Ex::MapExNum &values_){
    value = (TYPEREEL)expr[step_].subs_numerical(values_,true);
    return value;
}

/// Permet de simuler un comportement de nombre reel (Ex: TYPEREEL a = 2*E;)
Parameter::operator TYPEREEL() const{
    return value;
}

/// Fonction de debug
void Parameter::affiche() const {
    std::cerr << "********************************************************************" << std::endl;
    std::cerr << "symbol : " << symbol << std::endl;
    std::cerr << "value : " << value << std::endl;
    std::cerr << "main_parameter : " << main_parameter << std::endl;
    std::cerr << "owner : " << (unsigned long) owner << std::endl;
    std::cerr << "expressions : " << expressions.size() << std::endl;
    for(int i = 0; i < expressions.size(); i++)
        std::cerr << i << " : " << expressions[i] << std::endl;
    std::cerr << "********************************************************************" << std::endl << std::endl;
}
*/
/****************************************************************************/
/*                         PARAMETERGROUP                                   */
/****************************************************************************/
/*
ParameterGroup::ParameterGroup():
parents(),parameters(),main_parameters()
{
    //if(this != &global_constants)
    //    add_parent(&global_constants);
}

ParameterGroup::~ParameterGroup(){
    free();
}

void ParameterGroup::free(){
}

void ParameterGroup::add_parent(ParameterGroup* group_){
    parents.push_back(group_);
}

Vec<ParameterGroup*> ParameterGroup::getAllParents(){
    Vec<ParameterGroup*> all_parents = parents;
    for(int i = 0; i < parents.size(); i++){
        Vec<ParameterGroup*> grand_parents = parents[i]->getAllParents();
        for(int j = 0; j < grand_parents.size(); j++){
            const int k_max = all_parents.size();
            int k;
            for(k = 0; k < k_max; k++){
                if(grand_parents[j] == all_parents[k]){
                    break;
                }
            }
            if(k!=k_max){ /// Si le parent n'est pas deja dans la liste...
                all_parents.push_back(grand_parents[j]);
            }
        }
    }
    return all_parents;
}

void ParameterGroup::addParameter(Parameter* param_){
    parameters.push_back(param_);
}

void ParameterGroup::addMainParameter(MainParameter* param_){
    main_parameters.push_back(param_);
}

Vec<Parameter*> ParameterGroup::getParentsParameters(){
    Vec<Parameter*> tmp_param;
    Vec<ParameterGroup*> all_parents = getAllParents();
    for(int i = 0; i < all_parents.size(); i++){
        tmp_param.append(all_parents[i]->main_parameters);
        tmp_param.append(all_parents[i]->parameters);
        tmp_param.append(all_parents[i]->getParentsParameters());
    }
    return tmp_param;
}

Codegen::Ex::MapExNum ParameterGroup::getValuesMap(){
    Codegen::Ex::MapExNum values;
    /// Recuperation des ParameterS du parent
    Vec<Parameter*> parent_parameters = getParentsParameters();
    for(int i = 0; i < parent_parameters.size(); i++){
        values[parent_parameters[i]->self_ex] = parent_parameters[i]->value;
    }
    return values;
}

/// Prepare les expressions de tous les parametres (non principaux) du groupe (c.f. Parameter::prepare_expressions)
void ParameterGroup::prepare_all(){
    std::vector<Codegen::Ex> symbols;
    /// Recuperation des ParameterS du parent
    Vec<Parameter*> parent_parameters = getParentsParameters();
    for(int i = 0; i < parent_parameters.size(); i++){
        symbols.push_back(parent_parameters[i]->self_ex);
    }
    /// Ajout des parametres principaux
    for(int i = 0; i < main_parameters.size(); i++){
        symbols.push_back(main_parameters[i]->self_ex);
    }
    /// Preparation des parametres
    for(int i = 0; i < parameters.size(); i++){
        parameters[i]->prepare_expressions(symbols);
    }
}

/// Met a jour tous les parametres (non principaux) du groupe (c.f. Parameter::update_value). 'step' indique l'index des expressions a evaluer.
void ParameterGroup::update_all(int step){
    affiche();
    Codegen::Ex::MapExNum values;
    /// Recuperation des ParameterS du parent
    Vec<Parameter*> parent_parameters = getParentsParameters();
    for(int i = 0; i < parent_parameters.size(); i++){
        values[parent_parameters[i]->self_ex] = parent_parameters[i]->value;
    }
    /// Ajout des parametres principaux
    for(int i = 0; i < main_parameters.size(); i++){
        values[main_parameters[i]->self_ex] = main_parameters[i]->value;
    }
    /// Preparation des parametres
    for(int i = 0; i < parameters.size(); i++){
        parameters[i]->update_value(step,values);
    }
}

void ParameterGroup::affiche() const {
    std::cerr << "********************************************************************" << std::endl;
    std::cerr << "Parents : " << parents.size() << std::endl;
    for(int i = 0; i < parents.size(); i++)
        std::cerr << i << " : " << (unsigned long) parents[i] << std::endl;
    std::cerr << "Main parameters : " << main_parameters.size() << std::endl;
    for(int i = 0; i < main_parameters.size(); i++)
        std::cerr << i << " : " << (unsigned long) main_parameters[i] << "(" << main_parameters[i]->symbol << "," << main_parameters[i]->value << ")" << std::endl;
    std::cerr << "Parameters : " << parameters.size() << std::endl;
    for(int i = 0; i < parameters.size(); i++)
        std::cerr << i << " : " << (unsigned long) parameters[i] << "(" << parameters[i]->symbol << "," << parameters[i]->value << ")" << std::endl;
    std::cerr << "********************************************************************" << std::endl << std::endl;
}


*/
/****************************************************************************/
/*                             ATTRIBUTS STATIQUES                          */
/****************************************************************************/
//ParameterGroup ParameterGroup::global_constants;    // A REVOIR
//MainParameter MainParameter::PI("PI");                      // A REVOIR
