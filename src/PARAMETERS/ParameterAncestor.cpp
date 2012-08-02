#include "ParameterAncestor.h"


/***********************************************************************************************
 *                                     ParameterAncestor                                       *
 ***********************************************************************************************/

/// Constructeur: 'symbol_' = le symbol (Sc2String et Codegen::Ex) associe au parametre, 'group_' = le ParameterGroup le "possedant" (0 par defaut i.e. aucun)
ParameterAncestor::ParameterAncestor(const Sc2String& symbol_, ParameterGroup* group_):
symbol(symbol_),
group(group_),
self_ex(Codegen::symbol(symbol_.c_str()))
{
}

/// Retourne le groupe du parametre (0 si aucun)
ParameterGroup* ParameterAncestor::getGroup() const{
    return group;
}

/// Accesseur de value_changed: Indique si la valeur du parametre a ete modifiee (variation relative de 1e-6 par rapport a la valeur precedente)
bool ParameterAncestor::hasChanged() const{
    return value_changed;
}

/// Pour qu'il se comporte comme un Scalar dans les calculs (Ex: Scalar mu = 0.5*E/(1.0+nu))
ParameterAncestor::operator Scalar() const {return value;}


/// Pour detecter les variations du parametre (utile pour la maj des operateurs)
Scalar ParameterAncestor::setValue(Scalar val){
    value_changed = std::abs(value-val)/value < epsilon();
    value = val;
}



/***********************************************************************************************
 *                                       MainParameter                                         *
 ***********************************************************************************************/

/// Constructeur: voir la Constructeur de ParameterAncestor
MainParameter::MainParameter(const Sc2String& symbol_, ParameterGroup* group_):
ParameterAncestor(symbol_,group_)
{
    if(group != 0){
        group->addParameter(this);
    }
}

/// permet de definir le groupe du parametre apres l'instanciation
void MainParameter::setGroup(ParameterGroup *group_){
    group = group_;
    group->addParameter(this);
}

/// Surcharge pour modifier l'attribut 'value'
Scalar MainParameter::operator=(Scalar val){
    setValue(val);
}

/// Affichage de debug dans std::cout
void MainParameter::affiche(){}



/***********************************************************************************************
 *                                       UserParameter                                         *
 ***********************************************************************************************/

/// Constructeur: voir la Constructeur de ParameterAncestor
UserParameter::UserParameter(const Sc2String& symbol_, ParameterGroup* group_):
ParameterAncestor(symbol_,group_),
str_expr("0"),
is_a_number(false)
{
    if(group != 0){
        group->addParameter(this);
    }
}

/// permet de definir le groupe du parametre apres l'instanciation
void UserParameter::setGroup(ParameterGroup *group_){
    group = group_;
    group->addParameter(this);
}

/// Stocke temporairement l'expression 'expression_' dans l'attribut 'str_expr'
void UserParameter::setExpression(const Sc2String& expression_){
    if(expression_ == ""){
        str_expr = "0";
    } else {
        str_expr = expression_;
    }
}

/// Transforme l'attribut 'srt_expr'(Sc2String) en 'ex_expr'(Codegen::Ex) connaissant 'symbols_'. Teste si l'expression est un nombre
void UserParameter::prepareExpression(std::vector<Codegen::Ex> &symbols_){
    ex_expr = Codegen::read_ex(str_expr,symbols_);
    is_a_number = ex_expr.is_a_number();
}

/// Evalue l'attribut 'ex_expr'(Codegen::Ex) connaissant 'values_'. Retourne le resultat de l'evaluation
Scalar UserParameter::updateValue(Codegen::Ex::MapExNum &values_){
    setValue((Scalar) ex_expr.subs_numerical(values_));
    return value;
}

/// Equivalent a 'setExpression(expression)'
void UserParameter::operator=(const Sc2String &expression){
    setExpression(expression);
}

/// Affichage pour le debug
void UserParameter::affiche(){}



/***********************************************************************************************
 *                                       ParameterGroup                                        *
 ***********************************************************************************************/

/// Creer un lien de dependance vers parent
void ParameterGroup::addParent(ParameterGroup* parent){
    parents.push_back(parent);
}

/// Ajoute un parametre de controle au group
void ParameterGroup::addParameter(MainParameter* param_){
    if(param_->getGroup() != this){
        param_->setGroup(this);
    }
    for(std::vector<MainParameter*>::iterator p = main_parameters.begin(); p!= main_parameters.end(); p++){
        if((*p) == param_){
            return ;    /// Le parametre fais deja parti du groupe
        }
    }
    main_parameters.push_back(param_);
}

/// Ajoute un parametre utilisateur au groupe
void ParameterGroup::addParameter(UserParameter* param_){
    if(param_->getGroup() != this){
        param_->setGroup(this);
        user_parameters.push_back(param_);
    }
    for(std::vector<UserParameter*>::iterator p = user_parameters.begin(); p!= user_parameters.end(); p++){
        if((*p) == param_){
            return ;    /// Le parametre fais deja parti du groupe
        }
    }
    user_parameters.push_back(param_);
}

/// Renvoi la liste de tous les parents, grand-parents, ... du groupe
std::vector<ParameterGroup*> ParameterGroup::getAllParents(){
    std::vector<ParameterGroup*> all_parents = parents;
    const int i_max = parents.size();
    for(int i = 0; i < i_max; i++){
        std::vector<ParameterGroup*> grand_parents = parents[i]->getAllParents();
        const int j_max = grand_parents.size();
        for(int j = 0; j < j_max; j++){
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

/// Renvoi la liste de tous les parametres contenus dans les parents, grand-parents,... du groupe
std::vector<ParameterAncestor*> ParameterGroup::getParentsParameters(){
    std::vector<ParameterAncestor*> tmp_param;
    std::vector<ParameterGroup*> all_parents = getAllParents();
    for(int i = 0; i < all_parents.size(); i++){
        const std::vector<MainParameter*> &mains = all_parents[i]->main_parameters;
        const std::vector<UserParameter*> &users = all_parents[i]->user_parameters;
        const std::vector<ParameterAncestor*> &others = all_parents[i]->getParentsParameters();
        const int tmp_size = tmp_param.size();
        const int mains_size = mains.size();
        const int users_size = users.size();
        const int others_size = others.size();
        tmp_param.reserve(tmp_size+mains_size+users_size+others_size);
        for(int i = 0; i < mains_size; i++){
            tmp_param.push_back(mains[i]);
        }
        for(int i = 0; i < users_size; i++){
            tmp_param.push_back(users[i]);
        }
        for(int i = 0; i < others_size; i++){
            tmp_param.push_back(others[i]);
        }
    }
    return tmp_param;
}

std::vector<Codegen::Ex> ParameterGroup::getParentsSymbols(){
    std::vector<Codegen::Ex> symbols;
    std::vector<ParameterAncestor*> parent_parameters = getParentsParameters();
    symbols.reserve(parent_parameters.size());
    for(std::vector<ParameterAncestor*>::iterator i = parent_parameters.begin(); i != parent_parameters.end(); i++){
        symbols.push_back((*i)->self_ex);
    }
    return symbols;
}

/// Prepare les expressions de tous les parametres utilisateur du groupe en vue de leur evaluation
void ParameterGroup::prepareParameters(){
    std::vector<Codegen::Ex> symbols = getParentsSymbols();
    symbols.reserve(symbols.size() + main_parameters.size());
    for(std::vector<MainParameter*>::iterator i = main_parameters.begin(); i != main_parameters.end(); i++){
        symbols.push_back((*i)->self_ex);
    }
    for(std::vector<UserParameter*>::iterator i = user_parameters.begin(); i != user_parameters.end(); i++){
        (*i)->prepareExpression(symbols);
    }
}

/// Retourne la table qui va permettre d'evaluer les expressions des parametres utilisateur. Ne contient PAS les parametres de controle (pour assignation manuelle)
Codegen::Ex::MapExNum ParameterGroup::getParentsValues(){
    Codegen::Ex::MapExNum values;
    std::vector<ParameterAncestor*> parent_parameters = getParentsParameters();
    for(int i = 0; i < parent_parameters.size(); i++){
        values[parent_parameters[i]->self_ex] = parent_parameters[i]->value;
    }
    return values;
}

/// Met a jour tous les parametres utilisateur du groupe
void ParameterGroup::updateParameters(){
    Codegen::Ex::MapExNum values = getParentsValues();
    for(std::vector<MainParameter*>::iterator i = main_parameters.begin(); i != main_parameters.end(); i++){
        values[(*i)->self_ex] = (*i)->value;
    }
    for(std::vector<UserParameter*>::iterator i = user_parameters.begin(); i != user_parameters.end(); i++){
        (*i)->updateValue(values);
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
    std::cerr << "Parameters : " << user_parameters.size() << std::endl;
    for(int i = 0; i < user_parameters.size(); i++)
        std::cerr << i << " : " << (unsigned long) user_parameters[i] << "(" << user_parameters[i]->symbol << "," << user_parameters[i]->value << ")" << std::endl;
    std::cerr << "********************************************************************" << std::endl << std::endl;
}