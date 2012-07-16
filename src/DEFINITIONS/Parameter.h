#ifndef PARAMETER_H
#define PARAMETER_H
/*
#include <vector>

#include "../all_classes.h"
#include "../../LMT/include/codegen/codegen.h"
#include "../../LMT/include/containers/basicops.h"

using namespace Codegen;


struct ParameterGroup;

struct MainParameter{
public:
    /// Constructeur. 'symbol_' est la chaine de caractere associe au parametre, 'owner_' son groupe (aucun i.e. 0 par defaut)
    MainParameter(const Sc2String& symbol_ = "",ParameterGroup* owner_ = 0);
    /// Constructeur de copie (si 'parameter' appartient a un ParameterGroup, il en sera retirer et remplacer par sa copie)
    MainParameter(const MainParameter& parameter);
    /// Destructeur
    ~MainParameter();
    /// Fonction de desallocation memoire
    void free();
    /// Reparente le parametre au ParameterGroup passe en argument
    void setOwner(ParameterGroup* owner_);
    /// Pour modifier la valeur stockee
    TYPEREEL operator=(TYPEREEL val);
    /// Permet de simuler un comportement de nombre reel (Ex: TYPEREEL a = 2*E;)
    operator TYPEREEL() const;
    
    /// Fonction de debug
    void affiche() const;
    
    Sc2String symbol;           /// Symbole du parametre
    Codegen::Ex self_ex;        /// Expression associe au parametre lui-meme
    Vec<Sc2String> expressions; /// Expressions (sous forme de chaines, utile pour le debug)
    Vec<Codegen::Ex> expr;      /// Expressions (sous forme pretraitee) du parametre a chaque step (multiresolution ou temporel)
    TYPEREEL value;             /// Valeur courante du parametre
    bool main_parameter;
    ParameterGroup* owner;
    
    //static Parameter PI;    /// Represente la constante Pi = 3.14159... (voir macro M_PI dans math.h)
};

struct Parameter{
public:
    /// Constructeur. 'symbol_' est la chaine de caractere associe au parametre et 'owner_' son groupe (aucun i.e. 0 par defaut)
    Parameter(const Sc2String& symbol_ = "",ParameterGroup* owner_ = 0);
    /// Constructeur de copie (si 'parameter' appartient a un ParameterGroup, il en sera retirer et remplacer par sa copie)
    Parameter(const Parameter& parameter);
    /// Destructeur
    ~Parameter();
    /// Fonction de desallocation memoire
    void free();
    /// Reparente le parametre au ParameterGroup passe en argument
    void setOwner(ParameterGroup* owner_);
    /// Ajoute une expression au parametre
    void add_expression(const Sc2String& expr_);
    /// Transforme les chaines de caracteres en expressions "pretes-a-evaluer"
    void prepare_expressions(std::vector<Codegen::Ex> &symbols_);
    /// Actualise la valeur du parametre a partir de l'expression numéro 'step', connaissant la carte de valeurs 'values'
    TYPEREEL update_value(int step_,Codegen::Ex::MapExNum &values_);
    /// Permet de simuler un comportement de nombre reel (Ex: TYPEREEL a = 2*E;)
    operator TYPEREEL() const;
    
    /// Fonction de debug
    void affiche() const;
    
    Sc2String symbol;           /// Symbole du parametre
    Codegen::Ex self_ex;        /// Expression associe au parametre lui-meme
    Vec<Sc2String> expressions; /// Expressions (sous forme de chaines, utile pour le debug)
    Vec<Codegen::Ex> expr;      /// Expressions (sous forme pretraitee) du parametre a chaque step (multiresolution ou temporel)
    TYPEREEL value;             /// Valeur courante du parametre
    bool main_parameter;
    ParameterGroup* owner;
};

struct ParameterGroup{
public:
    ParameterGroup();                               /// Constructeur
    ~ParameterGroup();                              /// Destructeur
    void free();                                    /// Fonction de desallocation memoire
    
    void add_parent(ParameterGroup* group_);
    Vec<ParameterGroup*> getAllParents();       /// Retrouve tous les parents et parents de parents... (elimine les doublons)
    
    void addParameter(Parameter* param_);       /// Ajoute 'param_' a la liste des parametres a evaluer
    void addMainParameter(MainParameter* param_);   /// Ajoute 'param_' a la liste des parametres de control
    Vec<Parameter*> getParentsParameters();     /// Retourne un Vec contenant les pointeurs de tous les parametres des parents du groupe.
    
    //std::vector<Codegen::Ex> getSymbolsList();  /// Retourne la liste des symboles PAS UTILE JE PENSE
    Codegen::Ex::MapExNum getValuesMap(); /// Retourne la carte des valeurs NE CONTENANT PAS les parametres de controle, utile pour les calculs de CL
    void prepare_all();         /// Prepare les expressions de tous les parametres (non principaux) du groupe (c.f. Parameter::prepare_expressions)
    void update_all(int step);  /// Met a jour tous les parametres (non principaux) du groupe (c.f. Parameter::update_value). 'step' indique l'index des expressions a evaluer.
    
    void affiche() const; /// Fonction de debug
    
    Vec<ParameterGroup*> parents;       /// Pointeurs vers les parents
    Vec<MainParameter*> main_parameters;    /// Contient les parametres de controles (x, y, z, t ou m), dont les valeurs sont directement fournises
    Vec<Parameter*> parameters;         /// Contient les parametres a evaluer
    
    //static ParameterGroup global_constants; /// Regroupe toutes les constantes globales. Tous les ParameterGroupS en l'ont en parent
};
*/
#endif //PARAMETER_H