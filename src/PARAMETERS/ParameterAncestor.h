#ifndef PARAMETERANCESTOR
#define PARAMETERANCESTOR


#include <vector>
#include "../DEFINITIONS/main_typedef.h"

#include "../../LMT/include/codegen/codegen.h"
#include "../../LMT/include/containers/basicops.h"
using namespace Codegen;

class ParameterGroup;

/** Classe parente (et abstraite) des MainParameterS et UserParameterS
 * Elle definie :
 *     - ParameterGroup* group     : Le groupe auquel il appartient
 *     - const Sc2String symbol    : Le symbol associe
 *     - const Codegen::Ex self_ex : L'expression associe au symbol
 *     - TYPEREEL value;           : La valeur courante
 *     - bool value_changed;       : Indique si la valeur du parametre a ete modifiee (variation relative de 1e-6 par rapport a la valeur precedente)
 */
class ParameterAncestor{
    friend class ParameterGroup;  /// Pour facilite les acces
public:
    /// Constructeur: 'symbol_' = le symbol (Sc2String et Codegen::Ex) associe au parametre, 'group_' = le ParameterGroup le "possedant" (0 par defaut i.e. aucun)
    ParameterAncestor(const Sc2String& symbol_, ParameterGroup* group_ = 0);
    
    /// Retourne le groupe du parametre (0 si aucun)
    ParameterGroup* getGroup() const;
    
    /// Accesseur de value_changed: Indique si la valeur du parametre a ete modifiee (variation relative de 1e-6 par rapport a la valeur precedente)
    bool hasChanged() const;
    
    /// Pour qu'il se comporte comme un TYPEREEL dans les calculs (Ex: TYPEREEL mu = 0.5*E/(1.0+nu))
    operator TYPEREEL() const;
    
    /// Pour le debug et pour s'assurer que seules les classes filles sont instanciables
    virtual void affiche() = 0;
    
//protected:
    /// Pour detecter les variations du parametre (utile pour la maj des operateurs)
    TYPEREEL setValue(TYPEREEL val);
    /// Voir setValue
    static TYPEREEL epsilon(){return 1.e-6;}
    
    ParameterGroup* group;      /// Le groupe auquel il appartient
    Sc2String symbol;     /// Le symbol associe
    Codegen::Ex self_ex;  /// L'expression associe au symbol
    TYPEREEL value;             /// La valeur courante
    bool value_changed;         /// Indique si la valeur du parametre a ete modifiee
};


/** Classe de convenience (heritee de ParameterAncestor) pour definir les parametres de controle des expressions : t, x, y, ...
 * La valeur stockee peut etre modifiee via l'operateur = (assignation)
 * Ex : t = 10.0;
 */
class MainParameter: public ParameterAncestor{
public:
    /// Constructeur: voir la Constructeur de ParameterAncestor
    MainParameter(const Sc2String& symbol_, ParameterGroup* group_ = 0);
    
    /// permet de definir le groupe du parametre apres l'instanciation
    void setGroup(ParameterGroup *group_);
    
    /// Surcharge pour modifier l'attribut 'value'
    TYPEREEL operator=(TYPEREEL val);
    
    /// Affichage de debug dans std::cout
    void affiche();
};


/** Classe fille de ParameterAncestor pour manipuler les expressions Codegen::Ex plus facilement :
 * On passe une chaine de caracteres a la fonction setExpression ou par l'operateur =
 * On transforme cette chaine en Codegen::Ex (prepareExpression).
 * On l'evalue (updateValue).
 */
class UserParameter: public ParameterAncestor{
public:
    /// Constructeur: voir la Constructeur de ParameterAncestor
    UserParameter(const Sc2String& symbol_ = "", ParameterGroup* group_ = 0);
    
    /// permet de definir le groupe du parametre apres l'instanciation
    void setGroup(ParameterGroup *group_);
    
    /// Stocke temporairement l'expression 'expression_' dans l'attribut 'str_expr'
    void setExpression(const Sc2String& expression_);
    
    /// Transforme l'attribut 'srt_expr'(Sc2String) en 'ex_expr'(Codegen::Ex) connaissant 'symbols_'
    void prepareExpression(std::vector<Codegen::Ex> &symbols_);
    
    /// Evalue l'attribut 'ex_expr'(Codegen::Ex) connaissant 'values_'. Retourne le resultat de l'evaluation
    Scalar updateValue(Codegen::Ex::MapExNum &values_);
    
    /// Equivalent a 'setExpression(expression)'
    void operator=(const Sc2String &expression);
    
    /// Affichage pour le debug
    void affiche();
    
//private:
    Sc2String str_expr;     /// L'expression du parametre sous forme de chaine de caracteres
    Codegen::Ex ex_expr;    /// L'expression du parametre sous forme de Codegen::Ex (a lire apres appel a prepareExpression)
    bool is_a_number;       /// Indique si le parametre est une valeur constante (mis a jour par prepareExpression)
};


/**Cette classe defini des containers pour les classes MainParameter (parametres de controle) et UserParameter (parametres utilisateur)
 * Des methode de groupe sont egalement presentent pour faciliter la mise a jour des parametres utilisateur:
 *     - Forme automatisee:
 *         void addParent(ParameterGroup*): on cree des lien de dependance entre les groupes de parametres
 *         void prepareParameters(): preparation des parametres utilisateur (c.f. UserParameter::prepareExpression)
 *         void updateParameters(): mise a jour des parametres utilisateur (c.f. UserParameter::updateValue)
 *     - Forme personnalisee (en cas de parametre dans aucun groupe par exemple):
 *         void addParent(ParameterGroup*): on cree des lien de dependance entre les groupes de parametres
 *         void getParentsSymbols: on recupere la liste des symbols de tous les parametre dans les parents a passer en parametre de UserParameter::prepareExpression
 *         void getParentsValues: idem pour la table de correspondance a fournir a UserParameter::updateValue
 * 
 */
class ParameterGroup{
public:
    /// Creer un lien de dependance vers parent
    void addParent(ParameterGroup* parent);
    
    /// Ajoute un parametre de controle au group
    void addParameter(MainParameter* param_);
    
    /// Ajoute un parametre utilisateur au groupe
    void addParameter(UserParameter* param_);
    
    /// Renvoi la liste de tous les parents, grand-parents, ... du groupe
    std::vector<ParameterGroup*> getAllParents();
    
    /// Renvoi la liste de tous les parametres contenus dans les parents, grand-parents,... du groupe
    std::vector<ParameterAncestor*> getParentsParameters();
    
    std::vector<Codegen::Ex> getParentsSymbols();
    
    /// Prepare les expressions de tous les parametres utilisateur du groupe en vue de leur evaluation
    void prepareParameters();
    
    /// Retourne la table qui va permettre d'evaluer les expressions des parametres utilisateur. Ne contient PAS les parametres de controle (pour assignation manuelle)
    Codegen::Ex::MapExNum getParentsValues();
    
    /// Met a jour tous les parametres utilisateur du groupe
    void updateParameters();
    
    /// Affichage pour le debug
    void affiche() const;
    
//private:
    std::vector<ParameterGroup*> parents;
    std::vector<MainParameter*> main_parameters;
    std::vector<UserParameter*> user_parameters;
};

#endif //PARAMETERANCESTOR


