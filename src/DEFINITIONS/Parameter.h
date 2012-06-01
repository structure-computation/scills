#include <vector>

#include "../UTILS/Sc2String.h"
#include "../../LMT/include/containers/vec.h"
#include "../../LMT/include/codegen/codegen.h"

using namespace LMT;
using namespace Codegen;

struct ParameterController;

class Parameter{
private:
    struct ParameterController{
        std::vector<Parameter&> multiresolution_parameters;
        std::vector<Ex>  multiresolution_symbols;
        Ex::MapExNum  multiresolution_values;
        
        std::vector<Parameter&> time_parameters;
        std::vector<Ex>  time_symbols;
        Ex::MapExNum  time_values;
        
        std::vector<Parameter&> space_parameters;
        std::vector<Ex>  space_symbols;
        Ex::MapExNum  space_values;
        
        ParameterController();
    };
    static ParameterController &controller;
    
    Sc2String __symbol;   /// Identifiant
    Sc2String __alias;    /// Nom utilisateur
    Type __type;          /// Type du parametre
    Ex __ex;              /// Expression associee
    std::vector<Ex> __expressions; /// Expressions utilisateur
    TYPEREEL __value;     /// Valeur numerique courante
    bool __has_changed;   /// Indique si la valeur du parametre a change
    
public:
    enum Type{Multiresolution,Time,Space};
    
    Parameter(Sc2String _symbol,Type _type);
    ~Parameter();
    
    void set_expressions(const Vec<Sc2String> &expr,std::vector<Ex>  &symbols);
    void update_value(unsigned i_step,TYPEREEL value);
    void has_changed() const;
    
    /*static void update_multiresolution_parameters(const Process &process);
    static void update_time_parameters();
    static void update_space_parameters();*/
    
    /// Operateurs permettant de manipuler comme des TYPEREEL (Exemples: cout << E; E = 200000; mu = E/(2.+2.*nu);)
    operator TYPEREEL() const;
    TYPEREEL operator=(TYPEREEL value);
    
    static TYPEREEL t;
    static TYPEREEL x;
    static TYPEREEL y;
    static TYPEREEL z;
};
