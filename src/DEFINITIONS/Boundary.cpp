#include "Boundary.h"
#include "../UTILITAIRES/utilitaires.h"


std::vector<Ex> Boundary::time_symbols;    ///< Vecteur contenant les symboles temporels
Ex::MapExNum Boundary::time_values;        ///< Map< Symbol, Valeur > des grandeurs symboliques (autres que spatiales) du calcul
std::vector<Ex> Boundary::space_symbols;   ///< Vecteur contenant les symboles spatiaux
Ex::MapExNum Boundary::space_values;       ///< Map< Symbol, Valeur > des grandeurs symboliques spatiales du calcul


void Boundary::buildTimeSymbols(const DataUser &data_user) {
    if (time_symbols.size() != 0){
        //std::cout << "WARNING buildTimeSymbols a deja ete appele" << std::endl;
        return;
    }
    
    time_symbols.push_back(symbol("t"));
    if(data_user.options.Multiresolution_on==1){
        ///ajout des variables de multiresolution aux symboles
        for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
            Sc2String var="V"; var<<i_par; ///nom de la variable de multiresolution
            time_symbols.push_back(symbol(var.c_str()));
        }
    }
}

void Boundary::updateTimeValues(const Process &process,const DataUser &data_user) {
    unsigned i_step=process.temps->step_cur;
    TYPEREEL ti=process.temps->current_time;
    /// Recuperation de l'instant precedent (pour le calcul des derivees)
    TYPEREEL told = ti - process.temps->time_step[i_step].dt;
    
    time_values[time_symbols[0]]=ti;
    if(data_user.options.Multiresolution_on==1){
        /// Evaluation des variables de multiresolution
        for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
            time_values[time_symbols[1+i_par]] = data_user.Multiresolution_parameters[i_par].current_value;
    }
}

/// Construit la liste des symboles spatiaux
void Boundary::buildSpaceSymbols() {
    if (space_symbols.size() != 0){
        //std::cout << "WARNING buildTimeSymbols a deja ete appele" << std::endl;
        return;
    }
    #if DIM==2
    space_symbols.push_back(symbol("x"));
    space_symbols.push_back(symbol("y"));
    #elif DIM==3
    space_symbols.push_back(symbol("x"));
    space_symbols.push_back(symbol("y"));
    space_symbols.push_back(symbol("z"));
    #endif
}

/// Actualise les valeurs spatiales aux coordonnees du point "space_point"
void Boundary::updateSpaceValues(const Vec<TYPEREEL,DIM> &space_point) {
    for(int i = 0; i < DIM; i++)
        space_values[space_symbols[i]] = space_point[i];
}

/// Construit les listes de symboles permettant d'evaluer les CL
void Boundary::buildBoundarySymbols(const DataUser &data_user){
    buildTimeSymbols(data_user);
    buildSpaceSymbols();
}

void Boundary::free() {/*
    for(unsigned i = 0; i < time_symbols.size(); i++)
        if(time_symbols[i]!=0){
            delete time_symbols[i];
            time_symbols[i] = 0;
        }
        for(unsigned i = 0; i < space_symbols.size(); i++)
            if(space_symbols[i]!=0){
                delete space_symbols[i];
                space_symbols[i] = 0;
            }*/
}


/** Evalue la valeur de la CL au step numero "i_step" en chaque noeud de "nodeeq" et stocke le resultat dans V
 * Si un vecteur de normales equivalentes "neqs" est fourni, le resultat sera projete.
 * ATTENTION!!! Cette fonction ne reactualise pas les valeurs dependant du temps!
 */
void Boundary::evaluate(unsigned i_step, Vec<TYPEREEL> &V, Vec<Vec<TYPEREEL,DIM> > &nodeeq,const Vec<TYPEREEL> &neqs){
    Ex expr;
    /// Evaluation des fonctions temporelles
    expr = read_ex(fcts_temporelles[i_step],time_symbols);
    for( unsigned dir=0;dir<DIM ;dir++ ){
        ft[dir]=(TYPEREEL)expr.subs_numerical(time_values);
    }
    std::cout << "***********************************" << std::endl;
    std::cout << "comp : " << comp << std::endl;
    std::cout << "ti = " << time_values[time_symbols[0]] << std::endl;
    std::cout << "ft = " << ft[0] << " , " << ft[1] << " , " << ft[2] <<std::endl;
    std::cout << "***********************************" << std::endl;
    /// Evaluation des fonctions spatiales
    for(unsigned i_node=0;i_node<nodeeq.size();++i_node) {
        Vec<TYPEREEL, DIM> data;
        Vec<TYPEREEL,DIM> neq;
        if(neqs.size()>0)   /// Cas d'une CL definie sur la normale
            neq = neqs[range(i_node*DIM,(i_node+1)*DIM)];
        updateSpaceValues(nodeeq[i_node]);
        for(unsigned dir=0;dir<DIM;++dir) { /// boucle les directions de l'espace
                expr = read_ex(fcts_spatiales[i_step][dir],space_symbols);
                data[dir] = (TYPEREEL)expr.subs_numerical(space_values);
                data[dir] *= ft[dir];
        }
        /// Assignation
        if(neqs.size() == 0){
            V[range(i_node*DIM,(i_node+1)*DIM)]=data;
        }else{   /// Cas d'une CL definie sur la normale
            Vec<TYPEREEL,DIM> temp=V[range(i_node*DIM,(i_node+1)*DIM)];
            V[range(i_node*DIM,(i_node+1)*DIM)]=ProjT(temp,neq)+ft[0]*data*neq;
        }
    }
    
}