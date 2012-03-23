#ifndef ALGEBRE_H
#define ALGEBRE_H

using namespace LMT;

template<class TT> int sign(TT a){
    int res;
    double eps=1e-16;
    if((fabs(a)-a)<=eps){
        res=1;
    }
    else {
        res=-1;
    }
    return res;
}

// Recherche valeurs propres matrices 2x2 et 3x3

template<class TT, class STO, class TYP > void eig_jacobi(Mat<TT,STO,TYP> &A,Vec<Vec<TT> > &V, Vec<TT> &D){
    
    typedef Mat<TT,STO,TYP> TTM;
    
    //std::cout << "Attention valable pour un matrice carree symetrique" << std::endl;
    TT eps=1e-8;
    TT pi=3.14159265358979;
    int n=A.nb_rows();
    int nc=A.nb_cols();
    if (nc!=n)
    {std::cout << "matrice non carree" <<std::endl;
    assert(0);}
    TTM B;
    B.resize(n,n);
    B=A;
    
    TT maxi=1.0;
    Vec<int ,2> pq(0,1);
    //int debut=1;
    
    TTM G;
    G.resize(n,n);
    G.set(0.0);
    G.diag()+=1.0;
    
    int nbiter=0;
    
    while (maxi>eps and nbiter<=10000) {
        // recherche de p et q
        maxi=0.0;
        for(unsigned i=0;i<(unsigned)n;++i){
            for(unsigned j=i+1;j<(unsigned)n;++j){
                if (std::abs((TT)A(i,j))>=maxi){
                    pq=Vec<unsigned,2>(i,j);
                    maxi=std::abs((TT)A(i,j));
                }
            }
        }
        
        // calcul de tantheta, cos et sin
        int q=pq[1];
        int p=pq[0];
        TT t;
        TT theta;
        
        if (std::abs((TT)A(p,q))<eps){
            t=0.0;
            theta=0.0;
        }
        else {
            TT v=(A(q,q)-A(p,p))/(2.0*A(p,q));
            
            if (std::abs(v)<eps){
                t=1.0;
                theta=pi/4.0;
            }
            else {
                
                t=-1.0*v + sign(v)*std::sqrt(1+pow(v,2));
                theta=1.0/2.0*(std::atan(1.0/v));
            }
        }
        
        TT c=1.0/std::sqrt(1.0+pow(t,2));
        TT s=c*t;
        c=std::cos(theta);
        s=std::sin(theta);
        
        
        // calcul de la matrice de rotation
        TTM Gnew;
        Gnew.resize(n,n);
        Gnew.set(0.0);
        Gnew.diag()+=1.0;
        
        Gnew(p,p)=c;
        Gnew(p,q)=s;
        Gnew(q,p)=-s;
        Gnew(q,q)=c;
        
        B=trans(Gnew)*A*Gnew;
        Gnew=G*Gnew;
        //assignation de G et A
        A=B;
        G=Gnew;
        nbiter+=1;
        
    }
    
    if (nbiter==10000){
        std::cout << "La methode n'a pas converge au bout de 10000 iterations, l'erreur est de : " <<  maxi<< std::endl;
    }
    
    D.resize(n);
    for(unsigned i=0;i<D.size();++i){
        D[i]=A(i,i);
    }
    
    V.resize(G.nb_rows());
    for(unsigned i=0;i<(unsigned)n;++i){
        V[i].resize(n);
        for(unsigned j=0;j<(unsigned)n;++j){
            V[i][j]=G(j,i); 
        } 
    }
    
};

template<class TT > void orthonormalisation_schmidt(Vec<Vec<TT> > &V){
    
    typedef Vec<TT>  TV;
    unsigned n=V.size();
    for(unsigned i=0;i<n;++i){
        TV vnorm;
        vnorm.resize(n);
        vnorm.set(0.0);
        for(unsigned j=0;j<i;++j){
            vnorm -= dot(V[j],V[i])*V[j];
        }
        V[i]=(V[i]-vnorm);
        V[i]/=std::sqrt(dot(V[i],V[i]));
    }
    
};

template<class TV1> void orthonormalisation_schmidt(TV1 &V){
    
    typedef typename TV1::template SubType<0>::T  TV;
    unsigned n=V.size();
    for(unsigned i=0;i<n;++i){
        TV vnorm;
        vnorm.resize(n);
        vnorm.set(0.0);
        for(unsigned j=0;j<i;++j){
            vnorm -= dot(V[j],V[i])*V[j];
        }
        V[i]=(V[i]-vnorm);
        V[i]/=std::sqrt(dot(V[i],V[i]));
    }
    
};


template<class TM> void affichsparse(TM &M){
    for(unsigned i=0;i<M.nb_rows();++i){
        //std::cout <<M.data[i].indices.size()  << std::endl;
        for(unsigned j=0;j<M.data[i].indices.size();++j){
            std::cout << "("<< i << ","<< M.data[i].indices[j] << ") "<< M.data[i].data[j]<< std::endl;
        }}
};

template<class TM> void affichsparse_row(TM &M, unsigned i){
    for(unsigned j=0;j<M.data[i].indices.size();++j)
        std::cout << "("<< i << ","<< M.data[i].indices[j] << ") "<< M.data[i].data[j]<< std::endl;
};

template<class TV,class TT> TV mini(TV &v, TT &minim){
    TT res=minim;
    if (v[1]<minim[1]){
        res=v;
    }
    else if(v[1]==minim[1]){
        if(v[0]<minim[0])
            res=v;
    }
    return res;
}

template<class TV,class TT> TV maxi(TV &v, TT &maxim){
    TT res=maxim;
    if (v[1]>maxim[1]){
        res=v;
    }
    else if(v[1]==maxim[1]){
        if(v[0]>maxim[0])
            res=v;
    }
    
    return res;
} 

template<class TV, class TT> void getminmax(TV &v, TT &minim, TT &maxim){
    minim=v[0];
    maxim=v[0];
    for(unsigned i=0;i<v.size();++i){
        minim=mini(v[i],minim);
        maxim=maxi(v[i],maxim);
    }
}


// operateur permettant de renvoyer un booleen =1 si p1<p2 d'apres la composante 1 : 
// s'utilise avec sort(v,Less_Vec_col1()); o v est un vecteur de vecteur 
struct Less_Vec_col1{
    template<class P1, class P2>   bool operator()(const P1 &p1,const P2 &p2) const {
        bool res=0;
        if (p1[1]<p2[1]){
            res=1;
        }
        else if(p1[1]==p2[1]){
            if (p1[0]<p2[0])
                res=1;
        }
        return res;
    }
};

// operateur permettant de renvoyer un booleen =1 si p1<p2 d'apres la composante 0 : 
// s'utilise avec sort(v,Less_Vec_col1()); o v est un vecteur de vecteur 
struct Less_Vec_col0{
    template<class P1, class P2>   bool operator()(const P1 &p1,const P2 &p2) const {
        bool res=0;
        if (p1[0]<p2[0]){
            res=1;
        }
        else if(p1[0]==p2[0]){
            if (p1[1]<p2[1])
                res=1;
        }
        return res;
    }
};

template<class T1, class T2, class T3> bool find_with_index(T1 &v1, const T2&v2, T3 &v3){
    bool res=0;
    for(unsigned i=0;i<v1.size();++i){
        if (v1[i]==v2){v3.push_back(i);res=1;}
    }
    return res;
};

//particular matrix
template<class T,class Str,class Sto,class OP,class TT> Mat<T,Str,Sto,OP> ones(TT &m, TT &n){
    Mat<T,Str,Sto,OP> res;
    res.resize(m,n);
    res.diag().set( (T)1 );
    return res;
};

template<class TM,class TT> TM ones(TT &m){
    TM res;
    res.resize(m);
    res.diag().set( 1 );
    return res;
};

template<class T,class TT>  Mat<T, Gen<>, SparseLine<> > spones(TT m, TT n){
    Mat<T, Gen<>, SparseLine<> > res;
    res.resize(m,n);
    for(unsigned i=0;i<m;++i)
        for(unsigned j=0;j<n;++j)
            res(i,j)=(T)1;
        return res;
}


//return 1 if all the vector component are equal to 1 otherwise 0
template<class TV> bool vec_is_all_true(const TV &v){
    bool res=1;
    for(unsigned i=0;i<v.size();++i){
        if(v[i]==0){
            res=0;
            break;
        }
    }
    return res;
}

//calcul de l'inverse d'une matrice seule
template<class T, class IO, class STO> Mat<T,IO,STO> inverse(Mat<T,IO,STO> &A){
    //creation de la matrice inverse :
    Inv<T,IO,STO> invA=inv(A);
    Mat<T,IO,STO> res;
    res.resize(A.nb_rows(),A.nb_cols());
    for(unsigned i=0;i<A.nb_rows();i++){
        Vec<T> B;
        B.resize(A.nb_rows());
        B.set(0);
        B[i]=1;
        Vec<T> X=invA*B;
        res.col(i)=X;
    }
    return res;
}

#endif //ALGEBRE_H
