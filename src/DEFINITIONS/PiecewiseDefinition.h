#ifndef PIECEWISEDEFINITION
#define PIECEWISEDEFINITION

#include <vector>

class PiecewiseDomain{
public:
    typedef std::vector<int> SubdomainsList;
    typedef std::vector<int>::iterator SubdomainsIterator;
    typedef std::vector<double> PointsList;
    typedef std::vector<double>::iterator PointsIterator;
    
    PiecewiseDomain(unsigned nb_subdomains = 0,unsigned nb_points = 0):
    total_subdomains(nb_subdomains),
    size_subdomain(nb_subdomains),
    subdomain_index(),
    subdomain_changed(true),
    value_changed(true){
    }
    
    virtual ~PiecewiseDomain(){
    }
    
protected:
    virtual generate_lists() = 0;
    
    unsigned total_subdomains;
    SubdomainsList size_subdomain;
    SubdomainsIterator subdomain_index;
    bool subdomain_changed;
    
    unsigned total_points;
    PointsList value_point;
    PointsIterator point_index;
    bool value_changed;
};

#endif  //PIECEWISEDEFINITION