#ifndef SCCSVREADER
#define SCCSVREADER

#include <string>
#include <vector>


/** ScCsvReader (class)
* \inbrief A simple CSV reader
* ScCsvReader is a simple CSV parser.
* It needs to know :
*     - the path to the file
*     - if the parameters are stored by rows or columns
*     - if headers are provided
* 
* It gives access to :
*     - 
* Errors thrown :
*     - std::ios_base::failure, if cannot open the file
*     - std::runtime_error, if parsing fails
*/
class ScCsvReader
{
public:
    /// Constructor
    /// Parse the file and fill __headers and __values
    /// Throw an error if cannot open the file or parsing fails (see parseLine)
    ScCsvReader(std::string filename, bool param_in_rows, bool has_headers);
    
    /// Special
    void setDefaultBaseHeader(const std::string& base_header) {__base_header = base_header;}

    /// Accessors
    unsigned nbParameters() const {
        if(__param_in_rows)
            return __values.size();
        else
            return __values[0].size();
    }
    unsigned nbValues(unsigned id_param) const {
        if(__param_in_rows)
            return __values[id_param].size() - (__has_headers?1,0);
        else
            return __values[0].size();
    }
    std::string getValues(unsigned id_param, unsigned id_value) const {
        return __values[(__param_in_rows?id_param,id_value+(__has_headers?1,0))][(__param_in_rows?id_value(__has_headers?1,0),id_param)];
    }
    std::string getHeader(unsigned id_param) const {
        if(__has_headers)
           return __values[(__param_in_rows?id_param,0)][(__param_in_rows?0,id_param)];
        else {
           std::stringstream tmp(__base_header);
           tmp << id_param;
           return tmp.str();
        }
    }
    std::string operator()(unsigned id_param, unsigned id_value) const {
        return getValues(id_param, id_value);
    }

protected:
    /// Parse the "line" with "val_sep" as value separator and check if "nb_cols" is respected
    /// Throw an error if fails
    void parseLine(const std::string& line, char val_sep, unsigned nb_cols);
    
    std::string __filename; /// Path to the file
    bool __param_in_rows;   /// Are the parameters stored by row? Otherwise, they are in column
    bool __has_headers;     /// Are there headers? Then first values of parameters will be assumed to be their names
    std::string __base_header; /// Will be used to build headers if __has_headers is false
    std::vector< std::vector<std::string> > __values;    /// Store the values, first index is row's index
}

#endif //SCCSVREADER
