#include "ScCsvReader.h"
#include <fstream>
#include <sstream>
#include <iostream>


ScCsvReader::ScCsvReader(std::string filename, bool param_in_rows, bool has_headers):
    __filename(filename),
    __param_in_rows(param_in_rows),
    __has_headers(has_headers),
    __base_header("Parameter_")
{
    std::ifstream in_file;
    
    /// Ouverture du fichier
    in_file.open(__filename.c_str());
    if(not in_file.is_open())
    {
        // Gestion erreur ouverture
        std::string error("Error! Unable to open ")
        error += __filename
        std::cerr << error;
        throw std::ifstream::failure(error);
    }
    
    /// Lecture de la 1ere ligne, detection des separateurs
    std::string line;
    std::getline(in_file, line, '\n');
    int nb_cols = 0;
    char value_separator = ';';
    /// Compte le nombre de ; dans la ligne
    for(std::string::iterator character = line.begin(); character != line.end(); character++)
        if(character == value_separator)
            nb_cols++;
    /// Si il ny en a pas, cest que le separateur est une virgule
    if(nb_cols == 0)
    {
        value_separator = ',';
        for(std::string::iterator character = line.begin(); character != line.end(); character++)
            if(character == value_separator)
                nb_cols++;
    }
    /// On analyse chaque ligne
    parseLine(line,value_separator, nb_cols);
    while(not in_file.eof())
    {
        std::getline(in_file, line, '\n');
        parseLine(line,value_separator, nb_cols);
        /// Rajouter la gestion d'erreur
    }

    /// Fermeture du fichier
    in_file.close()
}

void ScCsvReader::parseLine(const std::string& line, char val_sep, unsigned nb_cols) {
    std::istringstream in_stream(line);
    std::string value;
    std::vector<std::string> row;
    row.reserve(nb_cols);
    for(unsigned i = 0; i < nb_cols, i++) {
        std::getline(in_stream,value,val_sep);
        row.push_back(value);
        /// Rajouter la gestion d'erreur
    }
    __values.push_back(row);
}
