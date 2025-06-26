#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>

#ifndef IO_HPP
#define IO_HPP

void print_vector(std::vector<double> &vec){
    std::cout<<"Vector"<<std::endl;

    for (size_t i = 0; i < vec.size(); ++i){
        std::cout << vec[i] << std::endl;
    }
}

class CSV1D{

    public:
        CSV1D(SimpleMesh *mesh);
        ~CSV1D();
        void add_field(std::vector<double> *field, const std::string& name);
        void write_csv(const std::string& file) const;
    
    private:
        std::vector<const std::vector<double>*> fields_;
        std::vector<std::string> names_;
};

CSV1D::CSV1D(SimpleMesh* mesh) {
    names_.push_back("x [m]");
    fields_.push_back(&mesh->x());
}

CSV1D::~CSV1D() {
    // Nothing to clean up, since we're not owning memory
}

void CSV1D::add_field(std::vector<double>* field, const std::string& name) {
    if(fields_[0]->size() != (*field).size()) throw std::runtime_error("Sizes don't match"); 
    fields_.push_back(field);
    names_.push_back(name);
}

void CSV1D::write_csv(const std::string& file) const {
    std::ofstream out(file);
    if (!out) {
        std::cerr << "Error opening file for writing: " << file << std::endl;
        return;
    }

    // Write header
    for (size_t i = 0; i < names_.size(); ++i) {
        out << names_[i];
        if (i + 1 < names_.size()) out << ",";
    }
    out << "\n";

    // Assume all fields are the same size as x
    size_t n = fields_[0]->size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < fields_.size(); ++j) {
            out << (*fields_[j])[i];
            if (j + 1 < fields_.size()) out << ",";
        }
        out << "\n";
    }

    out.close();
}


#endif