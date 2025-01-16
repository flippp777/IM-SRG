#include <vector>
#include <string>

struct MatrixElement {
    int s1, s2, ss1, ss2;
    double me;
};

class MatrixElements {
private:
    std::vector<int> states;
    std::vector<MatrixElement> elements;
public:
    double E;
    std::vector<std::vector<double>> f;
    std::vector<std::vector<std::vector<std::vector<double>>>> Gamma;

    MatrixElements(std::string filename);
    const std::vector<int>& getStates() const { return states; }
    const std::vector<MatrixElement>& getElements() const { return elements; }
    std::vector<MatrixElement> truncateElements(int Nmax) const;
    void setTruncatedElements(int Nmax);
    void E0();
    void f0();
    void Gamma0();

};

