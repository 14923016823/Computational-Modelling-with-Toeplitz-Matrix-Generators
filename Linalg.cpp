#include <algorithm>
#include <initializer_list>
#include <stdexcept>
#include <vector>
#include <cstdio>

class Vector {
private:
    int length;
    double* vec;

public:
    // Default constructor
    Vector() : length(0), vec(nullptr) {}

    // Explicit size constructor
    explicit Vector(int length)
        : length(length), vec(new double[length])
    {
        for (int i = 0; i < length; ++i)
            vec[i] = 0.0;
    }

    // Constructor from initializer list
    Vector(const std::initializer_list<double>& vecList)
        : length(static_cast<int>(vecList.size())),
          vec(new double[length])  
    {
        std::copy(vecList.begin(), vecList.end(), vec);
    }

    // Destructor 
    ~Vector()
    {
        delete[] vec;
    }

    // Copy constructor: u = Vector(v)
    Vector(const Vector& v)
        : length(v.length), vec(new double[v.length])
    {
        std::copy(v.vec, v.vec + length, vec);
    }

    // Copy assignment: u = v
    Vector& operator=(const Vector& v)
    {
        if (this == &v)
            return *this;

        double* newData = new double[v.length];
        std::copy(v.vec, v.vec + v.length, newData);

        delete[] vec;
        vec    = newData;
        length = v.length;

        return *this;
    }

    // Size
    int size() const { return length; }

    // Element access
    double& operator[](int i) { return vec[i]; }
    double  operator[](int i) const { return vec[i]; }

    void PrintVector() const
    {
        for (int i = 0; i < length; ++i)
            std::printf("vec[%d] = %f\n", i, vec[i]);

        std::printf("done\n");
    }

    // u += v
    Vector& operator+=(const Vector& v)
    {
        if (length != v.length)
            throw std::invalid_argument("Vector lengths mismatch");

        for (int i = 0; i < length; ++i)
            vec[i] += v.vec[i];

        return *this;
    }

    // u + v
    Vector operator+(const Vector& v) const
    {
        Vector result(*this);
        result += v;
        return result;
    }

    // u -= v
    Vector& operator-=(const Vector& v) 
    {
        if (length != v.length)
            throw std::invalid_argument("Vector lengths mismatch");

        for (int i = 0; i < length; ++i)
            vec[i] -= v.vec[i];

        return *this;
    }

    // unary minus: -u
    Vector operator-() const
    {
        Vector result(*this);
        result *= -1.0;
        return result;
    }

    // u - v
    Vector operator-(const Vector& v) const
    {
        Vector result(*this);
        result -= v;
        return result;
    }

    // u *= a
    Vector& operator*=(double a)
    {
        for (int i = 0; i < length; ++i)
            vec[i] *= a;

        return *this;
    }

    // u * a
    Vector operator*(double a) const
    {
        Vector result(*this);
        result *= a;
        return result;
    }

    // a * u
    friend Vector operator*(double a, const Vector& v)
    {
        return v * a;
    }

    // u dot v: u & v
    double operator&(const Vector& v) const
    {
        if (length != v.length)
            throw std::invalid_argument("Vector lengths mismatch");

        double d = 0.0;
        for (int i = 0; i < length; ++i)
            d += vec[i] * v.vec[i];

        return d;
    }

    // u.dot(v)
    double dot(const Vector& v) const { return (*this & v); }
};

class ToeplitzMat {
private:
    int rows;
    int cols;
    int slots;      
    int* diags;     
    double* vals;   

    void sortDiagonals()
    {
        if (slots <= 1) return;

        // permutation indices
        std::vector<int> idx(slots);
        for (int i = 0; i < slots; ++i)
            idx[i] = i;

        // sort indices by diags[]
        std::sort(idx.begin(), idx.end(),
                  [this](int a, int b) {
                      return diags[a] < diags[b];
                  });

        int*    newDiags = new int[slots];
        double* newVals  = new double[slots];

        for (int i = 0; i < slots; ++i) {
            newDiags[i] = diags[idx[i]];
            newVals[i]  = vals[idx[i]];
        }

        delete[] diags;
        delete[] vals;

        diags = newDiags;
        vals  = newVals;
    }

public:
    // Default constructor
    ToeplitzMat()
        : rows(0), cols(0), slots(0),
          diags(nullptr), vals(nullptr)
    {}

    // Constructor from diagonal lists
    ToeplitzMat(int rows,
                int cols,
                std::initializer_list<int> diagList,
                std::initializer_list<double> valsList)
        : rows(rows), cols(cols),
          slots(static_cast<int>(diagList.size())),
          diags(nullptr), vals(nullptr)
    {
        if (rows <= 0 || cols <= 0)
            throw std::invalid_argument("ToeplitzMat: rows and cols must be positive");

        if (diagList.size() != valsList.size())
            throw std::invalid_argument("ToeplitzMat: diag and value list size mismatch");

        if (slots == 0) {
            diags = nullptr;
            vals  = nullptr;
            return;
        }

        diags = new int[slots];
        vals  = new double[slots];

        int i = 0;
        for (int d : diagList) {
            diags[i++] = d;
        }

        i = 0;
        for (double v : valsList) {
            vals[i++] = v;
        }

        // ensure invariant: diagonals are sorted ascending
        sortDiagonals();
    }

    // Destructor
    ~ToeplitzMat()
    {
        delete[] diags;
        delete[] vals;
    }

    // Copy constructor: A = ToeplitzMat(B)
    ToeplitzMat(const ToeplitzMat& B)
        : rows(B.rows),
          cols(B.cols),
          slots(B.slots),
          diags(nullptr),
          vals(nullptr)
    {
        if (slots > 0) {
            diags = new int[slots];
            vals  = new double[slots];
            std::copy(B.diags, B.diags + slots, diags);
            std::copy(B.vals,  B.vals  + slots, vals);
        }
    }

    // Copy assignment: A = B
    ToeplitzMat& operator=(const ToeplitzMat& B)
    {
        if (this == &B)
            return *this;

        int*    newDiags = nullptr;
        double* newVals  = nullptr;

        if (B.slots > 0) {
            newDiags = new int[B.slots];
            newVals  = new double[B.slots];
            std::copy(B.diags, B.diags + B.slots, newDiags);
            std::copy(B.vals,  B.vals  + B.slots, newVals);
        }

        delete[] diags;
        delete[] vals;

        rows  = B.rows;
        cols  = B.cols;
        slots = B.slots;
        diags = newDiags;
        vals  = newVals;

        return *this;
    }

    int numRows() const { return rows; }
    int numCols() const { return cols; }
    int numSlots() const { return slots; }

    void printStructure() const
    {   
        std::printf("Toeplitz matrix structure (%d x %d)\n", rows, cols);
        std::printf("Number of diagonals: %d\n", slots);
        std::printf("---------------------------------------\n");
        std::printf("  diag_index    value\n");
        std::printf("---------------------------------------\n");

        for (int k = 0; k < slots; ++k) {
            std::printf("%10d    %+12.6f\n", diags[k], vals[k]);
        }

        std::printf("---------------------------------------\n");
    }

    void printFull() const
    {
        std::printf("Toeplitz matrix (%d x %d):\n", rows, cols);
        std::printf("---------------------------------------\n");

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {

                int d = j - i;

                auto it = std::lower_bound(diags, diags + slots, d);

                if (it != diags + slots && *it == d) {
                    int k = static_cast<int>(it - diags);
                    std::printf("%8.3f ", vals[k]);
                } else {
                    std::printf("%8.3f ", 0.0);
                }
            }
            std::printf("\n");
        }
        std::printf("---------------------------------------\n");
    }

    // A += B 
    ToeplitzMat& operator+=(const ToeplitzMat& B)
    {
        if (rows != B.rows || cols != B.cols)
            throw std::invalid_argument("ToeplitzMat += : dimension mismatch");

        int i = 0;
        int j = 0;
        int maxNewSlots = slots + B.slots;

        int*    newDiags = new int[maxNewSlots];
        double* newVals  = new double[maxNewSlots];

        int k = 0;

        while (i < slots && j < B.slots)
        {
            if (diags[i] == B.diags[j]) {
                newDiags[k] = diags[i];
                newVals[k]  = vals[i] + B.vals[j];
                ++i; ++j; ++k;
            }
            else if (diags[i] < B.diags[j]) {
                newDiags[k] = diags[i];
                newVals[k]  = vals[i];
                ++i; ++k;
            }
            else { // diags[i] > B.diags[j]
                newDiags[k] = B.diags[j];
                newVals[k]  = B.vals[j];
                ++j; ++k;
            }
        }

        while (i < slots) {
            newDiags[k] = diags[i];
            newVals[k]  = vals[i];
            ++i; ++k;
        }

        while (j < B.slots) {
            newDiags[k] = B.diags[j];
            newVals[k]  = B.vals[j];
            ++j; ++k;
        }

        delete[] diags;
        delete[] vals;

        diags = newDiags;
        vals  = newVals;
        slots = k;

        return *this;
    }

    // A + B
    ToeplitzMat operator+(const ToeplitzMat& B) const
    {
        ToeplitzMat result(*this);
        result += B;
        return result;
    }

    // A *= a
    ToeplitzMat& operator*=(double a)
    {
        for (int i = 0; i < slots; ++i)
            vals[i] *= a;
        return *this;
    }

    // A * a
    ToeplitzMat operator*(double a) const
    {
        ToeplitzMat result(*this);
        result *= a;
        return result;
    }

    // a * A
    friend ToeplitzMat operator*(double a, const ToeplitzMat& A)
    {
        ToeplitzMat result(A);
        result *= a;                 
        return result;
    }

    // A -= B
    ToeplitzMat& operator-=(const ToeplitzMat& B)
    {
        return *this += (-1.0 * B);
    }

    // A - B
    ToeplitzMat operator-(const ToeplitzMat& B) const
    {
        ToeplitzMat result(*this);
        result -= B;
        return result;
    }

    // A & u  (matrix–vector product)
    Vector operator&(const Vector& u) const
    {
        if (u.size() != cols)
            throw std::invalid_argument("ToeplitzMat &: dimension mismatch");

        Vector v(rows);

        for (int i = 0; i < rows; ++i) {
            double sum = 0.0;
            for (int j = 0; j < cols; ++j) {
                int d = j - i;

                int* it = std::lower_bound(diags, diags + slots, d);
                if (it != diags + slots && *it == d) {
                    int k = static_cast<int>(it - diags);
                    sum += vals[k] * u[j];
                }
            }
            v[i] = sum;
        }

        return v;
    }
};

