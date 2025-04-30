#include <cmath>
#include <utility>
#include <stdexcept>
#include <string>
#include <sstream>

#define FOUR_PI 12.566370614359172


/**
 * Represents a four-dimensional vector in spacetime.
 *
 */
struct FourVector {
    double t, x, y, z;

    // Static method to create a unit vector in a specific direction
    static FourVector unit(int direction) {
        switch(direction) {
            case 0: return {1.0, 0.0, 0.0, 0.0};  // time-like unit vector
            case 1: return {0.0, 1.0, 0.0, 0.0};  // x-direction
            case 2: return {0.0, 0.0, 1.0, 0.0};  // y-direction
            case 3: return {0.0, 0.0, 0.0, 1.0};  // z-direction
            default: throw std::out_of_range("Direction must be 0, 1, 2, or 3");
        }
    }

    // Element access via [] operator
    double& operator[](size_t index) {
        switch(index) {
            case 0: return t;
            case 1: return x;
            case 2: return y;
            case 3: return z;
            default: throw std::out_of_range("Index out of range");
        }
    }

    const double& operator[](size_t index) const {
        switch(index) {
            case 0: return t;
            case 1: return x;
            case 2: return y;
            case 3: return z;
            default: throw std::out_of_range("Index out of range");
        }
    }

    // Unary plus (returns the same vector)
    FourVector operator+() const {
        return *this;
    }

    // Unary minus (negates all components)
    FourVector operator-() const {
        return {-t, -x, -y, -z};
    }

    // Vector addition
    FourVector operator+(const FourVector& other) const {
        return {t + other.t, x + other.x, y + other.y, z + other.z};
    }

    // Vector subtraction
    FourVector operator-(const FourVector& other) const {
        return {t - other.t, x - other.x, y - other.y, z - other.z};
    }

    // Scalar multiplication
    FourVector operator*(double scalar) const {
        return {t * scalar, x * scalar, y * scalar, z * scalar};
    }

    // Scalar division
    FourVector operator/(double scalar) const {
        return {t / scalar, x / scalar, y / scalar, z / scalar};
    }

    // Convert to string
    std::string to_string() const {
        return "FourVector(t=" + std::to_string(t) +
               ", x=" + std::to_string(x) +
               ", y=" + std::to_string(y) +
               ", z=" + std::to_string(z) + ")";
    }
};


// Scalar multiplication (scalar * vector)
inline FourVector operator*(double scalar, const FourVector& v) {
    return v * scalar;
}


/**
 * Represents a one-form in spacetime.
 *
 */
struct OneForm {
    double t, x, y, z;

    // Static method to create a unit one-form in a specific direction
    static OneForm unit(int direction) {
        switch(direction) {
            case 0: return {1.0, 0.0, 0.0, 0.0};  // time-like unit one-form
            case 1: return {0.0, 1.0, 0.0, 0.0};  // x-direction
            case 2: return {0.0, 0.0, 1.0, 0.0};  // y-direction
            case 3: return {0.0, 0.0, 0.0, 1.0};  // z-direction
            default: throw std::out_of_range("Direction must be 0, 1, 2, or 3");
        }
    }

    // Element access via [] operator
    double& operator[](size_t index) {
        switch(index) {
            case 0: return t;
            case 1: return x;
            case 2: return y;
            case 3: return z;
            default: throw std::out_of_range("Index out of range");
        }
    }

    const double& operator[](size_t index) const {
        switch(index) {
            case 0: return t;
            case 1: return x;
            case 2: return y;
            case 3: return z;
            default: throw std::out_of_range("Index out of range");
        }
    }

    // Unary plus (returns the same vector)
    OneForm operator+() const {
        return *this;
    }

    // Unary minus (negates all components)
    OneForm operator-() const {
        return {-t, -x, -y, -z};
    }

    // Vector addition
    OneForm operator+(const OneForm& other) const {
        return {t + other.t, x + other.x, y + other.y, z + other.z};
    }

    // Vector subtraction
    OneForm operator-(const OneForm& other) const {
        return {t - other.t, x - other.x, y - other.y, z - other.z};
    }

    // Scalar multiplication
    OneForm operator*(double scalar) const {
        return {t * scalar, x * scalar, y * scalar, z * scalar};
    }

    // Scalar division
    OneForm operator/(double scalar) const {
        return {t / scalar, x / scalar, y / scalar, z / scalar};
    }

    // Convert to string
    std::string to_string() const {
        return "OneForm(t=" + std::to_string(t) +
               ", x=" + std::to_string(x) +
               ", y=" + std::to_string(y) +
               ", z=" + std::to_string(z) + ")";
    }
};


// Scalar multiplication (scalar * one-form)
inline OneForm operator*(double scalar, const OneForm& w) {
    return w * scalar;
}


/**
 * Calculates the square of the spatial norm of a four-vector.
 * This is the Euclidean distance in the spatial components only.
 *
 * @param v The four-vector
 * @return x² + y² + z²
 */
inline double spatial_norm_squared(const FourVector& v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}


/**
 * Calculates the spatial norm of a four-vector.
 * This is the Euclidean distance in the spatial components only.
 *
 * @param v The four-vector
 * @return The spatial norm sqrt(x² + y² + z²)
 */
inline double spatial_norm(const FourVector& v) {
    return std::sqrt(spatial_norm_squared(v));
}


/**
 * Lowers the indices of the FourVector to create a OneForm.
 * Uses Minkowski metric with signature (1, -1, -1, -1) to convert
 * from contravariant to covariant representation.
 *
 * @return OneForm with components transformed by the metric
 */
OneForm lower(const FourVector& u) {
    return OneForm{u.t, -u.x, -u.y, -u.z};
}


/**
 * Raises the indices of the OneForm to create a FourVector.
 * Uses Minkowski metric with signature (1, -1, -1, -1) to convert
 * from covariant to contravariant representation.
 *
 * @return FourVector with components transformed by the metric
 */
FourVector raise(const OneForm& w) {
    return FourVector{w.t, -w.x, -w.y, -w.z};
}


/**
 * Computes the contraction of a one-form with a four-vector.
 * This represents the inner product using the Minkowski metric.
 *
 * @param w The one-form (covariant vector)
 * @param u The four-vector (contravariant vector)
 * @return The scalar result of the contraction
 */
double contract(const OneForm& w, const FourVector& u) {
    return w.t * u.t + w.x * u.x + w.y * u.y + w.z * u.z;
}


/**
 * Calculates the invariant spacetime interval between two four-vectors.
 * Uses the Minkowski metric with signature (1, -1, -1, -1).
 *
 * The interval is given by: ds² = c²(t₂-t₁)² - (x₂-x₁)² - (y₂-y₁)² - (z₂-z₁)²
 * For null-like intervals, ds² = 0
 * For time-like intervals, ds² > 0
 * For space-like intervals, ds² < 0
 *
 * @param p1 The first four-vector point in spacetime
 * @param p2 The second four-vector point in spacetime
 * @return The invariant interval between the two points
 */
inline double invariant_interval(const FourVector& p1, const FourVector& p2) {
    return contract(lower(p2 - p1), p2 - p1);
}


/**
 * Represents a two-form in spacetime (antisymmetric tensor)
 * Used to represent electromagnetic field tensor, etc.
 */
struct TwoForm {
    // Store components as a 4x4 matrix
    double components[4][4];

    // Constructor initializes all components to zero
    TwoForm() {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                components[i][j] = 0.0;
            }
        }
    }

    // Access elements with operator()(mu, nu)
    double& operator()(size_t mu, size_t nu) {
        if (mu >= 4 || nu >= 4) throw std::out_of_range("Index out of range");
        return components[mu][nu];
    }

    const double& operator()(size_t mu, size_t nu) const {
        if (mu >= 4 || nu >= 4) throw std::out_of_range("Index out of range");
        return components[mu][nu];
    }

    // Addition of two-forms
    TwoForm operator+(const TwoForm& other) const {
        TwoForm result;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                result.components[i][j] = components[i][j] + other.components[i][j];
            }
        }
        return result;
    }

    // Subtraction of two-forms
    TwoForm operator-(const TwoForm& other) const {
        TwoForm result;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                result.components[i][j] = components[i][j] - other.components[i][j];
            }
        }
        return result;
    }

    // Scalar multiplication
    TwoForm operator*(double scalar) const {
        TwoForm result;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                result.components[i][j] = components[i][j] * scalar;
            }
        }
        return result;
    }

    // Scalar division
    TwoForm operator/(double scalar) const {
        TwoForm result;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                result.components[i][j] = components[i][j] / scalar;
            }
        }
        return result;
    }

    // Convert to string
    std::string to_string() const {
        std::string result = "TwoForm(\n";
        for (int i = 0; i < 4; i++) {
            result += "  [";
            for (int j = 0; j < 4; j++) {
                result += std::to_string(components[i][j]);
                if (j < 3) result += ", ";
            }
            result += "]";
            if (i < 3) result += ",\n";
        }
        result += "\n)";
        return result;
    }
};


// Scalar multiplication (scalar * two-form)
inline TwoForm operator*(double scalar, const TwoForm& F) {
    return F * scalar;
}


// ============================================================================
// ============================================================================
class ParticleTrajectory {
public:
    /** Angular frequency of oscillation */
    double omega = 0.9;

    /** Amplitude of oscillation */
    double ell = 1.0;

    /**
     * Returns the position four-vector of the particle at time t
     *
     * The particle oscillates along the z-axis with amplitude ell
     * and angular frequency omega.
     *
     * @param t Time at which to calculate the position
     * @return position four-vector of the particle
     */
    FourVector position(double t) const {
        return {t, 0.0, 0.0, ell * std::cos(omega * t)};
    }

    /**
     * Returns the four-velocity of the particle at time t
     *
     * Calculated as the derivative with respect to proper
     * time of the position function.
     *
     * @param t Time at which to calculate the velocity
     * @return four-velocity vector of the particle
     */
    FourVector four_velocity(double t) const {
        double vz = -ell * omega * std::sin(omega * t);
        double u0 = 1.0 / std::sqrt(1.0 - vz * vz);
        return {u0, 0.0, 0.0, u0 * vz};
    }

    // Convert to string
    std::string to_string() const {
        return "ParticleTrajectory(omega=" + std::to_string(omega) +
               ", ell=" + std::to_string(ell) + ")";
    }
};

// Function to find the retarded time using numerical root finding
static double find_retarded_time(const ParticleTrajectory& trajectory, const FourVector& r) {
    const double tolerance = 1e-10;
    const int max_iterations = 100;

    auto condition = [&](double t) {
        return (r.t - t) - spatial_norm(r - trajectory.position(t));
    };

    double t0 = r.t - spatial_norm(r);
    double t1 = r.t - spatial_norm(r) * 2.0;
    double f0 = condition(t0);
    double f1 = condition(t1);

    for (int i = 0; i < max_iterations; i++) {
        double t2 = t1 - f1 * (t1 - t0) / (f1 - f0);
        double f2 = condition(t2);

        if (std::abs(f2) < tolerance || std::abs(t2 - t1) < tolerance) {
            return t2;
        }
        t0 = t1;
        f0 = f1;
        t1 = t2;
        f1 = f2;
    }
    throw std::runtime_error("find_retarded_time: max iterations reached");
}


/**
 * Calculate the Liénard-Wiechert potential for a moving charge
 *
 * The Liénard-Wiechert potentials describe the electromagnetic field of a moving point charge.
 * They account for the retardation effects due to the finite speed of light propagation.
 *
 * The scalar potential Φ and vector potential A are calculated as:
 *   Φ = 1/(4πε₀) * q/[|R|*(1-n·v/c)]
 *   A = (v/c) * Φ
 *
 * where:
 *   - R is the vector from the charge's retarded position to the field point
 *   - n is the unit vector in the direction of R
 *   - v is the charge's velocity at the retarded time
 *
 * @param trajectory The trajectory object describing the particle's motion
 * @param field_point The four-vector where the potentials are evaluated
 * @return A four-vector containing the scalar potential (time component)
 *           and vector potential (spatial parts)
 */
static FourVector lienard_wiechert_potential(const ParticleTrajectory& trajectory, const FourVector& field_point) {
    double t_ret = find_retarded_time(trajectory, field_point);
    FourVector source_position = trajectory.position(t_ret);
    FourVector R = field_point - source_position;
    FourVector u = trajectory.four_velocity(t_ret);

    double vx = u.x / u.t;
    double vy = u.y / u.t;
    double vz = u.z / u.t;

    double R_mag = spatial_norm(R);
    double nx = R.x / R_mag;
    double ny = R.y / R_mag;
    double nz = R.z / R_mag;

    double beta_dot_n = nx * vx + ny * vy + nz * vz;
    double phi = 1.0 / R_mag / (1.0 - beta_dot_n);
    double Ax = vx * phi;
    double Ay = vy * phi;
    double Az = vz * phi;

    return {phi, Ax, Ay, Az};
}


/**
 * Compute the derivative of the four-potential A with respect to spacetime coordinates
 *
 * This function calculates the partial derivative ∂A/∂x[nu] at a given field point
 * using central finite differencing.
 *
 * @param trajectory The trajectory object describing the particle's motion
 * @param field_point The four-vector where the derivative is evaluated
 * @param nu The direction of the derivative
 * @param epsilon The step size for finite differencing (default: 1e-6)
 * @return The derivative ∂A[mu]/∂x[nu]
 */
static FourVector vector_potential_derivative(
    const ParticleTrajectory& trajectory,
    const FourVector& field_point,
    int nu,
    int fd_order,
    double epsilon
) {
    // Implement different finite-difference schemes based on the order
    switch(fd_order) {
        case 2: { // 2nd order central difference O(h²)
            FourVector point_p = field_point + FourVector::unit(nu) * epsilon;
            FourVector point_m = field_point - FourVector::unit(nu) * epsilon;

            FourVector A_p = lienard_wiechert_potential(trajectory, point_p);
            FourVector A_m = lienard_wiechert_potential(trajectory, point_m);

            return (A_p - A_m) / (2.0 * epsilon);
        }
        case 4: { // 4th order central difference O(h⁴)
            FourVector point_p1 = field_point + FourVector::unit(nu) * epsilon;
            FourVector point_p2 = field_point + FourVector::unit(nu) * (2.0 * epsilon);
            FourVector point_m1 = field_point - FourVector::unit(nu) * epsilon;
            FourVector point_m2 = field_point - FourVector::unit(nu) * (2.0 * epsilon);

            FourVector A_p1 = lienard_wiechert_potential(trajectory, point_p1);
            FourVector A_p2 = lienard_wiechert_potential(trajectory, point_p2);
            FourVector A_m1 = lienard_wiechert_potential(trajectory, point_m1);
            FourVector A_m2 = lienard_wiechert_potential(trajectory, point_m2);

            return (A_m2 - 8.0 * A_m1 + 8.0 * A_p1 - A_p2) / (12.0 * epsilon);
        }
        case 6: { // 6th order central difference O(h⁶)
            FourVector point_p1 = field_point + FourVector::unit(nu) * epsilon;
            FourVector point_p2 = field_point + FourVector::unit(nu) * (2.0 * epsilon);
            FourVector point_p3 = field_point + FourVector::unit(nu) * (3.0 * epsilon);
            FourVector point_m1 = field_point - FourVector::unit(nu) * epsilon;
            FourVector point_m2 = field_point - FourVector::unit(nu) * (2.0 * epsilon);
            FourVector point_m3 = field_point - FourVector::unit(nu) * (3.0 * epsilon);

            FourVector A_p1 = lienard_wiechert_potential(trajectory, point_p1);
            FourVector A_p2 = lienard_wiechert_potential(trajectory, point_p2);
            FourVector A_p3 = lienard_wiechert_potential(trajectory, point_p3);
            FourVector A_m1 = lienard_wiechert_potential(trajectory, point_m1);
            FourVector A_m2 = lienard_wiechert_potential(trajectory, point_m2);
            FourVector A_m3 = lienard_wiechert_potential(trajectory, point_m3);

            return (-A_m3 + 9.0 * A_m2 - 45.0 * A_m1 + 45.0 * A_p1 - 9.0 * A_p2 + A_p3) / (60.0 * epsilon);
        }
        default: {
            throw std::invalid_argument("Invalid finite difference order. Supported orders: 2, 4, 6");
        }
    }
}


/**
 * Calculate the Faraday tensor from the four-potential using finite differencing
 *
 * The Faraday tensor Fμν is defined as the curl of the four-potential:
 *   Fμν = ∂μAν - ∂νAμ
 *
 * This function computes Fμν using central difference approximations for the
 * derivatives of the four-potential.
 *
 * @param trajectory The trajectory object describing the particle's motion
 * @param field_point The four-vector where the tensor is evaluated
 * @param epsilon The step size for finite differencing
 * @return A 4×4 array representing the Faraday tensor at the field point
 */
static TwoForm faraday_tensor(
    const ParticleTrajectory& trajectory,
    const FourVector& field_point,
    int fd_order,
    double epsilon
) {
    OneForm dA[4] = {
        lower(vector_potential_derivative(trajectory, field_point, 0, fd_order, epsilon)),
        lower(vector_potential_derivative(trajectory, field_point, 1, fd_order, epsilon)),
        lower(vector_potential_derivative(trajectory, field_point, 2, fd_order, epsilon)),
        lower(vector_potential_derivative(trajectory, field_point, 3, fd_order, epsilon))
    };
    TwoForm F;

    // Fμν = ∂μAν - ∂νAμ
    for (size_t mu = 0; mu < 4; mu++) {
        for (size_t nu = 0; nu < 4; nu++) {
            F(mu, nu) = dA[mu][nu] - dA[nu][mu];
        }
    }
    return F;
}


/**
 * Extracts the electric field components from the Faraday tensor
 *
 * In the Faraday tensor with lowered indices and mostly-minus metric convention:
 *   E_x = +F_{01} = +F_{t,x}
 *   E_y = +F_{02} = +F_{t,y}
 *   E_z = +F_{03} = +F_{t,z}
 *
 * @param F The Faraday tensor with lowered indices
 * @return A three-component vector (as FourVector with t=0) representing the electric field
 */
static FourVector electric_field(const TwoForm& F) {
    return {0.0, +F(0, 1), +F(0, 2), +F(0, 3)};
}


/**
 * Extracts the magnetic field components from the Faraday tensor
 *
 * In the Faraday tensor with lowered indices and mostly-minus metric convention:
 *   B_x = -F_{23} = -F_{y,z}
 *   B_y = -F_{31} = -F_{z,x}
 *   B_z = -F_{12} = -F_{x,y}
 *
 * @param F The Faraday tensor with lowered indices
 * @return A three-component vector (as FourVector with t=0) representing the magnetic field
 */
static FourVector magnetic_field(const TwoForm& F) {
    return {0.0, -F(2, 3), -F(3, 1), -F(1, 2)};
}


/**
 * Calculate the Poynting vector from the electric and magnetic fields
 *
 * The Poynting vector S represents the directional energy flux density of an
 * electromagnetic field. It is defined as:
 *   S = E × B
 * where E is the electric field and B is the magnetic field.
 *
 * @param F The Faraday tensor containing the electromagnetic field
 * @return A "FourVector" which is really the left-most column of the EM energy-momentum tensor
 */
static FourVector poynting_vector(const TwoForm& F) {
    FourVector E = electric_field(F);
    FourVector B = magnetic_field(F);

    return FourVector{
        0.5 * (spatial_norm_squared(E) + spatial_norm_squared(B)),
        E.y * B.z - E.z * B.y,  // (E × B)_x
        E.z * B.x - E.x * B.z,  // (E × B)_y
        E.x * B.y - E.y * B.x   // (E × B)_z
    } / FOUR_PI;
}




#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

namespace py = pybind11;
PYBIND11_MODULE(_core, m) {
    m.doc() = "Core C++ extension functions";
    py::class_<FourVector>(m, "FourVector")
        .def(py::init<>())
        .def(py::init<double, double, double, double>())
        .def_static("unit", &FourVector::unit, "Create a unit vector in a specific direction")
        .def_readwrite("t", &FourVector::t)
        .def_readwrite("x", &FourVector::x)
        .def_readwrite("y", &FourVector::y)
        .def_readwrite("z", &FourVector::z)
        .def("__getitem__", [](const FourVector &v, size_t i) { return v[i]; })
        .def("__setitem__", [](FourVector &v, size_t i, double value) { v[i] = value; })
        .def_property_readonly("spatial_norm", spatial_norm)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(+py::self)
        .def(-py::self)
        .def(py::self * double())
        .def(py::self / double())
        .def(double() * py::self)
        .def("__repr__", &FourVector::to_string)
        .def("__str__", &FourVector::to_string);

    py::class_<OneForm>(m, "OneForm")
        .def(py::init<>())
        .def(py::init<double, double, double, double>())
        .def_static("unit", &OneForm::unit, "Create a unit one-form in a specific direction")
        .def_readwrite("t", &OneForm::t)
        .def_readwrite("x", &OneForm::x)
        .def_readwrite("y", &OneForm::y)
        .def_readwrite("z", &OneForm::z)
        .def("__getitem__", [](const OneForm &w, size_t i) { return w[i]; })
        .def("__setitem__", [](OneForm &w, size_t i, double value) { w[i] = value; })
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(+py::self)
        .def(-py::self)
        .def(py::self * double())
        .def(py::self / double())
        .def(double() * py::self)
        .def("__repr__", &OneForm::to_string)
        .def("__str__", &OneForm::to_string);

    py::class_<TwoForm>(m, "TwoForm")
        .def(py::init<>())
        .def("__getitem__", [](const TwoForm &F, std::pair<size_t, size_t> indices) {
            return F(indices.first, indices.second);
        })
        .def("__setitem__", [](TwoForm &F, std::pair<size_t, size_t> indices, double value) {
            F(indices.first, indices.second) = value;
        })
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(py::self / double())
        .def(double() * py::self)
        .def("__repr__", &TwoForm::to_string)
        .def("__str__", &TwoForm::to_string);

    py::class_<ParticleTrajectory>(m, "ParticleTrajectory")
        .def(py::init<double, double>(), py::arg("omega") = 0.9, py::arg("ell") = 1.0)
        .def_readwrite("omega", &ParticleTrajectory::omega)
        .def_readwrite("ell", &ParticleTrajectory::ell)
        .def("position", &ParticleTrajectory::position)
        .def("four_velocity", &ParticleTrajectory::four_velocity)
        .def("__repr__", &ParticleTrajectory::to_string)
        .def("__str__", &ParticleTrajectory::to_string);

    m.def("spatial_norm_squared", &spatial_norm_squared, "Calculate the square of the spatial norm of a four-vector");
    m.def("spatial_norm", &spatial_norm, "Calculate the spatial norm of a four-vector");
    m.def("invariant_interval", &invariant_interval, "Calculate the invariant interval between two points");
    m.def("lower", &lower, "Lower indices of a four-vector to create a one-form");
    m.def("raise", &raise, "Raise indices of a one-form to create a four-vector");
    m.def("contract", &contract, "Contract a one-form with a four-vector");
    m.def("find_retarded_time", &find_retarded_time, "Find the retarded time for a field point");
    m.def("lienard_wiechert_potential", &lienard_wiechert_potential, "Calculate the Liénard-Wiechert potential");
    m.def("vector_potential_derivative", &vector_potential_derivative, "Calculate the derivative of the four-potential");
    m.def("faraday_tensor", &faraday_tensor, "Calculate the Faraday tensor from the four-potential");
    m.def("electric_field", &electric_field, "Extract the electric field from the Faraday tensor");
    m.def("magnetic_field", &magnetic_field, "Extract the magnetic field from the Faraday tensor");
    m.def("poynting_vector", &poynting_vector, "Calculate the Poynting vector from the electromagnetic field");
}
