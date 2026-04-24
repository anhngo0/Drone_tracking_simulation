#include <Eigen/Dense>
#include <cmath>
#include <numeric>

#include <CGAL/Side_of_triangle_mesh.h>
#include "camera.hpp"

/*
 * @brief normalize to  [-PI, PI]
 */
void normalize_angle(double& angle);

struct GaussState {
    Eigen::VectorXd mean;
    Eigen::MatrixXd covariance;
};

/*angle (rad) | 0.001 rad \approx 0.0573 deg, when convert to pixel, it approximates 80pix,
 0.01745 rad \approx 1 deg , 1396 pixel */
const double MEASUREMENT_DEVIATION =  0.026175; 
// const double MEASUREMENT_DEVIATION =  0.00873; 

class UnscentedKalmanFilter {
public:
    /**
     * @brief Khởi tạo bộ lọc UKF
     * @param dt timestep (second)
     * @param sigma_a acceleration covariance (m/s^2)
     * @param sigma_meas measurement (angle) covariance (rad)
     */
    UnscentedKalmanFilter(double sigma_a, double sigma_meas);

    void predict(double dt);

    void update(std::vector<Camera*>& cameras);

    Eigen::VectorXd getState() const { return x; }
    Eigen::MatrixXd getCovariance() const { return P; }
    Eigen::MatrixXd getFMatrix() const {return F;}
    Eigen::MatrixXd getRMatrix() const {return R_base;}

    Eigen::MatrixXd getQMatrix() const { return Q; }

    GaussState getLeadState(double t);
   
    void setState(const Eigen::VectorXd& new_x) { x = new_x; }
    void setCovariance (const Eigen::MatrixXd& new_Cov){P = new_Cov;}

private:

    int n_x;   // state dimensions (9)
    int n_sig; //  Sigma points (2*n_x + 1)
    // double dt;
    double sigma_a;

    Eigen::VectorXd x; // state: [px, vx, ax, py, vy, ay, pz, vz, az]
    Eigen::MatrixXd P; 
    Eigen::MatrixXd F;
    Eigen::MatrixXd Q; //(Process Noise)
    Eigen::MatrixXd R_base; // measurement noise (2x2)

    Eigen::MatrixXd X_sig_pred; // Sigma points are predicted though F

    /*UKF parameters*/
    double alpha;
    double beta;
    double kappa;
    double lambda;

    // weights for Mean and Covariance
    Eigen::VectorXd weights_m;
    Eigen::VectorXd weights_c;

};

/*============================================*/
/*=======WEIGHTED K MEANS CLUSTERING==========*/
/*============================================*/

const int CLUSTERING_NUMBER = 4;

struct GaussianComponent {
    Eigen::VectorXd mean;
    Eigen::MatrixXd cov;
    double weight;
};

/*============================================*/
/*============  PARTICLE FILTER ==============*/
/*============================================*/

struct Particle {
    Eigen::VectorXd x; // Trạng thái 9x1
    double weight;

     Particle() : x(Eigen::VectorXd::Zero(9)), weight(1.0) {}
};

class ParticleFilter {
private:
    int num_particles;
    std::vector<Particle> particles;
    std::default_random_engine gen;

public:
    ParticleFilter(int n) : num_particles(n) {
        particles.resize(num_particles);
    }

    void init_particles_from_ukf(
        const Eigen::VectorXd& ukf_x,
        const Eigen::MatrixXd& ukf_P
    );

    void predict(const double dt, const double sigma_a);

    void resample_if_needed();
    void resample_from_gmm(const std::vector<GaussianComponent>& gmm);

    Eigen::VectorXd get_estimated_state();

    std::vector<Particle> pf_get_lead_particle(
    int lead_number,
    const Eigen::MatrixXd& F);

    /*For rotating object hidden camera*/
    std::vector<Particle> get_particles(){return particles;};

    void update_weights_focus(
    std::vector<Camera*>& visible_cams,
    std::vector<Camera*>& hidden_cams);

    Eigen::MatrixXd compute_ellipsoid_cov();
};

/*============================================*/
/*============  PARTICLE FILTER ==============*/
/*============================================*/

std::vector<int> weighted_kmeans(
    const std::vector<Particle>& particles,
    int K,
    int max_iters,
    double tol
);

std::vector<Eigen::VectorXd> weighted_kmeans_return_centers(
    const std::vector<Particle>& particles,
    int K,
    int max_iters ,
    double tol);

std::vector<GaussianComponent> fit_gmm_from_particles(
    const std::vector<Particle>& particles,
    int K
);

std::vector<GaussianComponent> fit_gmm_em(
    const std::vector<Particle>& particles,
    int K,
    int max_iters);