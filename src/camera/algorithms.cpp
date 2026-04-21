#include "../../include/camera/algorithms.hpp"


double normal(double mean, double stddev) {
    static thread_local std::mt19937 gen(std::random_device{}());
    std::normal_distribution<double> dist(mean, stddev);
    return dist(gen);
}

inline double deg2rad(double deg) { return deg * M_PI / 180.0; }

inline double clamp(double v, double vmin, double vmax) {
    return std::max(vmin, std::min(vmax, v));
}

void normalize_angle(double& angle) {
    while (angle >  M_PI) angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
}

/*===============================================================*/
/*=================  UNSCENTED KALMAN FILTER =====================*/
/*===============================================================*/

UnscentedKalmanFilter::UnscentedKalmanFilter(double dt, double sigma_a, double sigma_meas) {
    n_x = 9;
    n_sig = 2 * n_x + 1;
    this->dt = dt;

    x = Eigen::VectorXd::Zero(n_x);
  
    x(1) = x(4) = x(7) = 10.0;
    x(3) = 600.0;             

    /*initialize P high to believe in measurement*/
    P = Eigen::MatrixXd::Identity(n_x, n_x) * 1000.0;

    /*initialize state matrix (Constant Acceleration)*/
    F = Eigen::MatrixXd::Identity(9, 9);
    double dt2_ = 0.5 * dt * dt;
    for (int i = 0; i < 3; ++i) { 
        F(i*3, i*3 + 1) = dt;
        F(i*3, i*3 + 2) = dt2_;
        F(i*3 + 1, i*3 + 2) = dt;
    }

    Q = Eigen::MatrixXd::Zero(n_x, n_x);
    Eigen::MatrixXd Q_block(3, 3);
    double dt2 = dt * dt, dt3 = dt2 * dt, dt4 = dt3 * dt;
    
    Q_block << dt4/4.0, dt3/2.0, dt2/2.0,
               dt3/2.0, dt2/1.0, dt,
               dt2/2.0, dt,      1.0;
    Q_block *= (sigma_a * sigma_a);
    
    for(int i=0; i<3; ++i) Q.block<3,3>(i*3, i*3) = Q_block;


    R_base = Eigen::MatrixXd::Zero(2, 2);
    R_base(0, 0) = std::pow(sigma_meas, 2); // Elevation noise
    R_base(1, 1) = std::pow(sigma_meas, 2); // Azimuth noise

    // UKF
    alpha = 1e-3; 
    beta = 2.0;
    kappa = 0.0;
    lambda = alpha * alpha * (n_x + kappa) - n_x;

    weights_m.resize(n_sig);
    weights_c.resize(n_sig);
    weights_m[0] = lambda / (n_x + lambda);
    weights_c[0] = weights_m[0] + (1 - alpha * alpha + beta);
    for (int i = 1; i < n_sig; i++) {
        weights_m[i] = weights_c[i] = 1.0 / (2.0 * (n_x + lambda));
    }
}

void UnscentedKalmanFilter::predict() {

    Eigen::MatrixXd X_sig = Eigen::MatrixXd(n_x, n_sig);
    Eigen::LLT<Eigen::MatrixXd> lltOfP(P);
    
    /*Comput Cholesky matrix of Covariance P for computing Sigma points later*/
    if (lltOfP.info() == Eigen::NumericalIssue) {
        std::cerr << "!!! [UKF Warning] P is not positive definite at frame! ..." << std::endl;
        P = (P + P.transpose()).eval() * 0.5;
        P += Eigen::MatrixXd::Identity(n_x, n_x) * 1e-6;
        lltOfP.compute(P);
    }

    Eigen::MatrixXd A = lltOfP.matrixL();
    
    if (A.array().isNaN().any()) {
        std::cerr << "!!! [UKF Error] Cholesky decomposition produced NaN! ..." << std::endl;
        P = Eigen::MatrixXd::Identity(n_x, n_x) * 500.0; 
        A = P.llt().matrixL();
    }

    double shift = std::sqrt(n_x + lambda);

    /*Compute sigma points*/
    X_sig.col(0) = x;
    for (int i = 0; i < n_x; i++) {
        X_sig.col(i + 1)       = x + shift * A.col(i);
        X_sig.col(i + 1 + n_x) = x - shift * A.col(i);
    }

    X_sig_pred.resize(n_x, n_sig);
    double dt2_ = 0.5 * dt * dt;

    /*Propagate Sigma points though motion model*/
    for (int i = 0; i < n_sig; i++) {
        for (int a = 0; a < 3; a++) { 
            int idx = a * 3;
            double p = X_sig(idx, i), v = X_sig(idx+1, i), acc = X_sig(idx+2, i);
            X_sig_pred(idx, i)   = p + v * dt + acc * dt2_; // Position
            X_sig_pred(idx+1, i) = v + acc * dt;            // Velocity
            X_sig_pred(idx+2, i) = acc;                     // Acceleration
        }
    }

    /*Recalculate mean x*/
    x.setZero();
    for (int i = 0; i < n_sig; i++) x += weights_m[i] * X_sig_pred.col(i);

    /*Recalculate covariance P*/
    P.setZero();
    for (int i = 0; i < n_sig; i++) {
        Eigen::VectorXd diff = X_sig_pred.col(i) - x;
        P += weights_c[i] * diff * diff.transpose();
    }
    P += Q;
}

/*return innovation for debugging*/
void UnscentedKalmanFilter::update(std::vector<Camera*>& cameras)
{
    int num_cams = cameras.size();

    int n_z = 2 * num_cams;

    Eigen::MatrixXd Z_sig(n_z, n_sig);

    for (int i = 0; i < n_sig; i++)
    {
        for (int c = 0; c < num_cams; c++)
        {
            Eigen::Vector3d cam_pos = cameras[c]->getPosition();

            double dx = X_sig_pred(0, i) - cam_pos[0];
            double dy = X_sig_pred(3, i) - cam_pos[1];
            double dz = X_sig_pred(6, i) - cam_pos[2];

            double r = std::sqrt(dx * dx + dz * dz + 1e-6);

            Z_sig(2 * c, i)     = std::atan2(dy, r);   // elevation
            Z_sig(2 * c + 1, i) = std::atan2(dx, dz);  // azimuth
        }
    }

    Eigen::VectorXd z_pred = Eigen::VectorXd::Zero(n_z);

    // for (int i = 0; i < n_sig; i++)
    //     z_pred += weights_m[i] * Z_sig.col(i);
    for (int c = 0; c < num_cams; c++)
    {
        double sum_el = 0.0;

        double sum_sin = 0.0;
        double sum_cos = 0.0;

        for (int i = 0; i < n_sig; i++)
        {
            double el = Z_sig(2*c, i);
            double az = Z_sig(2*c+1, i);

            sum_el += weights_m[i] * el;

            sum_sin += weights_m[i] * std::sin(az);
            sum_cos += weights_m[i] * std::cos(az);
        }

        z_pred(2*c)   = sum_el;
        z_pred(2*c+1) = std::atan2(sum_sin, sum_cos);
    }

    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n_z, n_z);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n_x, n_z);

    for (int i = 0; i < n_sig; i++)
    {
        Eigen::VectorXd z_diff = Z_sig.col(i) - z_pred;

        for (int c = 0; c < num_cams; ++c)
            normalize_angle(z_diff(2 * c + 1));

        Eigen::VectorXd x_diff = X_sig_pred.col(i) - x;

        S += weights_c[i] * z_diff * z_diff.transpose();
        T += weights_c[i] * x_diff * z_diff.transpose();
    }

    /* Add measurement noise */
    for (int c = 0; c < num_cams; c++)
        S.block<2, 2>(2 * c, 2 * c) += R_base;

    Eigen::VectorXd z_meas(n_z);

    for (int c = 0; c < num_cams; c++)
    {
        CameraAngles ang = cameras[c]->get_total_object_to_cam_angles();
        
        /* ADD MEASUREMENTS NOISE HERE*/
        ang.azimuth += normal(0.0, MEASUREMENT_DEVIATION);
        ang.elevation += normal(0.0, MEASUREMENT_DEVIATION);   

        z_meas(2 * c)     = ang.elevation;
        z_meas(2 * c + 1) = ang.azimuth;
    }

    Eigen::VectorXd y = z_meas - z_pred;
    std::vector<double> innovations;
    std::cout << "innovations of visible cameras are ";
    for (int c = 0; c < num_cams; c++)
    { 
        double inno = (z_meas(2 * c) - z_pred(2*c)) * (z_meas(2 * c) - z_pred(2*c)) + 
                        (z_meas(2*c+1) - z_pred(2*c+1)) * (z_meas(2*c+1) - z_pred(2*c+1));
        innovations.push_back(inno);
        std::cout << inno << " || ";
    }
    std::cout << "\n";
    // std::cout << "innovation of this frame is " << innovation_pos.norm() << std::endl;

    for (int c = 0; c < num_cams; ++c)
        normalize_angle(y(2 * c + 1));

    Eigen::MatrixXd K = T * S.inverse();
    x = x + K * y;
    P = P - K * S * K.transpose();
}

/*Same as predicting stage, except we don't save it*/
GaussState UnscentedKalmanFilter::getLeadState(double t, double sigma_a) {
    GaussState result;
    
    Eigen::MatrixXd X_sig = Eigen::MatrixXd(n_x, n_sig);
    Eigen::MatrixXd A = P.llt().matrixL();
    double shift = std::sqrt(n_x + lambda);

    X_sig.col(0) = x;
    for (int i = 0; i < n_x; i++) {
        X_sig.col(i + 1)       = x + shift * A.col(i);
        X_sig.col(i + 1 + n_x) = x - shift * A.col(i);
    }

    Eigen::MatrixXd X_sig_lead = Eigen::MatrixXd(n_x, n_sig);
    double t2 = 0.5 * t * t;

    for (int i = 0; i < n_sig; i++) {
        for (int axis = 0; axis < 3; axis++) {
            int idx = axis * 3;
            double p = X_sig(idx, i);
            double v = X_sig(idx + 1, i);
            double a = X_sig(idx + 2, i);

            X_sig_lead(idx, i)     = p + v * t + a * t2; 
            X_sig_lead(idx + 1, i) = v + a * t;          
            X_sig_lead(idx + 2, i) = a;                 
        }
    }

    result.mean.setZero(n_x);
    for (int i = 0; i < n_sig; i++) {
        result.mean += weights_m[i] * X_sig_lead.col(i);
    }

    result.covariance.setZero(n_x, n_x);
    for (int i = 0; i < n_sig; i++) {
        Eigen::VectorXd diff = X_sig_lead.col(i) - result.mean;
        result.covariance += weights_c[i] * diff * diff.transpose();
    }

    result.covariance += Q; 

    return result;
}

/*========================================================*/
/*=================  PARTICLE FILTER =====================*/
/*========================================================*/

void ParticleFilter::init_particles_from_ukf(
    const Eigen::VectorXd& ukf_x,
    const Eigen::MatrixXd& ukf_P)
{
    particles.clear();

    const int state_dim = ukf_x.size();

    // Cholesky decomposition: P = L * L^T
    Eigen::LLT<Eigen::MatrixXd> llt(ukf_P);
    Eigen::MatrixXd L = llt.matrixL();

    static std::mt19937 gen(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 1.0);

    for (int i = 0; i < num_particles; ++i) {

        // Sample standard normal vector
        Eigen::VectorXd z(state_dim);
        for (int j = 0; j < state_dim; ++j)
            z(j) = dist(gen);

        // Transform to target Gaussian
        Eigen::VectorXd sample = ukf_x + L * z;

        Particle p;
        p.x = sample;
        p.weight = 1.0 / num_particles;

        particles.push_back(p);
    }
}

/*propagate particles through modtion model*/
void ParticleFilter::predict(const Eigen::MatrixXd& F, const Eigen::MatrixXd& Q) {
    Eigen::LLT<Eigen::MatrixXd> lltOfQ(Q);
    Eigen::MatrixXd L = lltOfQ.matrixL();
    std::normal_distribution<double> dist(0.0, 1.0);

    for (auto& p : particles) {
        Eigen::VectorXd process_noise(9);
        for (int j = 0; j < 9; ++j) process_noise(j) = dist(gen);
        
        p.x = F * p.x + L * process_noise;
    }
}

Eigen::VectorXd ParticleFilter::get_estimated_state() {
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(9);
    for (const auto& p : particles) {
        mean += p.x * p.weight;
    }
    return mean;
}

std::vector<Particle> ParticleFilter::pf_get_lead_particle(
    int lead_number,
    const Eigen::MatrixXd& F
){
    // Tính F^lead
    Eigen::MatrixXd F_power = Eigen::MatrixXd::Identity(F.rows(), F.cols());

    for (int i = 0; i < lead_number; ++i)
        F_power = F * F_power;

    std::vector<Particle> prop_particles;

    for (const auto& p : particles)
    {
        Particle temp_particle;
        temp_particle.x = F_power * p.x;
        temp_particle.weight = p.weight;
        prop_particles.push_back(temp_particle);
    }

    return prop_particles;
}

/**
 * @param measurements: list {az, el} of n camera
 * @param cam_poses: position of n camera
 * @param R: measurement noise (2x2)
 */
void ParticleFilter::resample_if_needed() {
    double n_eff = 0;
    for (const auto& p : particles) n_eff += p.weight * p.weight;
    n_eff = 1.0 / n_eff;

    if (n_eff < num_particles / 2.0) {
        std::vector<Particle> new_particles;
        std::uniform_real_distribution<double> dist(0.0, 1.0 / num_particles);
        double r = dist(gen);
        double c = particles[0].weight;
        int i = 0;

        for (int m = 0; m < num_particles; ++m) {
            double u = r + (double)m / num_particles;
            while (u > c && i < num_particles - 1) {
                i++;
                c += particles[i].weight;
            }
            new_particles.push_back(particles[i]);
            new_particles.back().weight = 1.0 / num_particles;
        }
        particles = new_particles;
    }
}

void ParticleFilter::update_weights_focus(
    std::vector<Camera*>& visible_cams,
    std::vector<Camera*>& hidden_cams)
{
    const double sigma = MEASUREMENT_DEVIATION; 

    std::vector<double> log_weights(particles.size());

    for (size_t i = 0; i < particles.size(); i++)
    {
        auto& p = particles[i];

        Eigen::Vector3d pos(p.x(0), p.x(3), p.x(6));

        double log_w = std::log(p.weight + 1e-12);

        if(visible_cams.size() != 0){
            
            for (auto* cam : visible_cams)
            {
                // ===== 1. FOV constraint (secondary) =====
                bool inside = is_point_inside_polyhedron(
                    cam->get_FOV(), Point(pos.x(), pos.y(), pos.z()));
    
                if(!inside) 
                    // log_w += std::log(1e-3);
                    log_w += std::log(1e-4); /*gain panalty */
                else {

                // ===== 2. ANGLE LIKELIHOOD (MAIN) =====
                    Eigen::Vector3d cam_pos = cam->getPosition();
        
                    Eigen::Vector3d d = pos - cam_pos;
        
                    double dx = d.x();
                    double dy = d.y();
                    double dz = d.z();
        
                    double r = std::sqrt(dx*dx + dz*dz + 1e-6);
        
                    double pred_el = std::atan2(dy, r);
                    double pred_az = std::atan2(dx, dz);
        
                    auto ang = cam->get_total_object_to_cam_angles();
        
                    double d_el = pred_el - ang.elevation;
                    double d_az = pred_az - ang.azimuth;
        
                    normalize_angle(d_az);
        
                    double err2 = d_el*d_el + d_az*d_az;
        
                    double log_likelihood = -0.5 * err2 / (sigma*sigma);
        
                    log_w += log_likelihood;
                }
            }    
        }
        if(hidden_cams.size() != 0){
            // for (auto* cam : hidden_cams)
            // {
            //     bool inside = is_point_inside_polyhedron(
            //         cam->get_FOV(), Point(pos.x(), pos.y(), pos.z()));
    
            //     log_w += inside ? std::log(1e-2) : 0.0;
            // }
            /*only penaltize particle once*/
            bool is_in_any_hidden_fov = false;
            for (auto* cam : hidden_cams) {
                if (is_point_inside_polyhedron(cam->get_FOV(), Point(pos.x(), pos.y(), pos.z()))) {
                    is_in_any_hidden_fov = true;
                    break; 
                }
            }
            if (is_in_any_hidden_fov) log_w += std::log(0.1);
        }

        log_weights[i] = log_w;
    }

    // ===== normalize (log-sum-exp) =====
    double max_log = *std::max_element(log_weights.begin(), log_weights.end());

    double sum = 0;

    for (size_t i = 0; i < particles.size(); i++)
    {
        particles[i].weight = std::exp(log_weights[i] - max_log);
        sum += particles[i].weight;
    }

    for (auto& p : particles)
        p.weight /= (sum + 1e-12);
}

GaussState ParticleFilter::compute_ellipsoid_cov(Eigen::Vector3d target)
{
    GaussState res_state;

    int dim = particles[0].x.size();

    res_state.mean = Eigen::VectorXd::Zero(dim);

    for (auto& p : particles)
        res_state.mean += p.weight * p.x;

    res_state.mean(0) = target(0);
    res_state.mean(3) = target(1);
    res_state.mean(6) = target(2);

    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(dim, dim);

    for (auto& p : particles)
    {
        Eigen::VectorXd diff = p.x - res_state.mean;
        C += p.weight * diff * diff.transpose();
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(C);

    Eigen::MatrixXd V = solver.eigenvectors();
    Eigen::VectorXd D = solver.eigenvalues();

    double inflation = 4.0;

    res_state.covariance = V * (inflation * D.asDiagonal()) * V.transpose();

    // res_state.covariance += Eigen::MatrixXd::Identity(dim, dim) * 1e-6;
    return res_state;
}

/*=================================================================*/
/*================ For rotating camera ====================*/

void rotate_cameras(
    std::vector<Camera*> &cameras, 
    const Eigen::Vector3d target_pos)
{
    for(Camera* camera: cameras){
        Eigen::Vector3d dir = target_pos - camera->getPosition();
        camera->rotate_to_direction_vector(dir);
    }
}

/*============================================*/
/*=======WEIGHTED K MEANS CLUSTERING==========*/
/*============================================*/

void compute_gmm_stats(
    const std::vector<Particle>& particles,
    const std::vector<int>& labels,
    int K, /*K clusters/components*/
    std::vector<GaussianComponent>& gmm)
{
    int N = particles.size();
    int dim = particles[0].x.size();

    /*sum of weights in each cluster/component */
    std::vector<double> weight_sum(K, 0.0);

    // ===== Mean =====
    /*Loop in each particle to compute new gmm centroid position 
    gmm[k].mean = \dfrac{\sum_{i=0}^N w^{(i)}_k \chi{(i)}_k}{\sum_{i=0}^N w^{(i)}_k }    */
    for (int i = 0; i < N; ++i)
    {
        int k = labels[i]; /*get label of particle [i]*/
        double w = particles[i].weight;

        gmm[k].mean += w * particles[i].x; /* \sum_{i=0}^N w^{(i)}_k \chi{(i)}_k*/
        weight_sum[k] += w; /*\sum_{i=0}^N w^{(i)}_k */
    }

    for (int k = 0; k < K; ++k)
    {
        if (weight_sum[k] > 1e-12)
            gmm[k].mean /= weight_sum[k];
    }

    // ===== Covariance =====
    /*gmm[k].cov = \dfrac{\sum_{i=0}^N w^{(i)}_k (\chi{(i)}_k - \muy_k) (\chi{(i)}_k - \muy_k)^T}
                        {\sum_{i=0}^N w^{(i)}_k }    */
    for (int i = 0; i < N; ++i)
    {
        int k = labels[i];
        double w = particles[i].weight;

        Eigen::VectorXd diff = particles[i].x - gmm[k].mean;
        gmm[k].cov += w * diff * diff.transpose(); /*\sum_{i=0}^N w^{(i)}_k (\chi{(i)}_k - \muy_k) (\chi{(i)}_k - \muy_k)^T*/
    }

    for (int k = 0; k < K; ++k)
    {
        if (weight_sum[k] > 1e-12)
            gmm[k].cov /= weight_sum[k];

        gmm[k].cov += 1e-6 * Eigen::MatrixXd::Identity(dim, dim);
        gmm[k].weight = weight_sum[k];
    }

    // normalize weights
    double total = 0;
    for (auto& c : gmm) total += c.weight;

    for (auto& c : gmm)
        c.weight /= (total + 1e-12);
}

std::vector<int> weighted_kmeans(
    const std::vector<Particle>& particles,
    int K = CLUSTERING_NUMBER,
    int max_iters = 5,
    double tol = 1e-4)
{
    int N = particles.size();
    int dim = particles[0].x.size();

    std::vector<Eigen::VectorXd> centroids(K);
    std::vector<Eigen::MatrixXd> covs(K, Eigen::MatrixXd::Identity(dim, dim));

    // ===== Init centroids (top weight particles) =====
    std::vector<int> indices(N);
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(),
              [&](int a, int b) {
                  return particles[a].weight > particles[b].weight;
              });

    for (int k = 0; k < K; ++k)
        centroids[k] = particles[indices[k % N]].x;

    std::vector<int> labels(N, 0);

    for (int iter = 0; iter < max_iters; ++iter)
    {
        // ===== Precompute inverse covariance =====
        std::vector<Eigen::MatrixXd> cov_inv(K);
        for (int k = 0; k < K; ++k)
        {
            covs[k] += 1e-6 * Eigen::MatrixXd::Identity(dim, dim); // regularize
            cov_inv[k] = covs[k].inverse();
        }

        // ===== Assign step (Mahalanobis = (\chi^{(i)}_k - \muy_k)^T * P^{-1} * (\chi^{(i)}_k - \muy_k)) =====
        for (int i = 0; i < N; ++i)
        {
            double best_dist = 1e18;
            int best_k = 0;

            for (int k = 0; k < K; ++k)
            {
                Eigen::VectorXd diff = particles[i].x - centroids[k];
                double dist = diff.transpose() * cov_inv[k] * diff;

                if (dist < best_dist)
                {
                    best_dist = dist;
                    best_k = k;
                }
            }
            labels[i] = best_k;
        }

        // ===== Update centroids =====
        std::vector<Eigen::VectorXd> new_centroids(
            K, Eigen::VectorXd::Zero(dim));
        std::vector<double> weight_sum(K, 0.0);

        for (int i = 0; i < N; ++i)
        {
            int k = labels[i];
            double w = particles[i].weight;

            new_centroids[k] += w * particles[i].x;
            weight_sum[k] += w;
        }

        for (int k = 0; k < K; ++k)
        {
            if (weight_sum[k] > 1e-12)
                new_centroids[k] /= weight_sum[k];
            else
                new_centroids[k] = centroids[k];
        }

        // ===== Update covariance =====
        std::vector<Eigen::MatrixXd> new_covs(
            K, Eigen::MatrixXd::Zero(dim, dim));

        for (int i = 0; i < N; ++i)
        {
            int k = labels[i];
            double w = particles[i].weight;

            Eigen::VectorXd diff = particles[i].x - new_centroids[k];
            new_covs[k] += w * (diff * diff.transpose());
        }

        for (int k = 0; k < K; ++k)
        {
            if (weight_sum[k] > 1e-12)
                new_covs[k] /= weight_sum[k];

            // regularization để tránh singular
            new_covs[k] += 1e-6 * Eigen::MatrixXd::Identity(dim, dim);
        }

        // ===== Check convergence =====
        double max_shift = 0.0;
        for (int k = 0; k < K; ++k)
        {
            double shift = (centroids[k] - new_centroids[k]).norm();
            if (shift > max_shift)
                max_shift = shift;
        }

        centroids = new_centroids;
        covs = new_covs;

        if (max_shift < tol)
            break;
    }

    return labels;
}

std::vector<GaussianComponent> fit_gmm_from_particles(
    const std::vector<Particle>& particles,
    int K)
{
    int dim = particles[0].x.size();

    std::vector<int> labels = weighted_kmeans(particles, K);

    std::vector<GaussianComponent> gmm(K);
    for (int k = 0; k < K; ++k)
    {
        gmm[k].mean = Eigen::VectorXd::Zero(dim);
        gmm[k].cov = Eigen::MatrixXd::Zero(dim, dim);
        gmm[k].weight = 0.0;
    }

    compute_gmm_stats(particles, labels, K, gmm);

    return gmm;
}

std::vector<GaussianComponent> fit_gmm_em(
    const std::vector<Particle>& particles,
    int K,
    int max_iters = 5)
{
    int N = particles.size();
    int dim = particles[0].x.size();

    auto gmm = fit_gmm_from_particles(particles, K);

    Eigen::MatrixXd gamma(N, K);

    for (int iter = 0; iter < max_iters; ++iter)
    {
        // ===== Precompute =====
        // ===== Precompute (LLT + log det) =====
        std::vector<Eigen::LLT<Eigen::MatrixXd>> llt(K);
        std::vector<double> log_det(K);

        for (int k = 0; k < K; ++k)
        {
            llt[k].compute(gmm[k].cov);

            if (llt[k].info() != Eigen::Success)
            {
                // fallback: thêm regularization
                gmm[k].cov += 1e-6 * Eigen::MatrixXd::Identity(dim, dim);
                llt[k].compute(gmm[k].cov);
            }

            const auto& L = llt[k].matrixL();

            double logdet = 0.0;
            for (int d = 0; d < dim; ++d)
                logdet += std::log(L(d, d));

            log_det[k] = 2.0 * logdet;
        }

        // ===== E-step =====
       const double log_2pi = std::log(2.0 * M_PI);

        for (int i = 0; i < N; ++i)
        {
            double max_log = -1e18;

            // ===== Compute log-probabilities =====
            for (int k = 0; k < K; ++k)
            {
                Eigen::VectorXd diff = particles[i].x - gmm[k].mean;

                // Solve Σ^{-1} * diff using LLT
                Eigen::VectorXd sol = llt[k].solve(diff);

                // Mahalanobis distance
                double mahal = diff.dot(sol);

                // Log Gaussian probability
                double log_prob =
                    std::log(gmm[k].weight + 1e-12)
                    - 0.5 * (dim * log_2pi + log_det[k] + mahal);

                gamma(i, k) = log_prob;

                if (log_prob > max_log)
                    max_log = log_prob;
            }

            // ===== Log-Sum-Exp normalization =====
            double sum = 0.0;
            for (int k = 0; k < K; ++k)
            {
                gamma(i, k) = std::exp(gamma(i, k) - max_log);
                sum += gamma(i, k);
            }

            gamma.row(i) /= (sum + 1e-12);
        }

        // ===== M-step =====
        for (int k = 0; k < K; ++k)
        {
            double Nk = 0.0;
            Eigen::VectorXd mean = Eigen::VectorXd::Zero(dim);
            Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(dim, dim);

            for (int i = 0; i < N; ++i)
            {
                double g = gamma(i, k) * particles[i].weight;
                Nk += g;
                mean += g * particles[i].x;
            }

            if (Nk > 1e-12)
                mean /= Nk;

            for (int i = 0; i < N; ++i)
            {
                 double g = gamma(i, k) * particles[i].weight;
                Eigen::VectorXd diff = particles[i].x - mean;
                cov += g * diff * diff.transpose();
            }

            if (Nk > 1e-12)
                cov /= Nk;

            cov += 1e-6 * Eigen::MatrixXd::Identity(dim, dim);

            gmm[k].mean = mean;
            gmm[k].cov = cov;
            gmm[k].weight = Nk;
        }

        // normalize weights
        double total = 0;
        for (auto& c : gmm) total += c.weight;

        for (auto& c : gmm)
            c.weight /= (total + 1e-12);
    }

    return gmm;
}

void ParticleFilter::resample_from_gmm(
    const std::vector<GaussianComponent>& gmm)
{
    std::vector<Particle> new_particles;
    new_particles.reserve(num_particles);

    int dim = gmm[0].mean.size();

    // ===== random generator =====
    static std::mt19937 gen(std::random_device{}());
    std::normal_distribution<double> normal_dist(0.0, 1.0);

    // ===== tạo distribution chọn component =====
    std::vector<double> weights;
    for (const auto& comp : gmm)
        weights.push_back(comp.weight);

    std::discrete_distribution<int> comp_dist(weights.begin(), weights.end());

    // ===== precompute Cholesky =====
    std::vector<Eigen::MatrixXd> Ls(gmm.size());
    for (size_t k = 0; k < gmm.size(); ++k)
    {
        Eigen::LLT<Eigen::MatrixXd> llt(gmm[k].cov);
        Ls[k] = llt.matrixL();
    }

    // ===== sampling =====
    for (int i = 0; i < num_particles; ++i)
    {
        int k = comp_dist(gen);

        // sample z ~ N(0, I)
        Eigen::VectorXd z(dim);
        for (int d = 0; d < dim; ++d)
            z(d) = normal_dist(gen);

        // transform
        Eigen::VectorXd sample = gmm[k].mean + Ls[k] * z;

        Particle p;
        p.x = sample;
        p.weight = 1.0 / num_particles;

        new_particles.push_back(p);
    }

    particles = new_particles;
}