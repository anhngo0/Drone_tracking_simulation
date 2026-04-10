#include "../include/drone/time_based_traj.hpp"

Time_Based_Traj::Time_Based_Traj(const Trajectory& traj) {

    const auto& sx = traj.getSx();
    const auto& sy = traj.getSy();
    const auto& sz = traj.getSz();
    const auto& pts = traj.getInterpolationPoints();

    std::size_t n = sx.size();
    segments_.resize(n);

    /*get coefficients for each segment*/
    for (std::size_t i = 0; i < n; ++i) {
        segments_[i].sx = sx[i];
        segments_[i].sy = sy[i];
        segments_[i].sz = sz[i];
        segments_[i].u_max = traj.getCylinderSegmentLength();

        /*compute trajectory length*/
        segments_[i].length =
            arcLengthSimpson(segments_[i], segments_[i].u_max, PARTIAL_PATH_IN_EACH_SEGMENT);

        total_length_ += segments_[i].length;
    }
}

void Time_Based_Traj::setVelocity(double v) {
    speed_ = v;
}

void Time_Based_Traj::setTimeStep(double dt) {
    dt_ = dt;
}

void Time_Based_Traj::setAcceleration(double a) {
    accel_ = a;
}

void Time_Based_Traj::setMaxSpeed(double v_max) {
    max_speed_ = v_max;
}

/*path = speed * time*/
void Time_Based_Traj::update() {

    speed_ += accel_ * dt_;

    if (speed_ < 0.0) speed_ = 0.0;
    if (max_speed_ > 0.0 && speed_ > max_speed_)
        speed_ = max_speed_;

    s_ += speed_ * dt_;
    /*check if path's length exceeds trajectory'length*/
    if (s_ > total_length_) {
        s_ = total_length_;
        speed_ = 0.0;
        accel_ = 0.0;
    }
}

/*return position based on path's length*/
Eigen::Vector3d Time_Based_Traj::position() const {
    std::size_t seg_id;
    double u;
    findSegmentAndU(s_, seg_id, u);

    const Segment& seg = segments_[seg_id];

    return Eigen::Vector3d(
        seg.sx.d*u*u*u + seg.sx.c*u*u + seg.sx.b*u + seg.sx.a,
        seg.sy.d*u*u*u + seg.sy.c*u*u + seg.sy.b*u + seg.sy.a,
        seg.sz.d*u*u*u + seg.sz.c*u*u + seg.sz.b*u + seg.sz.a
    );
}

/*return velocity on path's length*/
/*(ds)/(dt) = v, (ds)/(du) = ||R'(u)|| => du / dt = v / (||R'(u)||)*/
Eigen::Vector3d Time_Based_Traj::velocity() const {
    std::size_t seg_id;
    double u;
    findSegmentAndU(s_, seg_id, u);

    const Segment& seg = segments_[seg_id];

    Eigen::Vector3d tangent(
        3*seg.sx.d*u*u + 2*seg.sx.c*u + seg.sx.b,
        3*seg.sy.d*u*u + 2*seg.sy.c*u + seg.sy.b,
        3*seg.sz.d*u*u + 2*seg.sz.c*u + seg.sz.b
    );

    double norm = tangent.norm();
    /*check if norm is 0 or not*/
    if (norm < 1e-8) {
        return Eigen::Vector3d::Zero();
    }
    return tangent * (speed_ / norm);

}

/*find arc length using simpson rule : \int_0^{u_max} ||R'(u)||du */
/*args: 
    -> u_max : maximum value of interpolated value  
    -> N     :  param that divides arc length into small path 
    */
double Time_Based_Traj::arcLengthSimpson(
    const Segment& seg,
    double u_max,
    int N
) const {
    if (N % 2 != 0) N++;
    double h = u_max / N;

    auto f = [&](double u) {
        Eigen::Vector3d d(
            3*seg.sx.d*u*u + 2*seg.sx.c*u + seg.sx.b,
            3*seg.sy.d*u*u + 2*seg.sy.c*u + seg.sy.b,
            3*seg.sz.d*u*u + 2*seg.sz.c*u + seg.sz.b
        );
        return d.norm();
    };

    double sum = f(0.0) + f(u_max);
    for (int i = 1; i < N; ++i) {
        double u = i * h;
        sum += (i % 2 ? 4.0 : 2.0) * f(u);
    }
    return (h / 3.0) * sum;
}

/*using bisection to find interpolated value
when comparing path'length using simpson, 
assign segment number to 50*/
double Time_Based_Traj::findUfromArc(
    const Segment& seg,
    double target_s
) const {
    double lo = 0.0, hi = seg.u_max;
    for (int i = 0; i < 30; ++i) {
        double mid = 0.5 * (lo + hi);
        if (arcLengthSimpson(seg, mid, 300) < target_s)
            lo = mid;
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

void Time_Based_Traj::findSegmentAndU(
    double s,
    std::size_t& seg_id,
    double& u) const {

    double acc = 0.0;

    for (std::size_t i = 0; i < segments_.size(); ++i) {
        if (acc + segments_[i].length >= s) {
            seg_id = i;
            u = findUfromArc(
                segments_[i],
                s - acc);
            return;
        }
        acc += segments_[i].length;
    }

    seg_id = segments_.size() - 1;
    u = segments_.back().u_max;
}

/*reset path*/
void Time_Based_Traj::reset() {
    s_ = 0.0;
}

double Time_Based_Traj::getTotalLength() const{
    return total_length_;
}

bool Time_Based_Traj::finished() const {
    return s_ >= total_length_;
}

TrajectoryData load_trajectory_from_json(const std::string& filename)
{
    std::ifstream ifs(filename);
    if (!ifs) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    json J;
    ifs >> J;

    TrajectoryData traj;

    for (const auto& f : J.at("frames")) {
        FrameData frame;
        frame.t = f.at("t").get<double>();

        auto p = f.at("pos");
        frame.pos = Eigen::Vector3d(
            p[0].get<double>(),
            p[1].get<double>(),
            p[2].get<double>()
        );

        traj.frames.push_back(frame);
    }

    return traj;
}


std::vector<FrameData> load_frames_from_json(const std::string& filename)
{
    std::ifstream ifs(filename);
    if (!ifs) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    json J;
    ifs >> J;

    std::vector<FrameData> framesData;

    for (const auto& f : J.at("frames")) {
        FrameData frame;
        frame.t = f.at("t").get<double>();

        auto p = f.at("pos");
        frame.pos = Eigen::Vector3d(
            p[0].get<double>(),
            p[1].get<double>(),
            p[2].get<double>()
        );

        framesData.push_back(frame);
    }

    return framesData;
}