#include "../include/drone/trajectory.hpp"

/*===================================================================================*/
/* Constructor */
/*===================================================================================*/
Trajectory::Trajectory(
const Eigen::Vector3d& q0_,
const Eigen::Vector3d& qn_,
double cylinder_radius_,
int interpolation_point_number_
)
: q0(q0_),
qn(qn_),
cylinder_radius(cylinder_radius_),
interpolation_point_number(interpolation_point_number_)
{
	cylinder_height = (qn - q0).norm();
	move_vector = q0;
	
	Eigen::Vector3d vec = qn - q0;
	
	azimuth = std::atan2((-1)*vec.y(),vec.x());
	elevator = std::atan2(vec.z(),
	std::sqrt(vec.x()*vec.x() + vec.y()*vec.y()));
}


/*===================================================================================*/
/*--------------------------FIND INTERPOLATION POINT --------------------------------*/
/*===================================================================================*/
void Trajectory::build()
{
	int n = interpolation_point_number - 1;      // number of segments
	int inner_pts = n - 1;                       // random points count
	
	auto rng = create_generator(MERSENE_TWISTER_SEED);
	std::vector<double> angles = generate_random_angles(rng, inner_pts);
	
	auto cylinder_pts = random_points_on_cylinder_surface(
	angles, cylinder_height, cylinder_radius
	);
	
	// auto rotated_pts = cylinder_axis_to_trajectory_axis(
	// azimuth, elevator, move_vector, cylinder_pts
	// );
	auto rotated_pts = cylinder_axis_to_trajectory_axis(
	move_vector, cylinder_pts
	);
	
	interpolation_points = cal_orbit_interpolation_points(
	q0, qn, rotated_pts
	);
	
	setPoints(interpolation_points);
}

std::mt19937 Trajectory::create_generator(unsigned int seed = 0){
	return std::mt19937(seed);
}

std::vector<double> Trajectory::generate_random_angles(std::mt19937& rng, size_t angles_count){
	static const double pi = std::acos(-1.0);
	std::uniform_real_distribution<double> dist(0.0, pi);
	
	std::vector<double> angles;
	angles.reserve(angles_count);
	
	for (size_t i = 0; i < angles_count; i++){
		angles.push_back(dist(rng));
	}
	return angles;
}

std::vector<Eigen::Vector3d> Trajectory::random_points_on_cylinder_surface
( 
const std::vector<double>& angles, 
double cylinder_height, 
double cylinder_radius
)
{
	std::vector<Eigen::Vector3d> random_points;
	random_points.reserve(angles.size());
	for(size_t i = 1; i <= angles.size(); i++){
		Eigen::Vector3d r_point(
		/*n+1 interpolation points -> n-1 random points (exclude starting point and ending point) 
		-> n-1 random angles / we have n spline segment*/
		i * cylinder_height / (angles.size() + 1),
		cylinder_radius * std::cos(angles[i]),
		cylinder_radius * std::sin(angles[i])
		);
		random_points.push_back(r_point);
	}
	return random_points;
}

std::vector<Eigen::Vector3d>
Trajectory::cylinder_axis_to_trajectory_axis(
    Eigen::Vector3d& move_vector,
    std::vector<Eigen::Vector3d>& cylinder_r_points
)
{
    Eigen::Vector3d ex(1,0,0);
    Eigen::Vector3d v = (qn - q0).normalized();

    Eigen::Vector3d axis = ex.cross(v);
    double cos_theta = ex.dot(v);
    double sin_theta = axis.norm();

    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();

    if (sin_theta > 1e-8) {
        axis.normalize();
        Eigen::Matrix3d K;
        K <<     0, -axis.z(),  axis.y(),
             axis.z(),      0, -axis.x(),
            -axis.y(), axis.x(),      0;

        R = Eigen::Matrix3d::Identity()
            + sin_theta * K
            + (1 - cos_theta) * K * K;
    }

    std::vector<Eigen::Vector3d> orbit_r_points;
    orbit_r_points.reserve(cylinder_r_points.size());

    for (const auto& p : cylinder_r_points) {
        orbit_r_points.push_back(R * p + move_vector);
    }

    return orbit_r_points;
}


std::vector<Eigen::Vector3d> Trajectory::cal_orbit_interpolation_points(
const Eigen::Vector3d& starting_point, 
const Eigen::Vector3d& ending_point, 
std::vector<Eigen::Vector3d>& orbit_random_points
)
{
	std::vector<Eigen::Vector3d> orbit_interpolation_points;
	orbit_interpolation_points.reserve(orbit_random_points.size() + 2);
	orbit_interpolation_points.push_back(starting_point);
	for(const auto& orbit_random_point : orbit_random_points){
		orbit_interpolation_points.push_back(orbit_random_point);
	}	
	orbit_interpolation_points.push_back(ending_point);
	return orbit_interpolation_points;
}

void Trajectory::compute1D(
const std::vector<Eigen::Vector3d>& points,
std::vector<Cubic_Interpolation_Coeff>& out, 
int dim) {
	
	int n = points.size() - 1;   /* segment count */
	if (points.size() < 2) return;
	
	int eq = 4*n;
	
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(eq,eq);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(eq);
	
	std::vector<double> t = tNodes;
	
	int row = 0; /*Matrix row*/
	
	for (int i = 0; i < n; i++)
	{
		
		/* 1. Boundary S_0(ti) = Pi */ 
		A(row, 4*i) = 1;         // a0
		b[row] = points[i][dim];
		row++;
		
		/* 2. interpolation condition + continous deriviation*/
		// double h = t[i+1] - t[i];
		double h = tNodes[i+1] - tNodes[i];
		// S_i(t_{i+1}) = P_{i+1}
		A(row, 4*i)     = 1.0;
		A(row, 4*i + 1) = h;
		A(row, 4*i + 2) = h*h;
		A(row, 4*i + 3) = h*h*h;
		b[row] = points[i+1][dim];
		row++;
		
		if (i < n - 1)
		{
			// 1st derivative continuity
			A(row, 4*i + 1) = 1.0;
			A(row, 4*i + 2) = 2.0*h;
			A(row, 4*i + 3) = 3.0*h*h;
			A(row, 4*(i+1) + 1) = -1.0;
			row++;
			
			// 2nd derivative continuity
			A(row, 4*i + 2) = 2.0;
			A(row, 4*i + 3) = 6.0*h;
			A(row, 4*(i+1) + 2) = -2.0;
			row++;
		}
	}
	
	/*3. Natural Boundary (c0 = 0)*/
	A(row, 2) = 2;          // c0 = 0
	b[row] = 0; row++;
	
	A(row, 4*(n-1) + 2) = 2.0; 
	A(row, 4*(n-1) + 3) = 6.0*(tNodes[n] - tNodes[n-1]);
	b[row] = 0.0; row++;
	
	// Sanity check: row should equal eq (4*n)
	if (row != eq) {
		std::cerr << "Row count mismatch: row=" << row << " expected eq=" << eq << std::endl;
		return ;
	}
	
	/*Solve Ax = b*/
	Eigen::VectorXd sol = A.fullPivLu().solve(b);
	
	/*assign coefficients to out result*/
	out.resize(n);
	for (int i = 0; i < n; i++)
	{
		out[i] = {
			sol[4*i + 0],   // a
			sol[4*i + 1],   // b
			sol[4*i + 2],   // c
			sol[4*i + 3],   // d
		};
	}
}

/*=======================================================================================*/
/*----------------------- FIND SPLINE COEFFICIENTS FOR EACH SEGMENT ---------------------*/
/*=======================================================================================*/
/*determine which segment t belongs to*/
int Trajectory::findSegment(double t) const {
	if (t <= tNodes.front()) return 0; /*if t < first t_i = t_i(x_i) value => 1st segment*/
	if (t >= tNodes.back()) return tNodes.size() - 2; /*if t > last value => last segment*/
	
	/*find t if it is in range of t[i] and t[i+1]*/
	for (int i = 0; i < (int)tNodes.size() - 1; i++) {
		if (t >= tNodes[i] && t <= tNodes[i+1]) 
		return i;
	}
	
	/*fallback, 5 interpolation points -> 4 segments -> return 3rd segment/last segment*/
	return tNodes.size() - 2;
}

double Trajectory::getCylinderSegmentLength() const{
	return (qn-q0).norm() / (interpolation_point_number - 1);
}

/*for X axis: X(t) = ai​ + bi *dt + ci ​dt^2 + di*dt^3*/
double Trajectory::eval(const Cubic_Interpolation_Coeff& s, double dt) const {
	return s.a + s.b*dt + s.c*dt*dt + s.d*dt*dt*dt;
}

double Trajectory::eval_derivative(const Cubic_Interpolation_Coeff& s, double dt) const{
	return s.b + 2*s.c*dt + 3*s.d*dt*dt;
}

/*Set interpolation points & sompute spline for segments of each axis X/Y/Z */
void Trajectory::setPoints(const std::vector<Eigen::Vector3d>& pts) {
	int n = pts.size() - 1; /*n+1 interpolation points -> n spline segments*/
	
	tNodes.resize(n + 1);
	for (int i = 0; i <= n; i++)
	tNodes[i] = cylinder_height * (double)i / n;  // t_i = x_i = i*(d/n)
	
	compute1D(pts,sx, 0); /*x axis spline x(t)*/
	compute1D(pts,sy, 1); /*y axis spline y(t)*/
	compute1D(pts,sz, 2); /*z axis spline z(t)*/
}

// Evaluate spline 3D at segment t
Eigen::Vector3d Trajectory::get_pos_by_interpolated_param(double t) const {
	/*get spline segment, if t2 < t< t3 then it is the 2nd segment*/
	int seg = findSegment(t);
	
	/*define dt = t - t_i*/
	double dt = t - tNodes[seg];
	
	return Eigen::Vector3d(
	eval(sx[seg], dt), /*X(t) = ai​ + bi *dt + ci ​dt^2 + di*dt^3*/
	eval(sy[seg], dt), /*Y(t) = ai​ + bi *dt + ci ​dt^2 + di*dt^3*/
	eval(sz[seg], dt)  /*Z(t) = ai​ + bi *dt + ci ​dt^2 + di*dt^3*/
	);
}

Eigen::Vector3d Trajectory::get_derivative_value_of_interpolated_traj_at_u(double u) const {
	/*get spline segment, if t2 < t< t3 then it is the 2nd segment*/
	int seg = findSegment(u);
	
	/*define dt = t - t_i*/
	double du = u - tNodes[seg];
	
	return Eigen::Vector3d(
	eval_derivative(sx[seg], du), /*X(t) = ai​ + bi *dt + ci ​dt^2 + di*dt^3*/
	eval_derivative(sy[seg], du), /*Y(t) = ai​ + bi *dt + ci ​dt^2 + di*dt^3*/
	eval_derivative(sz[seg], du)  /*Z(t) = ai​ + bi *dt + ci ​dt^2 + di*dt^3*/
	);
}

