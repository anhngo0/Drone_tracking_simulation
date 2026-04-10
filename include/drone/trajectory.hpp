#pragma once 
#include <iostream>
#include <random>
#include <vector>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>


/*=============================================================================*/
/*The spline describes the way object moves from created random point i to {i+1} 
* on the cylinder's surface is as following:
*
*        p_i(t) = a_i + b_i*(t-x_i) + c_i*(t-x_i)^2 + d_i*(t-x_i)^3 [1]
*
* where: 
*    -> a_i, b_i, c_i, d_i are coefficents we need to determine.
-> x_i (= i * d / n)  is the x coordinate of created random point p_i.
-> p_i(t) is the x/y/z of trajectory point at time t.
-> t is a variable that changes along the length d of the pillar axis, and
it can be used as a descriptive parameter to interpolate trajectory.

* when: t = t_i = x_i then p_i = a_i.
* 	 t = t_{i+1} then p_i(t_{i+1}) = a_i + b_i*(t_{i+1}-t_i) + c_i*(t_{i+1}-t_i)^2 + d_i*(t_{i+1}-t_i)   
*
* Boundary Conditions:
* 	  p_i(t_i) = p_i          [2] (partition's starting point)
* 	  p_i(t_{i+1}) = p_{i+1}      (partition's ending point)
* First Deriative Condition and Sencond Deriative Condition
* 	  p'_{i-1}(t_i) = p'_i(t_i) <=> b_{i-1} + 2*c_{i-1}(t_i-t_{i-1}) + 3*d_{i-1}(t_i-t_{i-1})^2 = b_i;
*	  p''_{i-1}(t_i) = p''_i(t_i) <=> 2*c_{i-1} + 6_{i-1}(t_i-t_{i-1}) = 2c_i
* 
*	  p'_i(t_{i+1}) = p'_{i+1}(t_{i+1}) <=> b_i + 2*c_i(t_{i+1}-t_i) + 3*d_i(t_{i+1}-t_i)^2 = b_{i+1};
*	  p''_i{t_{i+1}} = p''_{i+1}(t_{i+1}) <=> 2*c_i + 6*d_i*(t_{i+1} - t_i) = 2*c_{i+1}	

* We use these conditions to form an matrix equation A^{-1}*b = D  
*/
const int MERSENE_TWISTER_SEED = 5;
struct Cubic_Interpolation_Coeff {
	double a, b, c, d; // hệ số spline 1 chiều
};

class Trajectory {
	
	/*input: need to set q0,qn, cylinder_radius, interpolation_point_number*/
	public: 
	Trajectory(
	const Eigen::Vector3d& q0,
	const Eigen::Vector3d& qn,
	double cylinder_radius,
	int interpolation_point_number
	);
	
	/*Function: set interpolation points.
	Result: std::vector<Eigen::Vector3d> interpolation_points  = q0 + (n-1) random points + qn*/
	void build(); 
    
	Eigen::Vector3d get_pos_by_interpolated_param(double u) const;
    Eigen::Vector3d get_derivative_value_of_interpolated_traj_at_u(double u) const; /*u is the interpolated param*/

	double getCylinderSegmentLength() const;
	
    double getLength() const { return cylinder_height; }
	
	/*Getter functions for export file json (drawing by python)*/
	const std::vector<Eigen::Vector3d>& getInterpolationPoints() const { return interpolation_points; }
	const std::vector<Cubic_Interpolation_Coeff>& getSx() const { return sx; }
	const std::vector<Cubic_Interpolation_Coeff>& getSy() const { return sy; }
	const std::vector<Cubic_Interpolation_Coeff>& getSz() const { return sz; }
	private: 
	/*inputs*/
	const Eigen::Vector3d& q0;
	const Eigen::Vector3d& qn;
	double cylinder_radius;
	int interpolation_point_number; /*n+1*/
	
	/*internals*/
	double cylinder_height; /*or distance between starting point and ending point*/ 
	int n = interpolation_point_number - 1; /*segment count*/
	
	/*--------------------------FIND INTERPOLATION POINT --------------------------------*/
	/*> Initially, the cylinder has the center of upper round is at the origin of 
	* euler coordinate and its axis is on the x-axis of euler coordinate */
	
	/*Define the azimuth and elevator angle between cylinder axis 
	* and vector created by starting point and ending point*/
	/*the cylinder axis is (d,0,0) where d is cylinder height*/
	double azimuth, elevator;
	
	/*calculate move vector to move the center of top circular base of cylinder to starting point*/
	Eigen::Vector3d move_vector;
	
	/*interpolation_point list*/
	std::vector<Eigen::Vector3d> interpolation_points;
	
	/*we use mersene twister algoritm to create random angles in 
	* each sample round of cyclinder */
	std::mt19937 create_generator(unsigned int seed);
	std::vector<double> generate_random_angles(std::mt19937& rng, size_t angles_count);
	
	/*The randon angles and radius are used for calculating y and z coordinate of random points 
	on outer surface along the height of the cylinder, while The x-distance between 
	the points is equal, and the sum is equal to the length of the cylinder's height*/
	std::vector<Eigen::Vector3d> random_points_on_cylinder_surface
	( 
	const std::vector<double>& angles, 
	double cylinder_height, 
	double cylinder_radius
	);
	
	/*We change created random points on the initial cylinder postion to trajectory cylinder position:
	+ arg: azimuth and elevator are 2 angles that we use to rotate the cyclinder. After rotating,
	the cylinder's axis has the same direction with trajectory's starting point and ending point vector
	+ arg: move_vector is used to move the center of top circular base coordinate to starting point position*/
	std::vector<Eigen::Vector3d> cylinder_axis_to_trajectory_axis
	(
	Eigen::Vector3d &move_vector,
	std::vector<Eigen::Vector3d>& cylinder_r_points
	);

	/*We create a  interpolation list which has starting point + (n-1) random point + ending point */
	std::vector<Eigen::Vector3d> cal_orbit_interpolation_points(
	const Eigen::Vector3d& starting_point, 
	const Eigen::Vector3d& ending_point, 
	std::vector<Eigen::Vector3d>& orbit_random_points
	);
	
	
	
	/*----------------------- FIND SPLINE COEFFICIENTS FOR EACH SEGMENT ---------------------*/
	std::vector<double> tNodes;    /*t_i at the first segment list: t0,t1,...tn.*/
	std::vector<Cubic_Interpolation_Coeff> sx, sy, sz; /*spline coefficents for each axis x,y,z*/
	
	/*Calculate spline 1D for X; Y; Z*/ 
	/*args: pts -> interpolation points
	out -> list of spline coefficents 
	dim ->  0 for X, 1 for Y, 2 for Z */
	void compute1D(
	const std::vector<Eigen::Vector3d>& pts,
	std::vector<Cubic_Interpolation_Coeff>& s, 
	int dim
	);
	
	/*for X axis: X(t) = ai​ + bi *dt + ci ​dt^2 + di*dt^3*/
	double eval(const Cubic_Interpolation_Coeff& s, double dt) const;

    /*for X axis: X(t) =  bi + 2*ci ​dt + 3*di*dt^2*/
	double eval_derivative(const Cubic_Interpolation_Coeff& s, double dt) const;
	
	/*Set interpolation points & sompute spline for segments of each axis X/Y/Z */
	void setPoints(const std::vector<Eigen::Vector3d>& pts);
	
	/*determine which segment t belongs to*/
	int findSegment(double t) const;
	
};

