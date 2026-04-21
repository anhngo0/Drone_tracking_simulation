#include <fstream>
#include <Eigen/Dense>
#include "../../include/json_convert/drone_to_json.hpp"


/*================== TRAJECTORY ====================*/
struct segment_coeffs {
	std::array<double, 4> x;
	std::array<double, 4> y;
	std::array<double, 4> z;
};

struct trajectory_json_struct{
	std::vector<Eigen::Vector3d> points;
	std::vector<segment_coeffs> seg_coeffs_list;
};

void trajectory_to_json(
const std::vector<Eigen::Vector3d>& points,
const std::vector<Cubic_Interpolation_Coeff>& sx,
const std::vector<Cubic_Interpolation_Coeff>& sy,
const std::vector<Cubic_Interpolation_Coeff>& sz,
const std::string& filename)
{
	std::ofstream f(filename);
	trajectory_json_struct traj;
	traj.points = points;
	
	const size_t n = std::min({sx.size(), sy.size(), sz.size()});
	traj.seg_coeffs_list.reserve(n);
	
	for (size_t i = 0; i < n; ++i) {
		segment_coeffs seg;
		seg.x = {sx[i].a, sx[i].b, sx[i].c, sx[i].d};
		seg.y = {sy[i].a, sy[i].b, sy[i].c, sy[i].d};
		seg.z = {sz[i].a, sz[i].b, sz[i].c, sz[i].d};
		
		traj.seg_coeffs_list.push_back(seg);
	}
	if (!f.is_open()) {
		throw std::runtime_error("Cannot open file: " + filename);
	}
	
	f << std::fixed << std::setprecision(8);
	f << "{\n";
		
		/* points */
		f << "  \"points\": [\n";
		for (size_t i = 0; i < traj.points.size(); ++i) {
			const auto& p = traj.points[i];
			f << "    [" << p.x() << ", " << p.y() << ", " << p.z() << "]";
			if (i + 1 < traj.points.size()) f << ",";
			f << "\n";
		}
		f << "  ],\n";
		
		/* segments */
		f << "  \"segments\": [\n";
		for (size_t i = 0; i < traj.seg_coeffs_list.size(); ++i) {
			const auto& s = traj.seg_coeffs_list[i];
			
			f << "    {\n";
				f << "      \"x\": [" << s.x[0] << ", " << s.x[1] << ", " << s.x[2] << ", " << s.x[3] << "],\n";
				f << "      \"y\": [" << s.y[0] << ", " << s.y[1] << ", " << s.y[2] << ", " << s.y[3] << "],\n";
				f << "      \"z\": [" << s.z[0] << ", " << s.z[1] << ", " << s.z[2] << ", " << s.z[3] << "]\n";
				f << "    }";
			
			if (i + 1 < traj.seg_coeffs_list.size()) f << ",";
			f << "\n";
		}
		f << "  ]\n";
		
		f << "}\n";
	
}

/*============== TIME-BASED TRAJECTORY ===============*/

void drone_motion_to_json(Time_Based_Traj time_traj, const std::string& filename) {

    std::ofstream file(filename);

    file << "{\n";
    file << "  \"frames\": [\n";

    bool first = true;
    double t = 0.0;

    while (!time_traj.finished()) {
        if (!first) file << ",\n";
        first = false;

        Eigen::Vector3d p = time_traj.position();
        Eigen::Vector3d v = time_traj.velocity();

        file << "    {\n";
        file << "      \"t\": " << t << ",\n";
        file << "      \"pos\": [" << p.x() << ", " << p.y() << ", " << p.z() << "]\n";
        file << "    }";

        time_traj.update();
        t += time_traj.getTimeStep();
    }

    file << "\n  ]\n}\n";
}

/*for plotting estimated traj*/
void export_trajectories_to_json(
    const std::vector<Eigen::Vector3d>& est_traj, 
    const std::string& filename) 
{
    json j;

    auto vector_to_json = [](const std::vector<Eigen::Vector3d>& traj) {
        json j_list = json::array();
        for (const auto& v : traj) {
            j_list.push_back({v.x(), v.y(), v.z()});
        }
        return j_list;
    };

    j["est_pos"] = vector_to_json(est_traj);

    std::ofstream file(filename);
    if (file.is_open()) {
        file << j.dump(4); 
        file.close();
        // std::cout << "Successfully exported trajectories to " << filename << std::endl;
    } else {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
    }
}

