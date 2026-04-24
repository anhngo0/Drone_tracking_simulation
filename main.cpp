#include <iostream>
#include "../include/simulation/simulate.hpp"


int main()
{
    // Simulation program_simulation{};
    // program_simulation.do_simulate();
    // program_simulation.show_result();

    SimulationResult average_result;

    for(int i = 0; i < 100; i++){
        std::cout<< "start simulation loop " << i << std::endl;
        Simulation program_simulation{};
        program_simulation.do_simulate();
        SimulationResult sim_res = program_simulation.getResult();

        average_result.pf_penalty_count += sim_res.pf_penalty_count;
        average_result.number_recover_from_reacquire += sim_res.number_recover_from_reacquire;
        average_result.successful_redetection_percentage += sim_res.successful_redetection_percentage;

        for(const Rmse_based_frame& rmse_frame : sim_res.rmse_list_tracking){

            if(i == 0){
                /*First time, initialize rmse list*/
                average_result.rmse_list_tracking.push_back(rmse_frame);
            }
            else {
                bool is_rmse_frame_existed = false;
                /*Loop through rmse_list and increment count/value of same rmse_frame*/
                for(Rmse_based_frame& avr_res : average_result.rmse_list_tracking){
                    if(avr_res.frame == rmse_frame.frame){
                        avr_res.rmse_count++;
                        avr_res.rmse += rmse_frame.rmse;
                        is_rmse_frame_existed = true;
                        break;
                    }
                }
                /*Add no existed rmse_frame into list*/
                if(!is_rmse_frame_existed)
                    average_result.rmse_list_tracking.push_back(rmse_frame);
            }
        }

        for(const Rmse_based_frame& rmse_frame : sim_res.rmse_list_det_succ){

            if(i == 0){
                /*First time, initialize rmse list*/
                average_result.rmse_list_det_succ.push_back(rmse_frame);
            }
            else {
                bool is_rmse_frame_existed = false;
                /*Loop through rmse_list and increment count/value of same rmse_frame*/
                for(Rmse_based_frame& avr_res : average_result.rmse_list_det_succ){
                    if(avr_res.frame == rmse_frame.frame){
                        avr_res.rmse_count++;
                        avr_res.rmse += rmse_frame.rmse;
                        is_rmse_frame_existed = true;
                        break;
                    }
                }
                /*Add no existed rmse_frame into list*/
                if(!is_rmse_frame_existed)
                    average_result.rmse_list_det_succ.push_back(rmse_frame);
            }
        }

        for(const Rmse_based_frame& rmse_frame : sim_res.rmse_list_det_fail){

            if(i == 0){
                /*First time, initialize rmse list*/
                average_result.rmse_list_det_fail.push_back(rmse_frame);
            }
            else {
                bool is_rmse_frame_existed = false;
                /*Loop through rmse_list and increment count/value of same rmse_frame*/
                for(Rmse_based_frame& avr_res : average_result.rmse_list_det_fail){
                    if(avr_res.frame == rmse_frame.frame){
                        avr_res.rmse_count++;
                        avr_res.rmse += rmse_frame.rmse;
                        is_rmse_frame_existed = true;
                        break;
                    }
                }
                /*Add no existed rmse_frame into list*/
                if(!is_rmse_frame_existed)
                    average_result.rmse_list_det_fail.push_back(rmse_frame);
            }
        }
        program_simulation.show_result();
    }

    // average_result.number_recover_from_reacquire /= average_result.pf_penalty_count; 
    average_result.pf_penalty_count /= 100.0;    
    average_result.number_recover_from_reacquire /= 100.0;
    average_result.successful_redetection_percentage /= 100.0;

    // divide rmse each frame for number it is counted
    for(auto& avr_res : average_result.rmse_list_tracking) 
    {
        avr_res.rmse /= (avr_res.rmse_count*1.0);
    }

     for(auto& avr_res : average_result.rmse_list_det_succ) 
    {
        avr_res.rmse /= (avr_res.rmse_count*1.0);
    }

     for(auto& avr_res : average_result.rmse_list_det_fail) 
    {
        avr_res.rmse /= (avr_res.rmse_count*1.0);
    }

    std::cout << "Average scanning time is " << average_result.pf_penalty_count << std::endl;
    std::cout << "average successful reacquiring  : " 
              << average_result.number_recover_from_reacquire << " ( acounting " 
              << average_result.successful_redetection_percentage  << " \% lost time)\n";
    save_rmse_to_json( average_result.rmse_list_tracking,"rmse_based_frame_tracking_0_035_8.json", LATENCY);
    save_rmse_to_json( average_result.rmse_list_det_succ,"rmse_based_frame_det_succ_0_035_8.json", LATENCY);
    save_rmse_to_json( average_result.rmse_list_det_fail,"rmse_based_frame_det_fail_0_035_8.json", LATENCY);
    return 0;
}
