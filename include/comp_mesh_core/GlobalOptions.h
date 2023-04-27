#pragma once

#include <string>

namespace composite_mesher
{

/// @brief A simple struct to hold global options for the module.
struct GlobalOptions
{
    static GlobalOptions& get()
    {
        static GlobalOptions opts;
        return opts;
    }

    std::string output_dir         = ".";
    bool write_intermediate_meshes = false;
    bool write_debug_meshes = false;
    bool check_common_errors       = false;

    std::string curr_status_padding = "";

  private:
    GlobalOptions() = default;
};

}