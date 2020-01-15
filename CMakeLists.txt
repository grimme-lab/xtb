cmake_minimum_required(VERSION 3.9)

# Setup the XTB Project
project(xtb
        VERSION 6.2.2)
enable_testing()

# Define Options

# Add subdirectories (the actual XTB setup)
add_subdirectory(cmake)
