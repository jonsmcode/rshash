#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
#include <filesystem>

int count_unique_x_values(const std::filesystem::path &filepath) {
    // Open the file
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filepath << std::endl;
        return -1; // Return an error code if the file cannot be opened
    }

    int total_unique_count = 0; // Total count of unique x values across all lines
    std::string line;

    // Process each line in the file
    while (std::getline(file, line)) {
        std::unordered_set<int> unique_x_values; // Set to store unique x values for the current line
        std::istringstream line_stream(line);
        std::string tuple;

        // Extract tuples from the line
        while (line_stream >> tuple) {
            // Remove parentheses
            if (tuple.front() == '(' && tuple.back() == ')') {
                tuple = tuple.substr(1, tuple.size() - 2); // Remove '(' and ')'
            }

            // Split the tuple into x and y
            size_t comma_pos = tuple.find(',');
            if (comma_pos != std::string::npos) {
                int x = std::stoi(tuple.substr(0, comma_pos)); // Extract x
                unique_x_values.insert(x); // Add x to the set
            }
        }

        // Add the count of unique x values for this line to the total
        total_unique_count += unique_x_values.size();
    }

    // Close the file
    file.close();

    return total_unique_count; // Return the total count of unique x values
}

int main() {
    // Path to the input file
    std::filesystem::path filepath = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.positions";

    // Count the total number of unique x values across all lines
    int total_unique_count = count_unique_x_values(filepath);

    if (total_unique_count == -1) {
        return 1; // Exit with error if the file could not be opened
    }

    // Print the result
    std::cout << "Total count of unique x values across all lines: " << total_unique_count << std::endl;

    return 0;
}