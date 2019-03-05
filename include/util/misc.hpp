#pragma once

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <boost/filesystem/path.hpp>

#include "string.hpp"

// Join multiple paths into one path.
template <typename... T>
std::string JoinPaths(T const&... paths);

// Print first-order heading with over- and underscores to `std::cout`.
void PrintHeading1(const std::string& heading);

// Print second-order heading with underscores to `std::cout`.
void PrintHeading2(const std::string& heading);

// Check if vector contains elements.
template <typename T>
bool VectorContainsValue(const std::vector<T>& vector, const T value);

template <typename T>
bool VectorContainsDuplicateValues(const std::vector<T>& vector);

// Parse CSV line to a list of values.
template <typename T>
std::vector<T> CSVToVector(const std::string& csv);

// Concatenate values in list to comma-separated list.
template <typename T>
std::string VectorToCSV(const std::vector<T>& values);

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

template <typename... T>
std::string JoinPaths(T const&... paths) {
    boost::filesystem::path result;
    int unpack[]{0, (result = result / boost::filesystem::path(paths), 0)...};
    static_cast<void>(unpack);
    return result.string();
}

void PrintHeading1(const std::string& heading) {
    std::cout << std::endl << std::string(78, '=') << std::endl;
    std::cout << heading << std::endl;
    std::cout << std::string(78, '=') << std::endl << std::endl;
}

void PrintHeading2(const std::string& heading) {
    std::cout << std::endl << heading << std::endl;
    std::cout << std::string(std::min<int>(heading.size(), 78), '-') << std::endl;
}

template <typename T>
bool VectorContainsValue(const std::vector<T>& vector, const T value) {
    return std::find_if(vector.begin(), vector.end(), [value](const T element) {
        return element == value;
    }) != vector.end();
}

template <typename T>
bool VectorContainsDuplicateValues(const std::vector<T>& vector) {
    std::vector<T> unique_vector = vector;
    return std::unique(unique_vector.begin(), unique_vector.end()) !=
          unique_vector.end();
}

template <>
std::vector<std::string> CSVToVector(const std::string& csv) {
  auto elems = StringSplit(csv, ",;");
  std::vector<std::string> values;
  values.reserve(elems.size());
  for (auto& elem : elems) {
    StringTrim(&elem);
    if (elem.empty()) {
      continue;
    }
      values.push_back(elem);
  }
  return values;
}

template <>
std::vector<int> CSVToVector(const std::string& csv) {
  auto elems = StringSplit(csv, ",;");
  std::vector<int> values;
  values.reserve(elems.size());
  for (auto& elem : elems) {
    StringTrim(&elem);
    if (elem.empty()) {
      continue;
    }
    try {
      values.push_back(std::stoi(elem));
    } catch (std::exception) {
      return std::vector<int>(0);
    }
  }
  return values;
}

template <>
std::vector<float> CSVToVector(const std::string& csv) {
  auto elems = StringSplit(csv, ",;");
  std::vector<float> values;
  values.reserve(elems.size());
  for (auto& elem : elems) {
    StringTrim(&elem);
    if (elem.empty()) {
      continue;
    }
    try {
      values.push_back(std::stod(elem));
    } catch (std::exception) {
      return std::vector<float>(0);
    }
  }
  return values;
}

template <>
std::vector<double> CSVToVector(const std::string& csv) {
  auto elems = StringSplit(csv, ",;");
  std::vector<double> values;
  values.reserve(elems.size());
  for (auto& elem : elems) {
    StringTrim(&elem);
    if (elem.empty()) {
      continue;
    }
    try {
      values.push_back(std::stold(elem));
    } catch (std::exception) {
      return std::vector<double>(0);
    }
  }
  return values;
}

template <typename T>
std::string VectorToCSV(const std::vector<T>& values) {
  std::string string;
  for (const T value : values) {
    string += std::to_string(value) + ", ";
  }
  return string.substr(0, string.length() - 2);
}
