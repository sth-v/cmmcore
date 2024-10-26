//
// Created by Andrew Astakhov on 20.10.24.
//

#ifndef CMMCORE_MEMORYVIEW_H
#define CMMCORE_MEMORYVIEW_H
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <numeric>
#include <iostream>
#include <cassert>
#include <sstream>
#include <string>
#include <typeinfo>
#include <cstring> // For memcpy
#include <cstdio>  // For FILE*
namespace cmmcore
{


    // Define a Slice struct to represent slicing parameters
    struct Slice
    {
        // If single_index is true, the slice represents a single index
        bool single_index;
        // For single index slicing
        int index;

        // For range slicing
        bool has_start;
        int start;
        bool has_stop;
        int stop;
        int step;

        // Constructors for different slicing types
        // Entire range
        Slice() : single_index(false), has_start(false), has_stop(false), step(1)
        {
        }

        // Single index
        Slice(int idx) : single_index(true), index(idx), has_start(false), has_stop(false), step(1)
        {
        }

        // Range slice with start and stop
        Slice(int s, int e, int st = 1) : single_index(false), has_start(true), start(s),
                                          has_stop(true), stop(e), step(st)
        {
        }
    };

    // Helper function to compute the total number of elements from shape
    size_t compute_total_size(const std::vector<size_t>& shape)
    {
        if (shape.empty()) return 0;
        return std::accumulate(shape.begin(), shape.end(), static_cast<size_t>(1), std::multiplies<size_t>());
    }

    template <typename T>
    class MemoryView
    {
    public:
        // Constructors
        MemoryView(T* data, const std::vector<size_t>& shape)
            : data_(data), shape_(shape), offset_(0)
        {
            compute_strides();
        }

        // Internal constructor used for slicing
        MemoryView(T* data, const std::vector<size_t>& shape, const std::vector<size_t>& strides, size_t offset)
            : data_(data), shape_(shape), strides_(strides), offset_(offset)
        {
        }

        // Get the shape of the MemoryView
        std::vector<size_t> shape() const
        {
            return shape_;
        }

        // Get the strides of the MemoryView
        std::vector<size_t> strides() const
        {
            return strides_;
        }

        // Slice method
        MemoryView<T> slice(const std::vector<Slice>& slices) const
        {
            std::vector<size_t> new_shape;
            std::vector<size_t> new_strides;
            size_t new_offset = offset_;

            size_t ndim = shape_.size();
            size_t slice_ndim = slices.size();

            for (size_t i = 0; i < slice_ndim; ++i)
            {
                if (i >= ndim)
                {
                    throw std::invalid_argument("Number of slices exceeds number of dimensions");
                }
                const Slice& s = slices[i];
                if (s.single_index)
                {
                    // Handle single index
                    int idx = s.index;
                    // Handle negative indices
                    if (idx < 0)
                    {
                        idx += static_cast<int>(shape_[i]);
                    }
                    if (idx < 0 || idx >= static_cast<int>(shape_[i]))
                    {
                        throw std::out_of_range("Slice index out of range");
                    }
                    new_offset += idx * strides_[i];
                }
                else
                {
                    // Handle range slice
                    int start = s.has_start ? s.start : 0;
                    int stop = s.has_stop ? s.stop : static_cast<int>(shape_[i]);
                    int step = s.step;

                    // Handle negative indices
                    if (start < 0) start += static_cast<int>(shape_[i]);
                    if (stop < 0) stop += static_cast<int>(shape_[i]);

                    // Clamp start and stop
                    start = std::max(0, std::min(start, static_cast<int>(shape_[i])));
                    stop = std::max(0, std::min(stop, static_cast<int>(shape_[i])));

                    if (step == 0)
                    {
                        throw std::invalid_argument("Slice step cannot be zero");
                    }

                    // Compute the number of elements in this slice
                    int length = 0;
                    if (step > 0)
                    {
                        if (start < stop)
                        {
                            length = (stop - start + step - 1) / step;
                        }
                    }
                    else
                    {
                        if (start > stop)
                        {
                            length = (start - stop - step - 1) / (-step);
                        }
                    }
                    length = std::max(length, 0);
                    new_shape.push_back(static_cast<size_t>(length));
                    new_strides.push_back(strides_[i] * step);
                    new_offset += start * strides_[i];
                }
            }

            // Handle dimensions not sliced (default slicing)
            for (size_t i = slice_ndim; i < ndim; ++i)
            {
                new_shape.push_back(shape_[i]);
                new_strides.push_back(strides_[i]);
            }

            return MemoryView<T>(data_, new_shape, new_strides, new_offset);
        }

        // Reshape method
        void reshape(const std::vector<size_t>& new_shape)
        {
            size_t new_total = compute_total_size(new_shape);
            size_t current_total = compute_total_size(shape_);
            if (new_total != current_total)
            {
                throw std::invalid_argument("Total size must remain unchanged in reshape");
            }
            shape_ = new_shape;
            compute_strides();
        }

        // Overload the [] operator to slice the first dimension
        MemoryView<T> operator[](size_t index) const
        {
            if (shape_.empty())
            {
                throw std::out_of_range("Cannot index into an empty MemoryView");
            }
            if (index >= shape_[0])
            {
                throw std::out_of_range("Index out of range in first dimension");
            }
            size_t new_offset = offset_ + index * strides_[0];
            std::vector<size_t> new_shape(shape_.begin() + 1, shape_.end());
            std::vector<size_t> new_strides(strides_.begin() + 1, strides_.end());
            return MemoryView<T>(data_, new_shape, new_strides, new_offset);
        }

        // Element access
        T& operator()(const std::vector<size_t>& indices) const
        {
            if (indices.size() != shape_.size())
            {
                throw std::invalid_argument("Number of indices does not match number of dimensions");
            }
            size_t idx = offset_;
            for (size_t i = 0; i < indices.size(); ++i)
            {
                if (indices[i] >= shape_[i])
                {
                    throw std::out_of_range("Index out of range in dimension " + std::to_string(i));
                }
                idx += indices[i] * strides_[i];
            }
            return data_[idx];
        }

        // Begin and end iterators (only valid for contiguous memory)
        T* begin() const
        {
            return data_ + offset_;
        }

        T* end() const
        {
            return data_ + offset_ + compute_total_size(shape_);
        }

        // 1. String Representation Method
        std::string to_string(size_t preview_limit = 10) const
        {
            std::ostringstream oss;
            oss << "MemoryView<" << typeid(T).name() << "> "
                << "shape: (";
            for (size_t i = 0; i < shape_.size(); ++i)
            {
                oss << shape_[i];
                if (i != shape_.size() - 1) oss << ", ";
            }
            oss << "), strides: (";
            for (size_t i = 0; i < strides_.size(); ++i)
            {
                oss << strides_[i];
                if (i != strides_.size() - 1) oss << ", ";
            }
            oss << "), offset: " << offset_ << "\nData preview: [";

            size_t total_elements = compute_total_size(shape_);
            size_t display_count = std::min(preview_limit, total_elements);
            for (size_t i = 0; i < display_count; ++i)
            {
                oss << data_[offset_ + i * strides_[0]];
                if (i != display_count - 1) oss << ", ";
            }
            if (display_count < total_elements)
            {
                oss << ", ...";
            }
            oss << "]";

            return oss.str();
        }

        // 2. Stream Output Operator Overload
        friend std::ostream& operator<<(std::ostream& os, const MemoryView<T>& mv)
        {
            os << mv.to_string();
            return os;
        }

        // 3a. Serialization Method
        // Serializes the MemoryView into a std::vector<char>
        // Format:
        // [Type Name Length][Type Name][Number of Dimensions][Shape][Strides][Offset][Data]
        std::vector<char> serialize() const
        {
            std::vector<char> buffer;

            // Serialize Type Name
            std::string type_name = typeid(T).name();
            size_t type_name_length = type_name.size();
            buffer.insert(buffer.end(),
                          reinterpret_cast<const char*>(&type_name_length),
                          reinterpret_cast<const char*>(&type_name_length) + sizeof(size_t));
            buffer.insert(buffer.end(), type_name.begin(), type_name.end());

            // Serialize Number of Dimensions
            size_t ndim = shape_.size();
            buffer.insert(buffer.end(),
                          reinterpret_cast<const char*>(&ndim),
                          reinterpret_cast<const char*>(&ndim) + sizeof(size_t));

            // Serialize Shape
            buffer.insert(buffer.end(),
                          reinterpret_cast<const char*>(shape_.data()),
                          reinterpret_cast<const char*>(shape_.data()) + sizeof(size_t) * ndim);

            // Serialize Strides
            buffer.insert(buffer.end(),
                          reinterpret_cast<const char*>(strides_.data()),
                          reinterpret_cast<const char*>(strides_.data()) + sizeof(size_t) * strides_.size());

            // Serialize Offset
            buffer.insert(buffer.end(),
                          reinterpret_cast<const char*>(&offset_),
                          reinterpret_cast<const char*>(&offset_) + sizeof(size_t));

            // Serialize Data
            size_t data_size = compute_total_size(shape_) * sizeof(T);
            buffer.insert(buffer.end(),
                          reinterpret_cast<const char*>(data_ + offset_),
                          reinterpret_cast<const char*>(data_ + offset_) + data_size);

            return buffer;
        }

        // 3b. Deserialization Method
        // Reconstructs a MemoryView from a std::vector<char>
        // Assumes that the buffer was created by the serialize() method
        static MemoryView<T> deserialize(const std::vector<char>& buffer)
        {
            size_t pos = 0;

            // Deserialize Type Name
            if (pos + sizeof(size_t) > buffer.size())
            {
                throw std::invalid_argument("Buffer too small to contain type name length");
            }
            size_t type_name_length;
            std::memcpy(&type_name_length, buffer.data() + pos, sizeof(size_t));
            pos += sizeof(size_t);

            if (pos + type_name_length > buffer.size())
            {
                throw std::invalid_argument("Buffer too small to contain type name");
            }
            std::string type_name(buffer.data() + pos, type_name_length);
            pos += type_name_length;

            // Check if the type matches
            if (type_name != typeid(T).name())
            {
                throw std::invalid_argument("Type mismatch during deserialization");
            }

            // Deserialize Number of Dimensions
            if (pos + sizeof(size_t) > buffer.size())
            {
                throw std::invalid_argument("Buffer too small to contain number of dimensions");
            }
            size_t ndim;
            std::memcpy(&ndim, buffer.data() + pos, sizeof(size_t));
            pos += sizeof(size_t);

            // Deserialize Shape
            if (pos + sizeof(size_t) * ndim > buffer.size())
            {
                throw std::invalid_argument("Buffer too small to contain shape");
            }
            std::vector<size_t> shape(ndim);
            std::memcpy(shape.data(), buffer.data() + pos, sizeof(size_t) * ndim);
            pos += sizeof(size_t) * ndim;

            // Deserialize Strides
            if (pos + sizeof(size_t) * ndim > buffer.size())
            {
                throw std::invalid_argument("Buffer too small to contain strides");
            }
            std::vector<size_t> strides(ndim);
            std::memcpy(strides.data(), buffer.data() + pos, sizeof(size_t) * ndim);
            pos += sizeof(size_t) * ndim;

            // Deserialize Offset
            if (pos + sizeof(size_t) > buffer.size())
            {
                throw std::invalid_argument("Buffer too small to contain offset");
            }
            size_t offset;
            std::memcpy(&offset, buffer.data() + pos, sizeof(size_t));
            pos += sizeof(size_t);

            // Deserialize Data
            size_t total_elements = compute_total_size(shape);
            size_t data_size = total_elements * sizeof(T);
            if (pos + data_size > buffer.size())
            {
                throw std::invalid_argument("Buffer too small to contain data");
            }

            // Allocate memory for data
            T* data = new T[total_elements];
            std::memcpy(data, buffer.data() + pos, data_size);
            pos += data_size;

            // Create MemoryView
            return MemoryView<T>(data, shape, strides, 0);
        }

        // 4a. Write Serialized Data to FILE*
        void write_to_file(FILE* file) const
        {
            if (!file)
            {
                throw std::invalid_argument("Invalid FILE* provided for writing");
            }
            std::vector<char> buffer = serialize();
            size_t written = fwrite(buffer.data(), sizeof(char), buffer.size(), file);
            if (written != buffer.size())
            {
                throw std::runtime_error("Failed to write all data to file");
            }
        }

        // 4b. Read Serialized Data from FILE* and Reconstruct MemoryView
        static MemoryView<T> read_from_file(FILE* file)
        {
            if (!file)
            {
                throw std::invalid_argument("Invalid FILE* provided for reading");
            }

            // Determine file size
            if (fseek(file, 0, SEEK_END) != 0)
            {
                throw std::runtime_error("Failed to seek to end of file");
            }
            long file_size = ftell(file);
            if (file_size == -1L)
            {
                throw std::runtime_error("Failed to get file size");
            }
            rewind(file);

            // Read file into buffer
            std::vector<char> buffer(static_cast<size_t>(file_size));
            size_t read = fread(buffer.data(), sizeof(char), buffer.size(), file);
            if (read != buffer.size())
            {
                throw std::runtime_error("Failed to read all data from file");
            }

            // Deserialize
            return deserialize(buffer);
        }

    private:
        T* data_;
        std::vector<size_t> shape_;
        std::vector<size_t> strides_;
        size_t offset_;

        // Compute strides based on the current shape
        void compute_strides()
        {
            strides_.resize(shape_.size());
            if (shape_.empty()) return;
            strides_[shape_.size() - 1] = 1;
            for (int i = static_cast<int>(shape_.size()) - 2; i >= 0; --i)
            {
                strides_[i] = strides_[i + 1] * shape_[i + 1];
            }
        }
    };
}
#endif //CMMCORE_MEMORYVIEW_H
