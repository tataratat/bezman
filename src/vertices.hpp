/*
MIT License

Copyright (c) 2022 zwar@ilsb.tuwien.ac.at

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef SRC_DYNAMIC_POINT_HPP
#define SRC_DYNAMIC_POINT_HPP

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace bezman {

/*
 * VertexView class
 *
 * Contains the pointer to the Point of interest.
 */
template <typename BaseType = double>
struct VertexView {
 public:
  /// Provide Type for external use
  using Type_ = BaseType;
  using value_type = Type_;
  using ReturnVectorType_ = std::vector<BaseType>;
  using IndexType_ = std::size_t;
  using SmallIndexType_ = std::uint32_t;

  /// As this is a view, it can only be created specified dim
  /// This prohibits copy assignemnt
   Type_* const begin_;
   const IndexType_ dim_;

  // use this inplace of copy assignment
  template <typename IterableType>
  void Replace(IterableType const& rhs) {
    if constexpr (std::is_arithmetic_v<IterableType>) {
      assert(size() == 1);
      *begin_ = rhs;
    } else {
      assert(size() == rhs.size());
      for (SmallIndexType_ i{}; i < dim_; ++i) {
        begin_[i] = rhs[i];
      }
    }
  }

  Type_* begin() const { return begin_; }
  Type_* end() const { return begin_ + dim_; }
  const IndexType_& size() const { return dim_; }
  Type_& operator[](const IndexType_& i) const { return *(begin_ + i); }
  ReturnVectorType_ as_vector() const {
    ReturnVectorType_ return_vector;
    return_vector.reserve(dim_);
    for (SmallIndexType_ i{}; i < dim_; ++i) {
      return_vector.push_back(*(begin_ + i));
    }
    return return_vector;
  };

  /// Addition using RAII
  template <typename IterableType>
  constexpr ReturnVectorType_ operator+(const IterableType& rhs) const {
    assert(rhs.size() == size());

    ReturnVectorType_ sum;
    sum.reserve(dim_);

    for (SmallIndexType_ i{}; i < dim_; ++i) {
      sum.emplace_back(begin_[i] + rhs[i]);
    }

    return sum;
  }

  /// Inplace addition
  template <typename IterableType>
  constexpr VertexView& operator+=(const IterableType& rhs) {
    if constexpr (std::is_arithmetic_v<IterableType>) {
      assert(size() == 1);
      *begin_ += rhs;
    } else {
      if (rhs.size() != size()) {
        std::cout << "rhs " << rhs.size() << " size " << size() << std::endl;
      }

      assert(rhs.size() == size());

      for (SmallIndexType_ i{}; i < dim_; ++i) {
        begin_[i] += rhs[i];
      }
    }

    return (*this);
  }

  /// Subtraction using RAII
  template <typename IterableType>
  constexpr ReturnVectorType_ operator-(const IterableType& rhs) const {
    assert(rhs.size() == size());

    ReturnVectorType_ difference;
    difference.reserve(dim_);

    for (SmallIndexType_ i{}; i < dim_; ++i) {
      difference.emplace_back(begin_[i] - rhs[i]);
    }

    return difference;
  }

  /// Inplace Subtraction
  template <typename IterableType>
  constexpr VertexView& operator-=(const IterableType& rhs) {
    if constexpr (std::is_arithmetic_v<IterableType>) {
      assert(size() == 1);
      *begin_ -= rhs;
    } else {
      assert(rhs.size() == size());

      for (SmallIndexType_ i{}; i < dim_; i++) {
        begin_[i] -= rhs[i];
      }
    }
    return (*this);
  }

  /// Inversion
  constexpr ReturnVectorType_ operator-() const {
    ReturnVectorType_ inverted;
    inverted.reserve(dim_);

    for (SmallIndexType_ i{}; i < dim_; ++i) {
      inverted.emplace_back(begin_[i] * static_cast<Type_>(-1));
    }

    return inverted;
  }

  /// Multiplication
  constexpr ReturnVectorType_ operator*(const Type_& factor) const {
    ReturnVectorType_ multiplied;
    multiplied.reserve(dim_);

    for (SmallIndexType_ i{}; i < dim_; ++i) {
      multiplied.emplace_back(begin_[i] * factor);
    }

    return multiplied;
  }

  /// Inplace multiplication with a scalar
  constexpr VertexView& operator*=(const Type_& factor) {
    for (SmallIndexType_ i{}; i < dim_; ++i) {
      begin_[i] *= factor;
    }
    return (*this);
  }

  /// Friend injection for reversed arguments (facilitates readability)
  friend constexpr ReturnVectorType_ operator*(const Type_& factor,
                                               const VertexView& point) {
    return point * factor;
  }

  /// Division with a scalar
  constexpr VertexView& operator/=(const Type_& factor) {
    const Type_ inv_factor{static_cast<Type_>(1.) / factor};
    (*this) *= inv_factor;
    return (*this);
  }

  /// Division
  constexpr ReturnVectorType_ operator/(const Type_& factor) const {
    ReturnVectorType_ divided;
    divided.reserve(dim_);

    const Type_ inv_factor{static_cast<Type_>(1.) / factor};
    for (SmallIndexType_ i{}; i < dim_; ++i) {
      divided.emplace_back(begin_[i] * inv_factor);
    }

    return divided;
  }

  /// Friend injection for reversed arguments (facilitates readability)
  friend constexpr ReturnVectorType_ operator/(const Type_& factor,
                                               const VertexView& point) {
    ReturnVectorType_ divided;
    divided.reserve(point.dim_);

    for (SmallIndexType_ i{}; i < point.dim_; ++i)
      divided[i] = factor / point[i];
    return divided;
  }

  /// Scalar Product with another Point of same size
  template <typename IterableType>
  auto operator*(const IterableType& point) const {
    assert(point.size() == size());

    using ScalarReturn =
        decltype(Type_{} * typename IterableType::value_type{});

    ScalarReturn result{};
    for (SmallIndexType_ i{}; i < dim_; ++i) {
      result += begin_[i] * point[i];
    }

    return result;
  }

  /// Calculate norm
  Type_ SquaredEuclidianNorm() const {
    Type_ norm{};
    for (SmallIndexType_ i{}; i < dim_; ++i) {
      norm += begin_[i] * begin_[i];
    }
    return norm;
  }

  /// Calculate norm
  Type_ EuclidianNorm() const {
    return std::sqrt((*this).SquaredEuclidianNorm());
  }

  /// Facilitate User Output
  friend std::ostream& operator<<(std::ostream& os, const VertexView& p) {
    os << p.toString();
    return os;
  }

  /// Converts Point into a string object
  std::string toString(const int& output_precision = 6) const {
    std::ostringstream out{};
    out << "[" << std::setw(output_precision + 2)
        << std::setprecision(output_precision) << (*this)[0];
    for (SmallIndexType_ i{1}; i < dim_; ++i) {
      out << ", " << std::setw(output_precision + 2)
          << std::setprecision(output_precision) << (*this)[i];
    }
    out << "]";
    return out.str();
  }
};

template <typename BaseType = double>
struct ConstVertexView {
 public:
  /// Provide Type for external use
  using Type_ = BaseType;
  using value_type = Type_;
  using ReturnVectorType_ = std::vector<BaseType>;
  using IndexType_ = std::size_t;
  using SmallIndexType_ = std::uint32_t;

  /// As this is a view, it can only be created specified dim
  const Type_* begin_;
  const IndexType_ dim_;

  const Type_* begin() const { return begin_; }
  const Type_* end() const { return begin_ + dim_; }
  const IndexType_& size() const { return dim_; }
  const Type_& operator[](const IndexType_& i) const { return *(begin_ + i); }
  ReturnVectorType_ as_vector() const {
    ReturnVectorType_ return_vector;
    return_vector.reserve(dim_);
    for (SmallIndexType_ i{}; i < dim_; ++i) {
      return_vector.push_back(*(begin_ + i));
    }
    return return_vector;
  };

  /// Addition using RAII
  template <typename IterableType>
  constexpr ReturnVectorType_ operator+(const IterableType& rhs) const {
    assert(rhs.size() == size());

    ReturnVectorType_ sum;
    sum.reserve(dim_);

    for (SmallIndexType_ i{}; i < dim_; ++i) {
      sum.emplace_back(begin_[i] + rhs[i]);
    }

    return sum;
  }

  /// Subtraction using RAII
  template <typename IterableType>
  constexpr ReturnVectorType_ operator-(const IterableType& rhs) const {
    assert(rhs.size() == size());

    ReturnVectorType_ difference;
    difference.reserve(dim_);

    for (SmallIndexType_ i{}; i < dim_; ++i) {
      difference.emplace_back(begin_[i] - rhs[i]);
    }

    return difference;
  }

  /// Inversion
  constexpr ReturnVectorType_ operator-() const {
    ReturnVectorType_ inverted;
    inverted.reserve(dim_);

    for (SmallIndexType_ i{}; i < dim_; ++i) {
      inverted.emplace_back(begin_[i] * static_cast<Type_>(-1));
    }

    return inverted;
  }

  /// Multiplication
  constexpr ReturnVectorType_ operator*(const Type_& factor) const {
    ReturnVectorType_ multiplied;
    multiplied.reserve(dim_);

    for (SmallIndexType_ i{}; i < dim_; ++i) {
      multiplied.emplace_back(begin_[i] * factor);
    }

    return multiplied;
  }

  /// Friend injection for reversed arguments (facilitates readability)
  friend constexpr ReturnVectorType_ operator*(const Type_& factor,
                                               const ConstVertexView& point) {
    return point * factor;
  }

  /// Division
  constexpr ReturnVectorType_ operator/(const Type_& factor) const {
    ReturnVectorType_ divided;
    divided.reserve(dim_);

    const Type_ inv_factor{static_cast<Type_>(1.) / factor};
    for (SmallIndexType_ i{}; i < dim_; ++i) {
      divided.emplace_back(begin_[i] * inv_factor);
    }

    return divided;
  }

  /// Friend injection for reversed arguments (facilitates readability)
  friend constexpr ReturnVectorType_ operator/(const Type_& factor,
                                               const ConstVertexView& point) {
    ReturnVectorType_ divided;
    divided.reserve(point.dim_);

    for (SmallIndexType_ i{}; i < point.dim_; ++i)
      divided[i] = factor / point[i];
    return divided;
  }

  /// Scalar Product with another Point of same size
  template <typename IterableType>
  auto operator*(const IterableType& point) const {
    assert(point.size() == size());

    using ScalarReturn =
        decltype(Type_{} * typename IterableType::value_type{});

    ScalarReturn result{};
    for (SmallIndexType_ i{}; i < dim_; ++i) {
      result += begin_[i] * point[i];
    }

    return result;
  }

  /// Calculate norm
  Type_ SquaredEuclidianNorm() const {
    Type_ norm{};
    for (SmallIndexType_ i{}; i < dim_; ++i) {
      norm += begin_[i] * begin_[i];
    }
    return norm;
  }

  /// Calculate norm
  Type_ EuclidianNorm() const {
    return std::sqrt((*this).SquaredEuclidianNorm());
  }

  /// Facilitate User Output
  friend std::ostream& operator<<(std::ostream& os, const ConstVertexView& p) {
    os << p.toString();
    return os;
  }

  /// Converts Point into a string object
  std::string toString(const int& output_precision = 6) const {
    std::ostringstream out{};
    out << "[" << std::setw(output_precision + 2)
        << std::setprecision(output_precision) << (*this)[0];
    for (SmallIndexType_ i{1}; i < dim_; ++i) {
      out << ", " << std::setw(output_precision + 2)
          << std::setprecision(output_precision) << (*this)[i];
    }
    out << "]";
    return out.str();
  }
};

template <typename Type>
class Vertices {
 public:
  using Type_ = Type;
  using Scalar = Type_;  // for compatibility
  using IndexType_ = std::size_t;
  using VertexView_ = VertexView<Type_>;
  using ConstVertexView_ = ConstVertexView<Type_>;
  using ContainerType_ = std::vector<Type_>;
  using ReturnVectorType_ = typename VertexView_::ReturnVectorType_;

 protected:
  ContainerType_ vertices_;
  IndexType_ dimension_;

 public:
  Vertices() = default;
  explicit Vertices(ContainerType_ vertices, const IndexType_ dimension)
      : vertices_(std::move(vertices)) {
    SetDimension(dimension);
  }
  Vertices(Vertices const& other) = default;
  Vertices(Vertices&& other) noexcept = default;
  Vertices& operator=(Vertices const& rhs) = default;
  Vertices& operator=(Vertices&& rhs) noexcept = default;
  ~Vertices() = default;

  bool operator==(Vertices const& rhs) const {
    return rhs.GetVertices() == GetVertices();
  }

  void SetDimension(const IndexType_ dimension) {
    if ((vertices_.size() % dimension) != 0) {
      // although, this allows 0 size vertices for init
      throw std::runtime_error(
          "bezman::Vertices - Invalid dimension. Size of the vertices should "
          "be a multiple of the dimension.");
    }
    dimension_ = dimension;
  }

  IndexType_ GetDimension() const {
    if (dimension_ != 0) {
      return dimension_;
    } else {
      throw std::runtime_error("bezman::Vertices - Dimension not set.");
    }
  }

  ContainerType_& GetVertices() { return vertices_; }
  const ContainerType_& GetVertices() const { return vertices_; }

  VertexView_ operator[](const IndexType_ id) {
    const auto dim = GetDimension();

    return {GetVertices().data() + (id * dim), dim};
  };

  ConstVertexView_ operator[](const IndexType_ id) const {
    const auto dim = GetDimension();
    return {GetVertices().data() + (id * dim), dim};
  };

  IndexType_ GetNumberOfVertices() const {
    return vertices_.size() / GetDimension();
  }

  IndexType_ size() const { return GetNumberOfVertices(); }

  void resize(const IndexType_ new_size) {
    vertices_.resize(new_size * GetDimension());
  };
};

namespace vector_plus_vector {

template <typename T>
std::vector<T> operator+(std::vector<T> const& a, std::vector<T> const& b) {
  const auto a_size = a.size();
  assert(a.size() == b.size());

  std::vector<T> c;
  c.reserve(a_size);

  for (decltype(a.size()) i{}; i < a_size; ++i) {
    c.emplace_back(a[i] + b[i]);
  }

  return c;
}

}  // namespace vector_plus_vector

namespace scalar_times_vector {

template <typename T>
std::vector<T> operator*(T const& a, std::vector<T> const& b) {
  const auto b_size = b.size();
  std::vector<T> c;
  c.reserve(b_size);

  for (decltype(b.size()) i{}; i < b_size; ++i) {
    c.emplace_back(a * b[i]);
  }

  return c;
}

}  // namespace scalar_times_vector

namespace vector_times_scalar {

template <typename T>
std::vector<T> operator*(std::vector<T> const& a, T const& b) {
  const auto a_size = a.size();
  std::vector<T> c;
  c.reserve(a_size);

  for (decltype(a.size()) i{}; i < a_size; ++i) {
    c.emplace_back(a[i] * b);
  }

  return c;
}

}  // namespace vector_times_scalar

namespace vector_times_vector {

template <typename T>
std::vector<T> operator*(std::vector<T> const& a, std::vector<T> const& b) {
  const auto a_size = a.size();
  assert(a.size() == b.size());

  std::vector<T> c;
  c.reserve(a_size);

  for (decltype(a.size()) i{}; i < a_size; ++i) {
    c.emplace_back(a[i] * b[i]);
  }

  return c;
}

}  // namespace vector_times_vector

namespace vector_minus_vector_inplace {

template <typename T>
void operator-=(std::vector<T>& a, std::vector<T> const& b) {
  const auto a_size = a.size();
  assert(a.size() == b.size());

  for (decltype(a.size()) i{}; i < a_size; ++i) {
    a[i] -= b[i];
  }
}

}  // namespace vector_minus_vector_inplace

namespace vector_times_scalar_inplace {

template <typename T>
void operator*=(std::vector<T>& a, T const& b) {
  const auto a_size = a.size();

  for (decltype(a.size()) i{}; i < a_size; ++i) {
    a[i] *= b;
  }
}

}  // namespace vector_times_scalar_inplace

}  // namespace bezman

#endif  // SRC_DYNAMIC_POINT_HPP
