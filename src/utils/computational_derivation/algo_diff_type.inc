#ifndef UTILS_COMPUTATIONAL_DERIVATION_ALGO_DIFF_TYPE_INC
#define UTILS_COMPUTATIONAL_DERIVATION_ALGO_DIFF_TYPE_INC

// Addition
template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator+(
    const AlgoDiffType<Scalar>& b) const {
  ADT_ res{(*this)};
  res += b;
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator+=(
    const AlgoDiffType<Scalar>& b) {
  assert(b.GetNumberOfDerivatives() == GetNumberOfDerivatives());
  v_ += b.v_;
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    d_[i] += b.d_[i];
  }
  return (*this);
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator+(
    const Scalar& b) const {
  ADT_ res{(*this)};
  res += b;
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator+=(
    const Scalar& b) {
  v_ += b;
  return (*this);
}

// Substraction and negation
template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator-() const {
  ADT_ res{(*this)};
  res.v_ = -res.v_;
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    res.d_[i] = -res.d_[i];
  }
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator-(
    const AlgoDiffType<Scalar>& b) const {
  assert(b.GetNumberOfDerivatives() == GetNumberOfDerivatives());
  ADT_ res{(*this)};
  res -= b;
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator-=(
    const AlgoDiffType<Scalar>& b) {
  v_ -= b.v_;
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    d_[i] -= b.d_[i];
  }
  return (*this);
}
template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator-(
    const Scalar& b) const {
  ADT_ res{(*this)};
  res -= b;
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator-=(
    const Scalar& b) {
  v_ -= b;
  return (*this);
}

// Multiplication
template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator*(
    const AlgoDiffType<Scalar>& b) const {
  ADT_ res{(*this)};
  res *= b;
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator*=(
    const AlgoDiffType<Scalar>& b) {
  assert(b.GetNumberOfDerivatives() == GetNumberOfDerivatives());
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    d_[i] = b.d_[i] * v_ + b.v_ * d_[i];  // f'(x) = u'*v +  u*v'
  }
  // Must be at the end otherwise overwrite value
  v_ *= b.v_;
  return (*this);
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator*(
    const Scalar& b) const {
  ADT_ res{(*this)};
  res *= b;
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator*=(
    const Scalar& b) {
  v_ *= b;
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    d_[i] *= b;
  }
  return (*this);
}

// Division
template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator/(
    const AlgoDiffType<Scalar>& b) const {
  ADT_ res{(*this)};
  res /= b;
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator/=(
    const AlgoDiffType<Scalar>& b) {
  assert(b.GetNumberOfDerivatives() == GetNumberOfDerivatives());
  const Scalar_ inverse_b_squared = 1.0 / (b.v_ * b.v_);
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    // f'(x) = (u'*v -  u*v') / v^2
    d_[i] = (b.v_ * d_[i] - b.d_[i] * v_) * inverse_b_squared;
  }
  v_ /= b.v_;
  return (*this);
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator/(
    const Scalar& b) const {
  ADT_ res{(*this)};
  res /= b;
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator/=(
    const Scalar& b) {
  const Scalar inverse_b{1 / b};
  v_ *= inverse_b;
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    d_[i] *= inverse_b;
  }
  return (*this);
}

//////////////////////////////
// Friend Functions
//////////////////////////////

template <typename Scalar>
std::ostream& operator<<(std::ostream& os, const AlgoDiffType<Scalar>& a) {
  os << a.v_;
  return os;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> operator+(const Scalar& a,
                                         const AlgoDiffType<Scalar>& b) {
  return b + a;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> operator-(const Scalar& a,
                                         const AlgoDiffType<Scalar>& b) {
  return (-b) + a;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> operator*(const Scalar& a,
                                         const AlgoDiffType<Scalar>& b) {
  return b * a;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> operator/(const Scalar& a,
                                         const AlgoDiffType<Scalar>& b) {
  constexpr AlgoDiffType<Scalar> one_type{a, b.GetNumberOfDerivatives()};
  // A few too many operations (multiplication with 0), but typesafe
  return one_type / b;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> exp(const AlgoDiffType<Scalar>& exponent) {
  Scalar tmp;
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    tmp = std::exp(exponent.v_);
  } else {
    tmp = exp(exponent.v_);
  }
  AlgoDiffType<Scalar> res{exponent};
  res.v_ = Scalar{1.};
  res *= tmp;
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> abs(const AlgoDiffType<Scalar>& base) {
  return base.v_ > 0.0 ? base : (-base);
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> pow(const AlgoDiffType<Scalar>& base,
                                   const Scalar& power) {
  Scalar tmp;
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    tmp = std::pow(base.v_, power - 1.0);
  } else {
    tmp = pow(base.v_, power - 1.0);
  }
  AlgoDiffType<Scalar> res{base};
  res.v_ *= tmp;
  for (std::size_t i{}; i < base.GetNumberOfDerivatives(); i++) {
    res.d_[i] *= power * tmp;
  }
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> pow(const AlgoDiffType<Scalar>& base,
                                   const AlgoDiffType<Scalar>& power) {
  return exp(log(base) * power);
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> sqrt(const AlgoDiffType<Scalar>& radicand) {
  Scalar root;
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    root = std::sqrt(radicand.v_);
  } else {
    root = sqrt(radicand.v_);
  }
  const Scalar& half_inverse_root = 0.5 / root;
  AlgoDiffType<Scalar> res{radicand};
  res.v_ = root;
  for (std::size_t i{}; i < radicand.GetNumberOfDerivatives(); i++) {
    res.d_[i] *= half_inverse_root;
  }
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> log(const AlgoDiffType<Scalar>& xi) {
  using namespace std;               // required for intrinsic types
  const Scalar& temp = 1.0 / xi.v_;  // Multiplication is cheaper than division
  AlgoDiffType<Scalar> res{xi};
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    res.v_ = std::log(xi.v_);
  } else {
    res.v_ = log(xi.v_);
  }
  for (std::size_t i{}; i < xi.GetNumberOfDerivatives(); i++) {
    res.d_[i] *= temp;
  }
  return res;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> log10(const AlgoDiffType<Scalar>& a) {
  using namespace std;  // requ ired for intrinsic types
  const Scalar& temp = 1.0 / (a.v_ * std::log(10.));
  AlgoDiffType<Scalar> res{a};
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    res.v_ = std::log10(a.v_);
  } else {
    res.v_ = log10(a.v_);
  }
  for (std::size_t i{}; i < a.GetNumberOfDerivatives(); i++) {
    res.d_[i] *= temp;
  }
  return res;
};

template <typename Scalar>
constexpr bool operator>(const Scalar& scalar,
                         const AlgoDiffType<Scalar>& adt) {
  return adt < scalar;
};

template <typename Scalar>
constexpr bool operator>=(const Scalar& scalar,
                          const AlgoDiffType<Scalar>& adt) {
  return adt <= scalar;
};

template <typename Scalar>
constexpr bool operator<(const Scalar& scalar,
                         const AlgoDiffType<Scalar>& adt) {
  return adt > scalar;
};

template <typename Scalar>
constexpr bool operator<=(const Scalar& scalar,
                          const AlgoDiffType<Scalar>& adt) {
  return adt >= scalar;
};

template <typename Scalar>
constexpr bool operator==(const Scalar& scalar,
                          const AlgoDiffType<Scalar>& adt) {
  return adt == scalar;
};

template <typename Scalar>
constexpr bool operator!=(const Scalar& scalar,
                          const AlgoDiffType<Scalar>& adt) {
  return !(adt == scalar);
};

#endif  // UTILS_COMPUTATIONAL_DERIVATION_ALGO_DIFF_TYPE_INC
