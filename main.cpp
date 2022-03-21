#include <cmath>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace Geometry {

class Vector;

class IShape;
class Point;
class Segment;
class Ray;
class Line;
class Circle;
class Polygon;

class Vector {
 public:
  Vector();
  Vector(int64_t coord_x, int64_t coord_y);
  Vector(const Vector& vec);
  Vector& operator=(const Vector& vec);
  friend bool operator==(const Vector& vec1, const Vector& vec2);
  friend bool operator!=(const Vector& vec1, const Vector& vec2);
  friend int64_t operator*(const Vector& vec1, const Vector& vec2);
  friend int64_t operator^(const Vector& vec1, const Vector& vec2);
  Vector& operator-=(const Vector& vec);
  Vector& operator+=(const Vector& vec);
  friend Vector operator-(const Vector& vec1, const Vector& vec2);
  friend Vector operator+(const Vector& vec1, const Vector& vec2);
  Vector& operator*=(int64_t num);
  friend Vector operator*(int64_t num, const Vector& vec);
  friend Vector operator*(const Vector& vec, int64_t num);
  Vector operator-() const;

  [[nodiscard]] int64_t GetX() const;
  [[nodiscard]] int64_t GetY() const;
  void SetX(int64_t val);
  void SetY(int64_t val);
  [[nodiscard]] double Length() const;
  friend double Square(const Vector& vec1, const Vector& vec2);

 private:
  int64_t x_;
  int64_t y_;
};

class IShape {
 public:
  virtual IShape& Move(const Vector& vec) = 0;
  [[nodiscard]] virtual bool ContainsPoint(const Point& point) const = 0;
  [[nodiscard]] virtual bool CrossesSegment(const Segment& segment) const = 0;
  [[nodiscard]] virtual IShape* Clone() const = 0;
  [[nodiscard]] virtual std::string ToString() const = 0;
  virtual ~IShape() = default;
};

class Point : public IShape {
 public:
  Point();
  Point(int64_t coord_x, int64_t coord_y);
  Point(const Point& point);
  explicit Point(const Vector& vec);
  Point& operator=(const Point& point);
  friend bool operator==(const Point& point1, const Point& point2);
  friend Vector operator-(const Point& p1, const Point& p2);
  [[nodiscard]] Vector GetVec() const;
  IShape& Move(const Vector& vec) override;
  [[nodiscard]] bool ContainsPoint(const Point& point) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& segment) const override;
  [[nodiscard]] IShape* Clone() const override;
  [[nodiscard]] std::string ToString() const override;

 private:
  Vector vec_;
};

class Segment : public IShape {
 public:
  Segment();
  Segment(const Point& point1, const Point& point2);
  Segment(const Segment& segment);
  Segment(const Vector& vector1, const Vector& vector2);
  Segment& operator=(const Segment& segment);
  [[nodiscard]] Point GetP1() const;
  [[nodiscard]] Point GetP2() const;
  IShape& Move(const Vector& vec) override;
  [[nodiscard]] bool ContainsPoint(const Point& point) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& segment) const override;
  [[nodiscard]] IShape* Clone() const override;
  [[nodiscard]] std::string ToString() const override;

 private:
  Point p1_;
  Point p2_;
};

class Ray : public IShape {
 public:
  Ray();
  Ray(const Point& point1, const Point& point2);
  Ray(const Point& point, const Vector& vec);
  Ray(const Ray& ray);
  Ray& operator=(const Ray& ray);
  [[nodiscard]] Point GetP() const;
  [[nodiscard]] Vector GetVec() const;
  IShape& Move(const Vector& vec) override;
  [[nodiscard]] bool ContainsPoint(const Point& point) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& segment) const override;
  [[nodiscard]] IShape* Clone() const override;
  [[nodiscard]] std::string ToString() const override;

 private:
  Point p_;
  Vector vec_;
};

class Line : public IShape {
 public:
  Line();
  Line(const Point& point1, const Point& point2);
  Line(const Line& line);
  Line(int64_t a, int64_t b, int64_t c);
  Line& operator=(const Line& line);
  [[nodiscard]] int64_t GetA() const;
  [[nodiscard]] int64_t GetB() const;
  [[nodiscard]] int64_t GetC() const;
  [[nodiscard]] Vector GetVec() const;
  [[nodiscard]] static std::pair<int64_t, int64_t> LineIntersect(
      const Line& line1, const Line& line2);
  [[nodiscard]] std::pair<double, double> FindPoint() const;
  static double LineDistance(const Line& line1, const Line& line2);

  IShape& Move(const Vector& vec) override;
  [[nodiscard]] bool ContainsPoint(const Point& point) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& segment) const override;
  [[nodiscard]] IShape* Clone() const override;
  [[nodiscard]] std::string ToString() const override;

 private:
  int64_t A_;
  int64_t B_;
  int64_t C_;
};

class Polygon : public IShape {
 public:
  Polygon();
  explicit Polygon(std::vector<Point> vec);
  Polygon(const Polygon& polygon);
  Polygon& operator=(const Polygon& polygon);
  [[nodiscard]] std::vector<Point> GetVec() const;

  IShape& Move(const Vector& vec) override;
  [[nodiscard]] bool ContainsPoint(const Point& point) const override;
  [[nodiscard]] bool ContainsPointRandom(const Point& point,
                                         const Vector& vec) const;
  [[nodiscard]] bool CrossesSegment(const Segment& segment) const override;
  [[nodiscard]] IShape* Clone() const override;
  [[nodiscard]] std::string ToString() const override;

 private:
  std::vector<Point> vec_;
};

class Circle : public IShape {
 public:
  Circle();
  Circle(const Point& point, uint64_t radius);
  Circle(const Circle& circle);
  Circle& operator=(const Circle& circle);
  [[nodiscard]] Point GetP() const;
  [[nodiscard]] uint64_t GetR() const;

  IShape& Move(const Vector& vec) override;
  [[nodiscard]] bool ContainsPoint(const Point& point) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& segment) const override;
  [[nodiscard]] IShape* Clone() const override;
  [[nodiscard]] std::string ToString() const override;

 private:
  Point p_;
  uint64_t r_;
};

Vector::Vector() {
  x_ = 0;
  y_ = 0;
}

Vector::Vector(int64_t coord_x, int64_t coord_y) {
  x_ = coord_x;
  y_ = coord_y;
}

Vector::Vector(const Vector& vec) {
  x_ = vec.x_;
  y_ = vec.y_;
}

Vector& Vector::operator=(const Vector& vec) = default;

int64_t operator*(const Vector& vec1, const Vector& vec2) {
  return vec1.x_ * vec2.x_ + vec1.y_ * vec2.y_;
}

int64_t operator^(const Vector& vec1, const Vector& vec2) {
  return vec1.x_ * vec2.y_ - vec1.y_ * vec2.x_;
}

Vector& Vector::operator+=(const Vector& vec) {
  x_ += vec.x_;
  y_ += vec.y_;
  return *this;
}

Vector& Vector::operator-=(const Vector& vec) {
  x_ -= vec.x_;
  y_ -= vec.y_;
  return *this;
}

Vector operator+(const Vector& vec1, const Vector& vec2) {
  Vector result(vec1);
  return result += vec2;
}

Vector operator-(const Vector& vec1, const Vector& vec2) {
  Vector result(vec1);
  return result -= vec2;
}
Vector& Vector::operator*=(int64_t num) {
  x_ *= num;
  y_ *= num;
  return *this;
}

Vector operator*(const Vector& vec, int64_t num) {
  Vector result(vec);
  return result *= num;
}

Vector operator*(int64_t num, const Vector& vec) {
  Vector result(vec);
  return result *= num;
}

Vector Vector::operator-() const {
  Vector result;
  result.x_ = -this->x_;
  result.y_ = -this->y_;
  return result;
}
double Vector::Length() const {
  return sqrt(static_cast<double>(this->x_ * this->x_ + this->y_ * this->y_));
}

double Square(const Vector& vec1, const Vector& vec2) {
  return std::abs(static_cast<double>(vec1 ^ vec2) / 2);
}

int64_t Vector::GetX() const { return x_; }

int64_t Vector::GetY() const { return y_; }

void Vector::SetX(int64_t val) { x_ = val; }
void Vector::SetY(int64_t val) { y_ = val; }
bool operator==(const Vector& vec1, const Vector& vec2) {
  return (vec1.x_ == vec2.x_ && vec1.y_ == vec2.y_);
}
bool operator!=(const Vector& vec1, const Vector& vec2) {
  return !(vec1.x_ == vec2.x_ && vec1.y_ == vec2.y_);
}

////////////////////////////////////////////////////////////////////////

Point::Point() { vec_ = Vector(); }

Point::Point(int64_t coord_x, int64_t coord_y) {
  vec_ = Vector(coord_x, coord_y);
}
Point::Point(const Vector& vec) { vec_ = Vector(vec); }

Point::Point(const Point& point) : vec_(point.vec_) {}

Point& Point::operator=(const Point& point) = default;

Vector Point::GetVec() const { return vec_; }
IShape& Point::Move(const Vector& vec) {
  vec_ += vec;
  return *this;
}
bool Point::ContainsPoint(const Point& point) const { return point == *this; }

bool Point::CrossesSegment(const Segment& segment) const {
  return segment.ContainsPoint(*this);
}

Vector operator-(const Point& p1, const Point& p2) {
  return Vector(p1.GetVec() - p2.GetVec());
}
bool operator==(const Point& point1, const Point& point2) {
  return point1.GetVec() == point2.GetVec();
}
IShape* Point::Clone() const {
  auto* p = new Point(*this);
  return p;
}
std::string Point::ToString() const {
  std::string str = "Point(";
  str += std::to_string(vec_.GetX()) + ", " + std::to_string(vec_.GetY()) + ")";
  return str;
}

////////////////////////////////////////////////////////////////////////////////

Segment::Segment() : p1_(), p2_() {}

Segment::Segment(const Point& point1, const Point& point2)
    : p1_(point1), p2_(point2) {}
Segment::Segment(const Segment& segment) : p1_(segment.p1_), p2_(segment.p2_) {}

Segment::Segment(const Vector& vector1, const Vector& vector2)
    : p1_(vector1), p2_(vector2) {}

Segment& Segment::operator=(const Segment& segment) = default;

Point Segment::GetP1() const { return p1_; }

Point Segment::GetP2() const { return p2_; }

IShape& Segment::Move(const Vector& vec) {
  p1_.Move(vec);
  p2_.Move(vec);
  return *this;
}

bool Segment::ContainsPoint(const Point& point) const {
  Vector vec1 = (p2_ - point);
  Vector vec2 = (p1_ - point);
  int64_t min_x = std::min(p1_.GetVec().GetX(), p2_.GetVec().GetX());
  int64_t max_x = std::max(p1_.GetVec().GetX(), p2_.GetVec().GetX());
  int64_t min_y = std::min(p1_.GetVec().GetY(), p2_.GetVec().GetY());
  int64_t max_y = std::max(p1_.GetVec().GetY(), p2_.GetVec().GetY());
  return (vec1 ^ vec2) == 0 &&
         (min_x <= point.GetVec().GetX() && point.GetVec().GetX() <= max_x) &&
         min_y <= point.GetVec().GetY() && point.GetVec().GetY() <= max_y;
}

bool Segment::CrossesSegment(const Segment& segment) const {
  Vector vec1 = p2_ - p1_;
  Vector vec2 = segment.p2_ - segment.p1_;
  Vector vec3 = segment.p1_ - p1_;
  Vector vec4 = segment.p2_ - p1_;
  Vector vec5 = p2_ - segment.p1_;
  return ((vec1 ^ vec4) * (vec1 ^ vec3) < 0 &&
          ((vec2 ^ -vec3) * (vec2 ^ vec5) < 0)) ||
         ((segment.ContainsPoint(p1_) || segment.ContainsPoint(p2_) ||
           ContainsPoint(segment.p1_) || ContainsPoint(segment.p2_)));
}

IShape* Segment::Clone() const {
  auto p = new Segment(*this);
  return p;
}

std::string Segment::ToString() const {
  std::string str = "Segment(";
  str += p1_.ToString() + ", " + p2_.ToString() + ')';
  return str;
}

////////////////////////////////////////////////////////////////////////////////

Ray::Ray() : p_(), vec_() {}

Ray::Ray(const Point& point1, const Point& point2)
    : p_(point1), vec_(point2.GetVec() - point1.GetVec()) {}

Ray::Ray(const Point& point, const Vector& vec) : p_(point), vec_(vec) {}

Ray::Ray(const Ray& ray) : p_(ray.p_), vec_(ray.vec_) {}

Ray& Ray::operator=(const Ray& ray) = default;

Point Ray::GetP() const { return p_; }

Vector Ray::GetVec() const { return vec_; }
IShape& Ray::Move(const Vector& vec) {
  p_.Move(vec);
  return *this;
}
bool Ray::ContainsPoint(const Point& point) const {
  Vector vec = point - p_;
  return (vec ^ vec_) == 0 && (vec * vec_) >= 0;
}

bool Ray::CrossesSegment(const Segment& segment) const {
  Line line1(
      vec_.GetY(), -vec_.GetX(),
      (vec_.GetX() * p_.GetVec().GetY() - (vec_.GetY() * p_.GetVec().GetX())));
  if (line1.CrossesSegment(segment)) {
    Line line2(segment.GetP1(), segment.GetP2());
    if (line1.GetA() * line2.GetB() == line1.GetB() * line2.GetA() &&
        line1.GetB() * line2.GetC() == line1.GetC() * line2.GetB()) {
      return (ContainsPoint(segment.GetP1()) || ContainsPoint(segment.GetP2()));
    }
    auto point = Geometry::Line::LineIntersect(line1, line2);
    int64_t vec_x = point.first - p_.GetVec().GetX();
    int64_t vec_y = point.second - p_.GetVec().GetY();
    Vector vec(vec_x, vec_y);
    return vec * vec_ >= 0;
  }
  return false;
}

IShape* Ray::Clone() const {
  auto res = new Ray(*this);
  return res;
}

std::string Ray::ToString() const {
  std::string str = "Ray(";
  str += p_.ToString() + ", Vector(" + std::to_string(vec_.GetX()) + ", " +
         std::to_string(vec_.GetY()) + "))";
  return str;
}

////////////////////////////////////////////////////////////////////////////////

Line::Line() : A_(0), B_(0), C_(0) {}
Line::Line(const Point& point1, const Point& point2) {
  int64_t x1 = point1.GetVec().GetX();
  int64_t x2 = point2.GetVec().GetX();
  int64_t y1 = point1.GetVec().GetY();
  int64_t y2 = point2.GetVec().GetY();
  A_ = y2 - y1;
  B_ = x1 - x2;
  C_ = x2 * y1 - x1 * y2;
}
Line::Line(const Line& line) : A_(line.A_), B_(line.B_), C_(line.C_) {}

Line::Line(int64_t a, int64_t b, int64_t c) {
  A_ = a;
  B_ = b;
  C_ = c;
}

Line& Line::operator=(const Line& line) = default;

int64_t Line::GetA() const { return A_; }

int64_t Line::GetB() const { return B_; }

int64_t Line::GetC() const { return C_; }

IShape& Line::Move(const Vector& vec) {
  C_ = C_ - vec.GetX() * A_ - vec.GetY() * B_;
  return *this;
}
bool Line::ContainsPoint(const Point& point) const {
  return A_ * point.GetVec().GetX() + B_ * point.GetVec().GetY() + C_ == 0;
}
bool Line::CrossesSegment(const Segment& segment) const {
  int64_t p1 = A_ * segment.GetP1().GetVec().GetX() +
               B_ * segment.GetP1().GetVec().GetY() + C_;
  int64_t p2 = A_ * segment.GetP2().GetVec().GetX() +
               B_ * segment.GetP2().GetVec().GetY() + C_;
  return p1 * p2 <= 0;
}
IShape* Line::Clone() const {
  auto line = new Line(*this);
  return line;
}
std::string Line::ToString() const {
  std::string str = "Line(";
  str += std::to_string(A_) + ", " + std::to_string(B_) + ", " +
         std::to_string(C_) + ')';
  return str;
}
std::pair<int64_t, int64_t> Line::LineIntersect(const Line& line1,
                                                const Line& line2) {
  int64_t denom = line2.A_ * line1.B_ - line1.A_ * line2.B_;
  if (denom == 0) {
    return {0, 0};
  }
  int64_t x = (line1.C_ * line2.B_ - line1.B_ * line2.C_) / denom;
  int64_t y = (line1.A_ * line2.C_ - line1.C_ * line2.A_) / denom;
  return {x, y};
}
Vector Line::GetVec() const {
  Vector vec(-B_, A_);
  return vec;
}

std::pair<double, double> Line::FindPoint() const {
  if (A_ == 0) {
    double y = static_cast<double>(-C_) / static_cast<double>(B_);
    return {0, y};
  }
  double x = static_cast<double>(-C_) / static_cast<double>(A_);
  return {x, 0};
}
double Line::LineDistance(const Line& line1, const Line& line2) {
  std::pair<double, double> point = line1.FindPoint();
  double dist = std::abs(static_cast<double>(line2.A_) * point.first +
                         static_cast<double>(line2.B_) * point.second +
                         static_cast<double>(line2.C_)) /
                sqrt(static_cast<double>(line2.A_ * line2.A_) +
                     static_cast<double>(line2.B_ * line2.B_));
  return dist;
}

////////////////////////////////////////////////////////////////////////////////

Polygon::Polygon() { vec_.emplace_back(0, 0); }

Polygon::Polygon(std::vector<Point> vec) : vec_(std::move(vec)) {}

Polygon::Polygon(const Polygon& polygon) : vec_(polygon.vec_) {}

Polygon& Polygon::operator=(const Polygon& polygon) = default;

std::vector<Point> Polygon::GetVec() const { return vec_; }

IShape& Polygon::Move(const Vector& vec) {
  for (auto& i : vec_) {
    i.Move(vec);
  }
  return *this;
}
bool Polygon::ContainsPoint(const Point& point) const {
  return (ContainsPointRandom(point, Vector(1, 0)) ||
          ContainsPointRandom(point, Vector(rand() % 25 + 1, rand() % 25 + 1)));
}

bool Polygon::ContainsPointRandom(const Point& point, const Vector& vec) const {
  for (size_t i = 0; i < vec_.size() - 1; ++i) {
    if (Segment(vec_[i], vec_[i + 1]).ContainsPoint(point)) {
      return true;
    }
  }
  if (Segment(vec_[vec_.size() - 1], vec_[0]).ContainsPoint(point)) {
    return true;
  }
  Ray ray(point, vec);
  int counter = 0;
  for (size_t i = 0; i < vec_.size() - 1; ++i) {
    if (ray.CrossesSegment(Segment(vec_[i], vec_[i + 1]))) {
      ++counter;
    }
  }
  if (ray.CrossesSegment(Segment(vec_[vec_.size() - 1], vec_[0]))) {
    ++counter;
  }
  return counter % 2 == 1;
}

bool Polygon::CrossesSegment(const Segment& segment) const {
  for (size_t i = 0; i < vec_.size() - 1; ++i) {
    if (Segment(vec_[i], vec_[i + 1]).CrossesSegment(segment)) {
      return true;
    }
  }
  return Segment(vec_[vec_.size() - 1], vec_[0]).CrossesSegment(segment);
}
IShape* Polygon::Clone() const {
  auto res = new Polygon(*this);
  return res;
}
std::string Polygon::ToString() const {
  std::string str = "Polygon(";
  for (size_t i = 0; i < vec_.size(); ++i) {
    str += vec_[i].ToString();
    if (i < vec_.size() - 1) {
      str += ", ";
    }
  }
  str += ')';
  return str;
}

////////////////////////////////////////////////////////////////////////////////

Circle::Circle() {
  p_ = Point(0, 0);
  r_ = 0;
}
Circle::Circle(const Point& point, uint64_t radius) : p_(point), r_(radius) {}

Circle::Circle(const Circle& circle) : p_(circle.p_), r_(circle.r_) {}

Circle& Circle::operator=(const Circle& circle) = default;

Point Circle::GetP() const { return p_; }

uint64_t Circle::GetR() const { return r_; }

IShape& Circle::Move(const Vector& vec) {
  p_.Move(vec);
  return *this;
}
bool Circle::ContainsPoint(const Point& point) const {
  Vector vec = point - p_;
  int64_t dst = vec * vec;
  return dst <= static_cast<int64_t>(r_ * r_);
}

double DstPointPoint(const Point& point1, const Point& point2) {
  Vector vec = point1 - point2;
  return sqrt(
      static_cast<double>(vec.GetX() * vec.GetX() + vec.GetY() * vec.GetY()));
}

double DstPointSegment(const Point& point, const Segment& segment) {
  Point p1 = segment.GetP1();
  Point p2 = segment.GetP2();
  Vector vec1 = point - p1;
  Vector vec2 = point - p2;
  Vector vec3 = p2 - p1;
  if (vec1 * vec3 < 0 || vec2 * (-vec3) < 0) {
    return std::min(DstPointPoint(p1, point), DstPointPoint(p2, point));
  }
  Line line(segment.GetP1(), segment.GetP2());
  return static_cast<double>(std::abs(line.GetA() * point.GetVec().GetX() +
                                      line.GetB() * point.GetVec().GetY() +
                                      line.GetC())) /
         sqrt(static_cast<double>(line.GetA() * line.GetA() +
                                  line.GetB() * line.GetB()));
}

bool Circle::CrossesSegment(const Segment& segment) const {
  Vector vec1 = segment.GetP1() - p_;
  int64_t dst1 = vec1 * vec1;
  Vector vec2 = segment.GetP2() - p_;
  int64_t dst2 = vec2 * vec2;
  if (dst1 < static_cast<int64_t>(r_ * r_) &&
      dst2 < static_cast<int64_t>(r_ * r_)) {
    return false;
  }
  double dst = DstPointSegment(p_, segment);
  return dst <= static_cast<double>(r_);
}
IShape* Circle::Clone() const {
  auto res = new Circle(*this);
  return res;
}
std::string Circle::ToString() const {
  std::string str = "Circle(";
  str += p_.ToString() + ", " + std::to_string(r_) + ')';
  return str;
}

}  // namespace Geometry

template <class T>
void Delete(T* ptr) {
  delete ptr;
}

void CheckFunctions(const Geometry::IShape* shape_ptr,
                    const Geometry::Point& point_a,
                    const Geometry::Point& point_b) {
  std::cout << "Given shape "
            << (shape_ptr->ContainsPoint(point_a) ? "contains"
                                                  : "does not contain")
            << " point A\n";

  const auto kSegmentAb = Geometry::Segment(point_a, point_b);
  std::cout << "Given shape "
            << (shape_ptr->CrossesSegment(kSegmentAb) ? "crosses"
                                                      : "does not cross")
            << " segment AB\n";

  const auto kVectorAb = point_b - point_a;
  const auto kClonedShapePtr =
      shape_ptr->Clone();  // may return either raw or smart pointer
  std::cout << kClonedShapePtr->Move(kVectorAb).ToString();

  Delete(kClonedShapePtr);  // raw pointer compatibility
}

int main() {
  std::unique_ptr<Geometry::IShape> shape_ptr;

  std::string command;
  std::cin >> command;

  int x = 0;
  int y = 0;
  int a = 0;
  int b = 0;

  if (command == "point") {
    std::cin >> x >> y;
    shape_ptr = std::make_unique<Geometry::Point>(x, y);
  } else if (command == "segment") {
    std::cin >> x >> y >> a >> b;
    shape_ptr = std::make_unique<Geometry::Segment>(Geometry::Point(x, y),
                                                    Geometry::Point(a, b));
  } else if (command == "ray") {
    std::cin >> x >> y >> a >> b;
    shape_ptr = std::make_unique<Geometry::Ray>(Geometry::Point(x, y),
                                                Geometry::Point(a, b));
  } else if (command == "line") {
    std::cin >> x >> y >> a >> b;
    shape_ptr = std::make_unique<Geometry::Line>(Geometry::Point(x, y),
                                                 Geometry::Point(a, b));
  } else if (command == "polygon") {
    size_t n_points = 0;
    std::cin >> n_points;
    std::vector<Geometry::Point> points;
    points.reserve(n_points);
    for (size_t i = 0; i < n_points; ++i) {
      std::cin >> x >> y;
      points.emplace_back(x, y);
    }
    shape_ptr = std::make_unique<Geometry::Polygon>(std::move(points));
  } else if (command == "circle") {
    std::cin >> x >> y;
    const auto kCenter = Geometry::Point(x, y);
    uint64_t radius = 0;
    std::cin >> radius;
    shape_ptr = std::make_unique<Geometry::Circle>(kCenter, radius);
  } else {
    std::cerr << "Undefined command" << std::endl;
    return 1;
  }

  std::cin >> x >> y;
  const auto kPointA = Geometry::Point(x, y);
  std::cin >> x >> y;
  const auto kPointB = Geometry::Point(x, y);

  CheckFunctions(shape_ptr.get(), kPointA, kPointB);
  return 0;
}
