#ifndef VECTOR3_H
#define VECTOR3_H

#include <cmath>

class Vector3 {
public:
    float x, y, z;

    Vector3()
        : x(0), y(0), z(0) {}

    Vector3(float x, float y, float z)
        : x(x), y(y), z(z) {}

    // Addition of two vectors
    Vector3 operator+(const Vector3& other) const {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    // Subtraction of two vectors
    Vector3 operator-(const Vector3& other) const {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    // Assignment operator
    Vector3& operator=(const Vector3& other) {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }

    // Compound assignment operators (+= and -=)
    Vector3& operator+=(const Vector3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    Vector3& operator-=(const Vector3& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    // Scalar multiplication
    Vector3 operator*(float scalar) const {
        return Vector3(x * scalar, y * scalar, z * scalar);
    }

    // Scalar division
    Vector3 operator/(float scalar) const {
        return Vector3(x / scalar, y / scalar, z / scalar);
    }

    // Dot product of two vectors
    float dot(const Vector3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // Cross product of two vectors
    Vector3 cross(const Vector3& other) const {
        return Vector3(
            y * other.z - z * other.y, 
            z * other.x - x * other.z, 
            x * other.y - y * other.x
        );
    }

    // Magnitude (length) of the vector
    float magnitude() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    // Normalize the vector (make it unit length) or (unit vector)
    Vector3 normalize() const {
        float mag = magnitude();
        return Vector3(x / mag, y / mag, z / mag);
    }
};

using vec3 = Vector3;
using Vec3 = Vector3;

#endif //VECTOR3_H