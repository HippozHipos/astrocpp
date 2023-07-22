#ifndef VECTOR3_H
#define VECTOR3_H

#include <cmath>

class Vector4 {
public:
    float x, y, z, w;

    Vector4()
        : x(0), y(0), z(0), w(0) {}

    Vector4(float x, float y, float z, float w)
        : x(x), y(y), z(z), w(w) {}

    // Addition of two vectors
    Vector4 operator+(const Vector4& other) const {
        return Vector4(x + other.x, y + other.y, z + other.z, w + other.w);
    }

    // Subtraction of two vectors
    Vector4 operator-(const Vector3& other) const {
        return Vector4(x - other.x, y - other.y, z - other.z, w - other.w);
    }

    // Assignment operator
    Vector4& operator=(const Vector4& other) {
        x = other.x;
        y = other.y;
        z = other.z;
        w = other.w;
        return *this;
    }

    // Compound assignment operators (+= and -=)
    Vector4& operator+=(const Vector4& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        w += other.w;
        return *this;
    }

    Vector4& operator-=(const Vector4& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        w -= other.w;
        return *this;
    }

    // Scalar multiplication
    Vector4 operator*(float scalar) const {
        return Vector4(x * scalar, y * scalar, z * scalar, w * scalar);
    }

    // Scalar division
    Vector4 operator/(float scalar) const {
        return Vector3(x / scalar, y / scalar, z / scalar, w / scalar);
    }

    // Dot product of two vectors
    float dot(const Vector4& other) const {
        return x * other.x + y * other.y + z * other.z + w * other.w;
    }

    // Cross product of two vectors
    Vector4 cross(const Vector4& other) const {
        return Vector3(
            y * other.z - z * other.y, 
            z * other.w - w * other.z, 
            w * other.x - x * other.w, 
            x * other.y - y * other.x
        );
    }

    // Magnitude (length) of the vector
    float magnitude() const {
        return std::sqrt(x * x + y * y + z * z + w * w);
    }

    // Normalize the vector (make it unit length) or (unit vector)
    Vector4 normalize() const {
        float mag = magnitude();
        return Vector3(x / mag, y / mag, z / mag, w / mag);
    }
};

using vec4 = Vector4;
using Vec4 = Vector4;

#endif //VECTOR3_H