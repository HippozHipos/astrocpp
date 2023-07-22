#ifndef VECTOR2_H
#define VECTOR2_H

#include <cmath>

class Vector2 {
public:
    float x, y;

    Vector2()
        : x(0), y(0) {}

    Vector2(float x, float y)
        : x(x), y(y) {}

    // Addition of two vectors
    Vector2 operator+(const Vector2& other) const { return Vector2(x + other.x, y + other.y);   }

    // Subtraction of two vectors
    Vector2 operator-(const Vector2& other) const { return Vector2(x - other.x, y - other.y);   }

    // Assignment operator
    Vector2& operator=(const Vector2& other) {
        x = other.x;
        y = other.y;
        return *this;
    }

    // Compound assignment operators (+= and -=)
    Vector2& operator+=(const Vector2& other) {
        x += other.x;
        y += other.y;
        return *this;
    }

    Vector2& operator-=(const Vector2& other) {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    // Scalar multiplication
    Vector2 operator*(float scalar) const {
        return Vector3(x * scalar, y * scalar, z * scalar);
    }

    // Scalar division
    Vector2 operator/(float scalar) const {
        return Vector3(x / scalar, y / scalar, z / scalar);
    }

    // Dot product of two vectors
    float dot(const Vector2& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // Cross product of two vectors
    Vector2 cross(const Vector2& other) const {
        return Vector2(
            y * other.x - x * other.y, 
            x * other.y - y * other.x
        );
    }

    // Magnitude (length) of the vector (sqrt(x^2 + y^2))
    float magnitude() const {
        return std::sqrt(x * x + y * y);
    }

    // Normalize the vector (make it unit length) or (unit vector)
    Vector3 normalize() const {
        float mag = magnitude();
        return Vector3(x / mag, y / mag);
    }
};

using vec2 = Vector2;
using Vec2 = Vector2;

#endif