#ifndef NUMCPP_H
#define NUMCPP_H


namespace numcpp
{



	template <class T, int X, int Y>
	class Matrix 
	{
	public:
		//Big three
		Matrix() = default;
		Matrix(const Matrix& o);
		Matrix()~;

		//operators
		Matrix& operator= (const Matrix& o);
		T& operator[] (const int i) { return this->type[i]; }

		//Functions
		void fillZeros() {for (int i = 0; i < data.size(); ++i) { data[i] = 0; }

        // Addition of two vectors
        Matrix operator+(const Matrix & other) const {
            return Matrix(x + other.x, y + other.y, z + other.z);
        }

        // Subtraction of two vectors
        Matrix operator-(const Matrix & other) const {
            return Matrix(x - other.x, y - other.y, z - other.z);
        }

        // Assignment operator
        Matrix& operator=(const Matrix& other) {
            for (int i : data) data[i] = other.i;
            return *this;
        }

        // Compound assignment operators (+= and -=)
        Matrix& operator+=(const Matrix& other) {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        Matrix& operator-=(const Matrix& other) {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        // Scalar multiplication
        Matrix operator*(float scalar) const {
            return Matrix(x * scalar, y * scalar, z * scalar);
        }

        // Scalar division
        Matrix operator/(float scalar) const {
            return Matrix(x / scalar, y / scalar, z / scalar);
        }

        // Dot product of two vectors
        float dot(const Vector3 & other) const {
            return x * other.x + y * other.y + z * other.z;
        }

        // Cross product of two vectors
        Matrix cross(const Matrix& other) const {
            return Matrix(
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
        Matrix normalize() const {
            float mag = magnitude();
            return Matrix(x / mag, y / mag, z / mag);
        }
        };


	private:
		const T[X * Y] data;
	};
	
	template <class T>
	using Vector3<T> = Matrix<T, 3, 3>;

};

#endif //NUMCPP_H