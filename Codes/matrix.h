#pragma once

#include <vector>
#include <ostream>

template<typename T>
class Matrix{

public:
	Matrix(unsigned int width = 0, unsigned int height = 0)
		: width_(width)
		, height_(height)
		, content_(width* height)
	{}

	T& operator()(unsigned int x, unsigned int y) {
		return content_[x + y*width_];
	}

	const T& operator() (unsigned int x, unsigned int y) const {
		return content_[x + y*width_];
	}

	T& operator()(const std::pair<unsigned int,unsigned int>& pos) {
		return content_[pos.first + pos.second*width_];
	}

	const T& operator() (const std::pair<unsigned int,unsigned int>& pos) const {
		return content_[pos.first + pos.second*width_];
	}

	void operator= (const Matrix<T>& other) {
		height_ = other.height_;
		width_ = other.width_;
		content_ = other.content_;
	}

	unsigned int getHeight() const {
		return height_;
	}

	unsigned int getWidth() const {
		return width_;
	}

protected:
	std::vector<T> content_;
	unsigned int height_, width_;
};

template<typename T>
std::ostream& operator<<(std::ostream& stream, Matrix<T> const& obj)
{
	for(unsigned int i=0; i < obj.getHeight(); ++i) {
		for(unsigned int j=0; j < obj.getWidth(); ++j) {
			stream << obj(j,i) << ", ";
		}
		stream << std::endl;
	}
	return stream;
}