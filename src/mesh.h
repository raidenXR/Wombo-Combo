#include <iostream>
#include <vector>

// using Vec3 = double[3];
struct Vec3 {double x, y, z;};

template <typename T> struct Slice {
    T* ptr;
    size_t len;
    size_t size;
};

template <typename T> Slice<T> alloc (int n) {
    return Slice<T>{
        .ptr = (T*)malloc(n * sizeof(T)),
        .len = (size_t)n,
        .size = (size_t)(n * sizeof(T)),
    };    
}

struct Mesh
{
    std::vector<Vec3> vertices;
    std::vector<uint> indices;

    Mesh (int vertices_len, int indices_len) {
        vertices.reserve(vertices_len);
        indices.reserve(indices_len);
    };

    ~Mesh () {
        vertices.clear();
        indices.clear();
    }
};

