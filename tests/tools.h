#ifndef RBDYN_TESTS_TOOLS_H
#define RBDYN_TESTS_TOOLS_H

template <class T>
inline void resetVector(std::vector<T> & v, unsigned size, const T& t)
{
	v.resize(size);
	for (unsigned i=0;i<size;++i)
		v[i] = t;
}

#endif //RBDYN_TESTS_TOOLS_H

