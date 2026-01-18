/*
Copyright (C) 2020-2022 Michele Martone

This file is part of librsb.

librsb is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

librsb is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with librsb; see the file COPYING.
If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef OPENMP_ALLOCATOR_HPP_INCLUDED
#define OPENMP_ALLOCATOR_HPP_INCLUDED

template <class T> class OpenMP_Allocator
{
	public:
	typedef T* pointer; // since C++17, unnecessary
	typedef const T* const_pointer;
	typedef T& reference;
	typedef const T& const_reference;
	typedef size_t size_type;
	typedef ptrdiff_t difference_type;
	typedef T value_type;
	OpenMP_Allocator(void) {}
	OpenMP_Allocator(const OpenMP_Allocator&) {}
	~OpenMP_Allocator(void) {}

	pointer allocate(size_type numObjects)
	{
		size_type len = numObjects * sizeof(value_type);
		pointer p = static_cast<T*>(std::malloc(len));
#ifndef _OPENMP
#define omp_in_parallel() 0
#endif
		if(!omp_in_parallel())
#pragma omp parallel for schedule(static)
			for(size_type i=0; i<numObjects ; i++)
				p[i]={};
		return p;
#ifndef _OPENMP
#undef omp_in_parallel
#endif
	}

	void deallocate(pointer ptrToMemory, size_type)
	{
		std::free(ptrToMemory);
	}

	void construct(pointer p, const value_type& x)
	{
		new(p) value_type(x);
	}

	void destroy(pointer p)
	{
		p->~value_type();
	}
	private:
	void operator=(const OpenMP_Allocator&) {}
	public:
	bool operator!=(const OpenMP_Allocator&) const { return true; }
};

#endif /* OPENMP_ALLOCATOR_HPP_INCLUDED */
