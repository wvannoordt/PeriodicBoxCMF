#ifndef MD_ARRAY_H
#define MD_ARRAY_H

#include <type_traits>

template <const int... vals> struct Dims
{
	template <const int depth, typename ...T> static constexpr size_t prodr(void) {return 1;}
	template <const int depth, const int T, const int... Tss> static constexpr size_t prodr(void) {return (depth==0)?T:T*prodr<depth-1, Tss...>();}
	template <const int depth, const int... Tss> static constexpr size_t prod(void) {return (depth==0)?1:prodr<depth-1, Tss...>();}
	static constexpr size_t size(void)
	{
		return prod<sizeof...(vals), vals...>();
	}
	template <const int depth, typename index> static inline size_t idxR(index i)
	{
		static_assert(std::is_integral<index>::value, "Integral value for index required.");
		return i*prod<sizeof...(vals)-1, vals...>();
	}
	template <const int depth, typename index, typename... indices> static inline size_t idxR(index i, indices... is)
	{
		static_assert(std::is_integral<index>::value, "Integral value for index required.");
		return prod<depth, vals...>()*i + idxR<depth+1>(is...);
	}
	template <typename... indices> static inline size_t offset(indices... is)
	{
		return idxR<0, indices...>(is...);
	}
};
template <class BaseType, class ArrayDims, class BaseContainer=BaseType[ArrayDims::size()]> struct DataView
{
	BaseContainer data;
	DataView(void){}
	template <class rhsType> DataView(const rhsType rhs)
	{
		data = rhs;
	}
	template <typename... indices> inline BaseType& operator () (indices... is)
	{
		return data[ArrayDims::offset(is...)];
	}
};

#endif