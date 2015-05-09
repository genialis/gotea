#ifndef __PROCESSOR_H__
#define __PROCESSOR_H__

#include<stddef.h>
#include<unordered_set>
#include<vector>
#include<algorithm>

/* string library; hashing, serialization */

struct StrComp
{
	bool operator()(const char *str1, const char *str2) const
	{
		return strcmp(str1,str2)<0;
	}
};

struct StrHash
{
	size_t operator()(const char *s) const
	{
		return StrHash::fakehash(s);
	}
	static size_t fakehash(const char *s)
	{
		return ((size_t*)s)[-1];
	}

	static size_t realhash(const char *s)
	{
		/*FNV-1a 32- or 64bit*/
		static const size_t prime=(sizeof(size_t)==8)?1099511628211LU:16777619LU;
		static const size_t offset=(sizeof(size_t)==8)?14695981039346656037LU:2166136261LU;
		size_t accum=offset;
		for (int i=0;s[i]!=0;i++)
		{
			accum^=s[i];
			accum*=prime;
		}
		return accum;
	}
};
struct PtrHash
{
	size_t operator()(const char *s) const { return (size_t)s; }
};

struct StrEq
{
	bool operator()(const char *str1, const char *str2) const
	{
		return strcmp(str1,str2)==0;
	}
};
struct PtrEq
{
	bool operator()(const char *str1, const char *str2) const { return str1==str2; }
};

static inline constexpr int str_serialize_padsize() { return 2*sizeof(size_t); }
static inline int str_serialize_len(const char *s)
{
	return strlen(s)+str_serialize_padsize();
}
static inline void str_serialize(char *buf, const char *s)
{
	size_t *nums=(size_t*)buf;
	nums[-1]=StrHash::realhash(s);
	nums[-2]=strlen(s);
	strcpy(buf,s);
}

static inline size_t strlen_bin(const char *s)
{
	return ((size_t*)s)[-2];
}

/*elements must always be pushed in ascending order*/
struct quick_set_t : public std::vector<const char*>
{
	quick_set_t(size_t res)
	{
		reserve(res);
	}
	bool push_gene(const char *gene)
	{
		if (size()==0 || back()!=gene)
		{
			push_back(gene);
			return true;
		}
		return false;
	}
	bool has_gene(const char *gene) const
	{
		return std::binary_search(begin(),end(),gene,[] (const char *g1, const char *g2) -> bool { return g1<g2; });
	}
};

/*struct quick_set_t : public std::unordered_set<const char*, PtrHash, PtrEq>
{
	quick_set_t(size_t) {}
	bool push_gene(const char *gene)
	{
		return insert(gene).second;
	}
	bool has_gene(const char *gene) const
	{
		return count(gene)!=0;
	}
};*/

typedef double float_type;

struct int_term_t
{
	union
	{
		struct int_term_t **ptr;
		ptrdiff_t idx;
	} parents,children;
	int nparents,nchildren;
	union
	{
		const char *str;
		ptrdiff_t idx;
	} id,name,prerender;
	struct quick_set_t genes,intersect;
	float_type pval,score;
	int prerender_len,genes_len,intersect_len,dumped;
	long color;
	char *genebuf;
};

struct int_link_t
{
	union
	{
		const char *str;
		ptrdiff_t idx;
	} gene,term;
};

#endif /* __PROCESSOR_H__ */
