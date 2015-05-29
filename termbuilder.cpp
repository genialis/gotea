#include<string.h>
#include<stdlib.h>
#include<stddef.h>
#include<syslog.h>

#include<stdio.h>
#include<map>
using std::map;

#include"preprocessor.h"
#include"processor.h"

/*local allocator*/

static struct term_t *pool;
static int pool_size;

struct term_t *new_term()
{
	if (pool_size==0)
	{
		pool_size=100000;
		pool=(struct term_t*)malloc(sizeof(*pool)*pool_size);
	}
	pool_size--;
	new (pool) term_t();
	return pool++;
}

/*tree management*/

static map<const char*, struct term_t*, ltstr> tree;
static map<const char*, ptrdiff_t, ltstr> stringset;

void add_term(struct term_t *term)
{
	if (term==NULL) { return; }
	if (term==NULL || term->obsolete==true || tree.count(term->name)>0) { return; }
	tree[term->id]=term;
}

int link_terms()
{
	/*preprocess isa_for cache first; the inverses of these _should_ be subsets of the respective isa's*/
	for (auto &i : tree)
	{
		struct term_t *t=i.second;
		for (auto &j : t->isa_for)
		{
			if (tree.count(j)==0) { syslog(LOG_ERR,"term %s has nonexistent child '%s'",t->id,j); return 0; }
			tree[j]->isa.insert(i.first);
		}
	}

	/*now link term_t* parents and children by names*/
	for (auto &i : tree)
	{
		struct term_t *t=i.second,*p;
		for (auto &j : t->isa)
		{
			if (tree.count(j)==0) { syslog(LOG_ERR,"term %s has nonexistent parent '%s'",t->id,j); return 0; }
			p=tree[j];
			t->parents.insert(p);
			p->children.insert(t);
		}
		stringset[t->id]=0;
		stringset[t->name]=0;
	}
	return 1;
}

int coagulate_terms(char **buf, size_t *buf_len)
{
	int strlens=2*sizeof(int),linklens,termlens,renderlens;
	char *renderbuf;

	/*string index*/
	for (auto &i : stringset)
	{
		if (strlens%sizeof(size_t)) { strlens+=sizeof(size_t)-strlens%sizeof(size_t); }
		i.second=strlens+str_serialize_padsize();
		strlens+=str_serialize_len(i.first);
	}
	strlens=align(strlens,4096);

	/*term linkage arrays*/
	linklens=strlens;
	for (auto &i : tree)
	{
		struct term_t *t=i.second;
		t->i_parents=linklens;
		linklens+=t->parents.size()*sizeof(ptrdiff_t);
		t->i_children=linklens;
		linklens+=t->children.size()*sizeof(ptrdiff_t);
	}
	linklens=align(linklens,4096);

	/*term indexing*/
	termlens=linklens;
	int arr_index=0;
	for (auto &i : tree)
	{
		struct term_t *t=i.second;
		t->i_this=termlens;
		termlens+=sizeof(struct int_term_t);
		t->arr_index=arr_index++;
	}
	termlens=align(termlens,4096);

	/*get buffer size for the json pre-render*/
	renderlens=termlens;
	for (auto &i : tree)
	{
		struct term_t *t=i.second;
		renderlens += 200 + strlen(t->id) + strlen(t->name); //200 should encompass all the below rendered text
	}

	/*build the actual buffer*/

	*buf_len=renderlens;
	*buf=(char*)malloc(*buf_len);
	memset(*buf,0,*buf_len);

	((int*)(*buf))[0]=linklens;
	((int*)(*buf))[1]=tree.size();

	/*build the string section*/
	for (auto &i : stringset)
	{
		str_serialize((*buf)+i.second,i.first);
	}

	/*build the term linkage tables*/
	for (auto &i : tree)
	{
		ptrdiff_t *arr;
		struct term_t *t=i.second;

		arr=(ptrdiff_t*)((*buf)+t->i_parents);
		for (auto &j : t->parents)
		{
			*(arr++)=j->i_this;
		}

		arr=(ptrdiff_t*)((*buf)+t->i_children);
		for (auto &j : t->children)
		{
			*(arr++)=j->i_this;
		}
	}

	/*build the actual term list*/
	struct int_term_t *it=(struct int_term_t*)((*buf)+linklens);
	renderbuf=*buf+termlens;
	static const char *intplaceholder="        "; // 8 chars
	static const char *floatplaceholder="                        "; // 24 chars
	for (auto &i : tree)
	{
		struct term_t *t=i.second;
		it->parents.idx=t->i_parents;
		it->nparents=t->parents.size();
		it->children.idx=t->i_children;
		it->nchildren=t->children.size();
		it->id.idx=stringset[t->id];
		it->name.idx=stringset[t->name];

		int renderlen;
		it->prerender.idx=renderbuf-*buf;
		renderlen=sprintf(renderbuf,
			"{\"matched\":     %s" // +16
			",\"pval\":%s"         // +32
			",\"score\":       %s" // +72
			",\"total\":       %s" // +112
			",\"term_id\":\"%s\",\"term_name\":\"%s\",\"gene_ids\":[",
			intplaceholder,
			floatplaceholder,
			floatplaceholder,
			intplaceholder,
			t->id,t->name
		);
		it->prerender_len=renderlen;
		renderbuf+=align(renderlen,64);

		it++;
	}

	return 1;
}

