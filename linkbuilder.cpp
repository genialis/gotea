#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<math.h>

#include<map>
#include<vector>
#include<utility>

using std::multimap;
using std::map;
using std::vector;
using std::make_pair;

#include"preprocessor.h"
#include"processor.h"

static multimap<const char*, const char*, ltstr> links;
static map<const char*, ptrdiff_t, ltstr> geneset,termset;
static vector<const char*> sorted_links;

void add_link(const char *gene_name, const char *term_name)
{
	links.insert(make_pair(gene_name,term_name));
	geneset[gene_name]=0;
	termset[term_name]=0;
}

int coagulate_links(char **buf, size_t *buf_len)
{
	int tsize=3*sizeof(int),linksize,j;
	struct int_link_t *linktable;
	float_type *logs;

	for (auto &i : geneset)
	{
		if (tsize%sizeof(size_t)) { tsize+=sizeof(size_t)-tsize%sizeof(size_t); }
		i.second=tsize+str_serialize_padsize(); /*string file offset*/
		tsize+=str_serialize_len(i.first);
	}
	for (auto &i : termset)
	{
		if (tsize%sizeof(size_t)) { tsize+=sizeof(size_t)-tsize%sizeof(size_t); }
		i.second=tsize+str_serialize_padsize(); /*string file offset*/
		tsize+=str_serialize_len(i.first);
	}
	if (tsize%4096) { tsize+=4096-tsize%4096; }

	for_each(geneset.begin(),geneset.end(),[] (decltype(geneset)::value_type &v) {
		sorted_links.push_back(v.first);
	});
	sort(sorted_links.begin(),sorted_links.end(),[] (const char *g1, const char *g2) -> bool {
		return geneset[g1]<geneset[g2];
	});

	linksize=(links.size()+geneset.size())*sizeof(struct int_link_t);
	*buf_len=tsize + linksize + (geneset.size()+1)*sizeof(float_type);
	*buf=(char*)malloc(*buf_len);
	logs=(float_type*)((*buf)+tsize+linksize);
	((int*)(*buf))[0]=tsize;
	((int*)(*buf))[1]=geneset.size();
	((int*)(*buf))[2]=links.size();

	for (auto &i : geneset)
	{
		str_serialize((*buf)+i.second,i.first);
	}
	for (auto &i : termset)
	{
		str_serialize((*buf)+i.second,i.first);
	}

	linktable=(struct int_link_t*)((*buf)+tsize);
	j=0;
	for (auto &i : geneset)
	{
		linktable[j].gene.idx=geneset[i.first];
		linktable[j].term.idx=0;
		j++;
	}
	for (auto &link : sorted_links)
	{
		auto range=links.equal_range(link);
		for (auto i=range.first;i!=range.second;++i)
		{
			linktable[j].gene.idx=geneset[i->first];
			linktable[j].term.idx=termset[i->second];
			j++;
		}
	}

	logs[0]=0.0;
	for (size_t i=1;i<=geneset.size();i++)
	{
		logs[i]=logs[i-1]+log((float_type)i);
	}
	return 1;
}

