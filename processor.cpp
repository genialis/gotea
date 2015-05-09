#include<stdio.h>
#include<string.h>
#include<sys/mman.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>
#include<fcntl.h>
#include<syslog.h>
#include<errno.h>
#include<math.h>

#include<omp.h>

#include<unordered_map>
#include<vector>
#include<set>
#include<algorithm>
#include<queue>

#include"processor.h"

#include"formatters.cpp"

//#define DEBUG

using std::unordered_set;
using std::unordered_map;
using std::vector;
using std::set;
using std::pair;
using std::make_pair;
using std::max;
using std::queue;

/* global data */

static struct int_term_t *terms;
static int nterms;

static struct int_link_t *links;
static int nlinks;

char *gene_buffer;

int input_size=0;
float_type pval_threshold;
unsigned int gene_count_threshold;

StrComp comp;

static unordered_map<const char*, struct int_term_t*, StrHash, StrEq> ids;
//using vector_map = pair<const char*, struct int_term_t*>;
//static vector<vector_map> ids;
//static unordered_set<const char*, StrHash, StrEq> gene_names;
static vector<const char*> gene_names;
static vector<const char*> input_genes;
static vector<struct int_term_t*> roots;

/* statistics functions */
float_type *logs;
int nlogs;

static inline float_type hypergeometric(int k, int n, int K, int N)
{
	return exp((logs[K] + logs[N-K] + logs[n] + logs[N-n]) -
		(logs[k] + logs[K-k] + logs[n-k] + logs[N-K-n+k] + logs[N]));
}

static inline float_type score(int matched_genes, int assoc_genes, int input_genes_num, int all_genes)
{
	return ((float_type)matched_genes/(float_type)input_genes_num) / ((float_type)assoc_genes/(float_type)all_genes);
}

/* output thingy */
char *outbuf;
int outcursor,outlen;

static inline void print_char(char c)
{
	outbuf[outcursor++]=c;
}

#define print_buf(BUF,LEN) do {\
		memcpy(outbuf+outcursor,(BUF),(LEN));\
		outcursor+=(LEN);\
	} while (0)

static void print_flush()
{
	if (outlen-outcursor<4096)
	{
		fwrite(outbuf,outcursor,1,stdout);
		outcursor=0;
	}
}

/* file opening */

struct mapped_file_t
{
	int fd;
	size_t size;
	const char *data;
} obo,gaf;

static int open_file(const char *name, struct mapped_file_t *file)
{
	int fd;
	struct stat st;

	fd=open(name,O_RDONLY|O_NOCTTY);
	if (fd<0)
	{
		syslog(LOG_ERR,"can't open file '%s' (%d - %s)",name,errno,strerror(errno));
		return 0;
	}

	if (fstat(fd,&st)<0)
	{
		syslog(LOG_ERR,"can't stat file '%s' (%d - %s)",name,errno,strerror(errno));
		return 0;
	}
	file->size=st.st_size;

	file->data=(char*)mmap(NULL,st.st_size,PROT_READ|PROT_WRITE,MAP_FILE|MAP_PRIVATE,fd,0);
	if (!file->data)
	{
		syslog(LOG_ERR,"can't mmap file '%s' (%d - %s)",name,errno,strerror(errno));
		return 0;
	}
	madvise((void*)file->data,st.st_size,MADV_SEQUENTIAL);
	return 1;
}

/* file initialization */

static void load_obo()
{
	terms=(struct int_term_t*)(obo.data+((int*)obo.data)[0]);
	nterms=((int*)obo.data)[1];

	ids.reserve(nterms);

	#pragma omp for
	for (int i=0;i<nterms;i++)
	{
		terms[i].id.str=obo.data+terms[i].id.idx;
		terms[i].name.str=obo.data+terms[i].name.idx;
		terms[i].prerender.str=obo.data+terms[i].prerender.idx;
		new (&terms[i].genes) decltype(terms[i].genes)(1024);
		new (&terms[i].intersect) decltype(terms[i].intersect)(16);
		terms[i].parents.ptr=(struct int_term_t**)(obo.data+terms[i].parents.idx);
		terms[i].children.ptr=(struct int_term_t**)(obo.data+terms[i].children.idx);
		for (int j=0;j<terms[i].nparents;j++)
		{
			terms[i].parents.ptr[j]=(struct int_term_t*)(obo.data+(ptrdiff_t)terms[i].parents.ptr[j]);
		}
		for (int j=0;j<terms[i].nchildren;j++)
		{
			terms[i].children.ptr[j]=(struct int_term_t*)(obo.data+(ptrdiff_t)terms[i].children.ptr[j]);
		}
	}
	
	for (int i=0;i<nterms;i++)
	{
		ids[terms[i].id.str]=terms+i;
		//ids.push_back(make_pair(terms[i].id.str,terms+i));
		if (terms[i].nparents==0)
		{
			roots.push_back(terms+i);
		}
		if (!strcmp(terms[i].id.str,"GO:0003676"))
		{
			for (int j=0;j<terms[i].nchildren;j++) { fprintf(stderr,"%s: child [%s]\n",terms[i].id.str,terms[i].children.ptr[j]->id.str); }
		}
		if (!strcmp(terms[i].id.str,"GO:0001067"))
		{
			for (int j=0;j<terms[i].nparents;j++) { fprintf(stderr,"%s: parent [%s]\n",terms[i].id.str,terms[i].parents.ptr[j]->id.str); }
		}
	}
}

int clashcounter,allcounter;
static void recurse_term(const char *gene, struct int_term_t *root)
{
	bool ret=root->genes.push_gene(gene);
	clashcounter+=(ret==false);
	allcounter++;
	if (!ret)
	{
		//fprintf(stderr,"gene %s is already on term %s\n",gene,root->id.str);
		return;
	}
	root->genes_len+=strlen_bin(gene)+3;
	for (int i=0;i<root->nparents;i++)
	{
		recurse_term(gene,root->parents.ptr[i]);
	}
}

static void load_gaf()
{
	const int ngenes=((int*)gaf.data)[1];
	links=(struct int_link_t*)(gaf.data+((int*)gaf.data)[0]);
	nlinks=((int*)gaf.data)[2];
	logs=(float_type*)(gaf.data+((int*)gaf.data)[0]+(nlinks+ngenes)*sizeof(struct int_link_t));

	gene_names.reserve(ngenes);

	for (int i=0;i<nlinks+ngenes;i++)
	{
		links[i].gene.str=gaf.data+links[i].gene.idx;
		links[i].term.str=gaf.data+links[i].term.idx;
	}

	for (int i=0;i<ngenes;i++)
	{
		gene_names.push_back(links[i].gene.str);
	}
	links+=ngenes;
	for (int i=0;i<nlinks;i++)
	{
		auto it=ids.find(links[i].term.str);
		//vector_map val=make_pair(links[i].term.str,nullptr);
		//auto it=equal_range(ids.begin(),ids.end(),val,[] (const vector_map &g1, const vector_map &g2) -> bool {
		//	return comp(g1.first,g2.first);
		//});
		if (it==ids.end())
		//if (it.first==it.second)
		{
			//syslog(LOG_WARNING,"tried to attach gene [%s] to unknown term [%s]",links[i].gene.str,links[i].term.str);
			continue;
		}
		recurse_term(links[i].gene.str,(*it).second);
	}
}

static void load_input(FILE *inp)
{
	char buf[4096];
	size_t *hash=(size_t*)buf;
	char *gene=buf+sizeof(size_t);
	gene[4095]=0;
	while (fscanf(inp,"%4095s",gene)==1)
	{
		//*hash=StrHash::realhash(gene);
		//auto i=gene_names.find(gene);
		input_size+=strlen(gene)+20;
		auto i=equal_range(gene_names.begin(),gene_names.end(),gene,comp);
		//if (i==gene_names.end())
		if (i.first==i.second)
		{
			syslog(LOG_WARNING,"unknown gene [%s] in association",gene);
			continue;
		}
		input_genes.push_back(*i.first);
	}
	sort(input_genes.begin(),input_genes.end(),[] (const char *g1, const char *g2) -> bool { return g1<g2; });
}

/* main processing implementations */

#ifdef DEBUG
void gene_finder(struct int_term_t *root)
{
	bool hasit=false;
	for (auto &i : root->genes)
	{
		if (!strcmp(i,"DDB_G0290751"))
		{
			hasit=true;
			break;
		}
	}
	if (!hasit) { return; }
	printf(">%s\n",root->id.str);
	for (int i=0;i<root->nchildren;i++)
	{
		printf("<%s - %s\n",root->id.str,root->children.ptr[i]->id.str);
		gene_finder(root->children.ptr[i]);
	}
}
#endif

/* try to match the gene to the term and recursively propagate to children if matched */
static void propagate_intersections(const char *gene, struct int_term_t *root)
{
	if (!root->genes.has_gene(gene)) { return; }
	if (root->intersect.push_gene(gene))
	{
		root->intersect_len+=strlen_bin(gene)+3;
	}
	for (int i=0;i<root->nchildren;i++)
	{
		propagate_intersections(gene,root->children.ptr[i]);
	}
}

/* is the germ eligible for rendering to output? */
static inline bool term_eligible(struct int_term_t *root)
{
	return (root->pval<=pval_threshold && root->intersect.size()>=gene_count_threshold);
}

/*
 * old-style tree dumper; structure the output dict into the full tree
 */

/* traverse the immediate children tree to find all eligible children */
static void termtree_find_eligible_children(struct int_term_t *root, queue<struct int_term_t*> &q, decltype(root->color) color)
{
	for (int i=root->nchildren-1;i>=0;i--)
	{
		struct int_term_t *t=root->children.ptr[i];
		if (t->color==color) { continue; }
		t->color=color;
		if (term_eligible(t)) { q.push(t); } else { termtree_find_eligible_children(t,q,color); }
	}
}

static void termtree_dump_term(struct int_term_t *root);
static void termtree_dump_children(struct int_term_t *root)
{
	queue<struct int_term_t*> children;

	termtree_find_eligible_children(root,children,(decltype(root->color))root);
	bool first=true;
	while (!children.empty())
	{
		struct int_term_t *t=children.front(); children.pop();
		if (!first) { print_char(','); }
		first=false;
		termtree_dump_term(t);
	}
}

static void termtree_dump_term(struct int_term_t *root)
{
	root->dumped=1;
	print_buf(root->prerender.str,root->prerender_len);
	print_buf(root->genebuf,root->intersect_len);
	print_buf(",\"children\":[",13);
	termtree_dump_children(root);
	print_char(']');
	print_char('}');
	print_flush();
}

/* prepare as much as possible for term dumping, in parallel */
static void term_preprocess(struct int_term_t *root)
{
	if (!term_eligible(root)) { return; }

	size_to_string((char*)root->prerender.str+16,root->intersect.size());
	double_to_string((char*)root->prerender.str+32,root->pval);
	double_to_string((char*)root->prerender.str+72,root->score);
	size_to_string((char*)root->prerender.str+112,root->genes.size());

	/*stringify matched gene list*/
	char *tbuf=root->genebuf;
	for (const char *g : root->intersect)
	{
		*(tbuf++)='"';
		memcpy(tbuf,g,strlen_bin(g));
		tbuf+=strlen_bin(g);
		*(tbuf++)='"';
		*(tbuf++)=',';
	}
	if (root->intersect.size()) { tbuf[-1]=' '; }
	*(tbuf++)=']';

	/*stringify whole associated gene list*/
	*(tbuf++)='"';
	memcpy(tbuf,root->id.str,strlen_bin(root->id.str));
	tbuf+=strlen_bin(root->id.str);
	*(tbuf++)='"';
	*(tbuf++)=':';
	*(tbuf++)='[';
	for (const char *g : root->genes)
	{
		*(tbuf++)='"';
		memcpy(tbuf,g,strlen_bin(g));
		tbuf+=strlen_bin(g);
		*(tbuf++)='"';
		*(tbuf++)=',';
	}
	if (root->genes.size()) { tbuf[-1]=' '; }
	*(tbuf++)=']';
}

/* main driver */

int main(int argc, const char *argv[])
{
	if (argc<6)
	{
		fprintf(stderr,"invalid usage, need 5 parameters\n");
		return 1; /* usage error */
	}

	openlog("gene_tree_processor",LOG_CONS|LOG_NDELAY|LOG_PID|LOG_PERROR,LOG_USER);

	{
		double tmp;
		sscanf(argv[1],"%lf",&tmp);
		pval_threshold=tmp;
		sscanf(argv[2],"%u",&gene_count_threshold);
	}

	if (!open_file(argv[3],&obo)) { return 2; /* system error */ }
	if (!open_file(argv[4],&gaf)) { return 2; /* system error */ }

	load_obo();
	load_gaf();

	FILE *inp=fopen(argv[5],"rb");
	load_input(inp);

	fprintf(stderr,"loaded %d terms (%ld roots) and processed %d gene associations\n",nterms,roots.size(),nlinks);
	fprintf(stderr,"gene linkage finished with %d clashes out of %d total inserts\n",clashcounter,allcounter);
	fprintf(stderr,"loaded %ld input gene names, total %ld genes\n",input_genes.size(),gene_names.size());

	/* match input genes to the term tree */
	const int numinps=input_genes.size();
	const int numroots=roots.size();
	//#pragma omp for collapse(2)
	for (int i=0;i<numinps;i++)
	{
		for (int j=0;j<numroots;j++)
		{
			propagate_intersections(input_genes[i],roots[j]);
		}
	}

	/* prepare buffer for per-term gene lists */
	int gbsize=0;
	for (int i=0;i<nterms;i++)
	{
		terms[i].genes_len+=1+strlen_bin(terms[i].id.str)+4; // closing ] and term id
		terms[i].intersect_len++; // closing ]
		gbsize += terms[i].genes_len
			+ terms[i].intersect_len;
	}
	gene_buffer=(char*)malloc(gbsize);
	terms[0].genebuf=gene_buffer;
	for (int i=1;i<nterms;i++)
	{
		terms[i].genebuf=terms[i-1].genebuf
			+ terms[i-1].intersect_len
			+ terms[i-1].genes_len;
	}

	/* calculate pval (hypergeometric probability) for every term
	 * and pre-render part of the output stringification
	 */
	int ngenes=gene_names.size();
	#pragma omp for
	for (int i=0;i<nterms;i++)
	{
		terms[i].pval=hypergeometric(
			terms[i].intersect.size(),
			numinps,
			terms[i].genes.size(),
			ngenes);
		if (terms[i].genes.size()>1)
		{
			terms[i].score=score(
				terms[i].intersect.size(),
				terms[i].genes.size(),
				numinps,
				ngenes);
		}
		term_preprocess(terms+i);
	}

	/* generate the output monstrosity */
	for (int i=0;i<nterms;i++)
	{
		outlen+=terms[i].intersect_len+terms[i].genes_len+terms[i].prerender_len;
	}
	outbuf=(char*)malloc(outlen*10);
	if (!outbuf)
	{
		syslog(LOG_ERR,"can't allocate output buffer (%d - %s) (estimate %d)",errno,strerror(errno),outlen);
		return 2;
	}

	print_buf("{\"total_genes\":",15);
	outcursor+=size_to_string(outbuf+outcursor,gene_names.size());

	print_buf(",\"tree\":{",9);
	bool prevroot=false;
	for (int i=(int)roots.size()-1;i>=0;i--)
	{
		if (!strcmp(roots[i]->id.str,"GO:0008150"))
		{
			if (prevroot) { outbuf[outcursor++]=','; }
			prevroot=true;
			print_buf("\"BP\":[",6);
		}
		else if (!strcmp(roots[i]->id.str,"GO:0005575"))
		{
			if (prevroot) { outbuf[outcursor++]=','; }
			prevroot=true;
			print_buf("\"CC\":[",6);
		}
		else if (!strcmp(roots[i]->id.str,"GO:0003674"))
		{
			if (prevroot) { outbuf[outcursor++]=','; }
			prevroot=true;
			print_buf("\"MF\":[",6);
		}
		else
		{
			continue;
		}
		//dump_term(roots[i]);
		//if (i) { print_char(','); }
		bool prev=false;
		for (int j=0;j<roots[i]->nchildren;j++)
		{
			if (!term_eligible(roots[i]->children.ptr[j])) { continue; }
			if (prev) { outbuf[outcursor++]=','; }
			prev=true;
			termtree_dump_term(roots[i]->children.ptr[j]);
		}
		outbuf[outcursor++]=']';
	}
	print_char('}');

	print_buf(",\"gene_associations\":{",22);
	bool first=true;
	for (int i=0;i<nterms;i++)
	{
		if (!terms[i].dumped) { continue; }
		if (!first) { print_char(','); }
		first=false;
		print_buf(terms[i].genebuf+terms[i].intersect_len,terms[i].genes_len);
		print_flush();
	}

	print_char('}');
	print_char('}');
#ifdef DEBUG
	print_char('\n');
	print_char('\n');
	for (int i=0;i<nterms;i++)
	{
		if (
			//strcmp(terms[i].id.str,"GO:0065007") &&
			//strcmp(terms[i].id.str,"GO:0008194") &&
			//strcmp(terms[i].id.str,"GO:0097359") &&
			strcmp(terms[i].id.str,"GO:0009225") &&
			//strcmp(terms[i].id.str,"GO:0009100") &&
			true
		) { continue; }
		dump_term(&terms[i]);
		print_char('\n');
		print_char('\n');
	}
	const char *gene="DDB_G0283005";
	int matches=0;
	for (int i=0;i<nterms;i++)
	{
		for (auto &j : terms[i].genes)
		{
			if (!strcmp(j,gene))
			{
				printf(" -> %s | %s\n",terms[i].id.str,terms[i].name.str);
				matches++;
				break;
			}
		}
	}
	printf("matches: %d\n",matches);
	printf("gene DDB_G0290751:\n");
	for (auto &i : roots)
	{
		gene_finder(i);
	}
#endif

	outlen=0;
	print_flush();

	/*for (int i=0;i<nterms;i++)
	{
		printf(" %ld ",terms[i].intersect.size());
	}
	printf("\n");*/

	return 0;
}

