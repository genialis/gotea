#ifndef __PREPROCESSOR_H__
#define __PREPROCESSOR_H__

#include<string.h>
#include<stddef.h>

#include<set>

struct ltstr
{
	bool operator()(const char *str1, const char *str2) const
	{
		return strcmp(str1,str2)<0;
	}
};

struct mapped_file_t
{
	const char *data;
	off_t size,cursor;
};

extern struct mapped_file_t obo,gaf;

extern char *strings;
extern off_t strings_len;

static inline int align(int num, int align)
{
	if (num%align) { num+=align-num%align; }
	return num;
}

/*OBO parser*/

int parse_line_obo(char **section, char **key, char **value);

struct term_t
{
	std::set<struct term_t*> parents,children;
	std::set<const char*, ltstr> isa,isa_for; /*isa_for is parent->child link cache in order to fill out isa's*/
	const char *id,*name;
	bool obsolete;

	ptrdiff_t i_parents,i_children,i_this,i_render;
	int arr_index;
};

struct term_t *new_term();
void add_term(struct term_t *term);
int link_terms();
int coagulate_terms(char **buf, size_t *buf_len);

/*GAF parser*/

int parse_line_gaf(char **fields, int max_fields);

void add_link(const char *gene_name, const char *term_name);
int coagulate_links(char **buf, size_t *buf_len);

#endif /* __PREPROCESSOR_H__ */
