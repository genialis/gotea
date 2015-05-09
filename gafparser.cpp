#include<stdlib.h> /* for off_t and NULL*/

#include"preprocessor.h"

static inline int iswhite(char c) { return c==' ' || c=='\r'; }

static char *consume_string()
{
	off_t prevlen=strings_len;

	for (; gaf.cursor<gaf.size && iswhite(gaf.data[gaf.cursor]); gaf.cursor++) ;
	if (gaf.cursor==gaf.size) { return NULL; }

	for (; gaf.cursor<gaf.size && gaf.data[gaf.cursor]!='\n' && gaf.data[gaf.cursor]!='\t'; gaf.cursor++)
	{
		strings[strings_len++]=gaf.data[gaf.cursor];
	}
	if (gaf.cursor!=gaf.size && gaf.data[gaf.cursor]!='\n') { gaf.cursor++; }

	for (; iswhite(strings[strings_len-1]); strings_len--) ;
	if (strings_len==prevlen) { return NULL; }
	strings[strings_len++]=0;
	return strings+prevlen;
}

int parse_line_gaf(char **fields, int max_fields)
{
	int i;
	char *str;

	for (; gaf.cursor<gaf.size && iswhite(gaf.data[gaf.cursor]); gaf.cursor++) ;
	if (gaf.cursor==gaf.size) { return 1; } /* input exhausted */

	switch (gaf.data[gaf.cursor])
	{
	case '!':
	case '\n':
		for (gaf.cursor++; gaf.cursor<gaf.size && gaf.data[gaf.cursor]!='\n'; gaf.cursor++) ;
		gaf.cursor++;
		return 2; /* empty */
	default:
		for (i=0;;i++)
		{
			str=consume_string();
			if (i==max_fields && !(gaf.data[gaf.cursor]=='\n' || gaf.cursor>=gaf.size)) { return 3; /* malformed input file; too many fields */ }
			if (gaf.data[gaf.cursor]=='\n' || gaf.cursor>=gaf.size) { break; }
			if (i==max_fields) { break; }
			fields[i]=str;
		}
		for (;i<max_fields;i++) { fields[i]=NULL; }
		gaf.cursor++;
		return 0;
	}
	return 0;
}

