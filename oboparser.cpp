#include<stdlib.h> /*for off_t and NULL*/

#include"preprocessor.h"

static inline int iswhite(char c) { return c==' ' || c=='\t' || c=='\r'; }

static inline void consume_line()
{
	obo.cursor++;
	if (obo.data[obo.cursor-1]=='\n') { return; }
	for (; obo.cursor<obo.size && (obo.data[obo.cursor]!='\n' || obo.data[obo.cursor-1]=='\\'); obo.cursor++) ;
	if (obo.cursor<obo.size) { obo.cursor++; }
}

static int consume_string(char **parsed, char terminator)
{
	*parsed=strings+strings_len;
	for (; obo.cursor<obo.size && iswhite(obo.data[obo.cursor]); obo.cursor++) ;
	if (obo.cursor>=obo.size || ( obo.data[obo.cursor]=='\n' || obo.data[obo.cursor]==terminator || (obo.data[obo.cursor]=='!' && terminator=='\n') ))
	{
		return 1; /* empty */
	}
	for (;; obo.cursor++)
	{
		if (obo.cursor>=obo.size || ( (obo.data[obo.cursor]==terminator || obo.data[obo.cursor]=='\n' || (obo.data[obo.cursor]=='!' && terminator=='\n')) && obo.data[obo.cursor-1]!='\\' ))
		{
			for (; iswhite(strings[strings_len-1]); strings_len--) ;
			strings[strings_len++]=0;
			if ((obo.cursor==obo.size || obo.data[obo.cursor]!=terminator) && terminator!='\n') { return 2; } /* unterminated */
			return 0; /* OK */
		}
		if ((obo.data[obo.cursor]==terminator || obo.data[obo.cursor]=='\n' || (obo.data[obo.cursor]=='!' && terminator=='\n')) && obo.data[obo.cursor-1]=='\\')
		{
			strings[strings_len-1]=obo.data[obo.cursor];
			continue;
		}
		strings[strings_len++]=obo.data[obo.cursor];
	}
}

int parse_line_obo(char **section, char **key, char **value)
{
	static char *last_section,*last_key,*last_value;
	int ret;

	for (; obo.cursor<obo.size && iswhite(obo.data[obo.cursor]); obo.cursor++) ;
	if (obo.cursor==obo.size) { return 1; } /* input exhausted */

	switch (obo.data[obo.cursor])
	{
	case '!':
	case '\n':
		consume_line();
		return 2; /* empty */
	case '[':
		obo.cursor++;
		ret=consume_string(&last_section,']');
		if (ret) { return 3; } /* malformed input file */
		obo.cursor++;
		consume_line();
		*section=last_section;
		*key=NULL;
		*value=NULL;
		return 0; /* ok */
	default:
		ret=consume_string(&last_key,':');
		if (ret) { return 3; } /* malformed input file */
		obo.cursor++;
		consume_string(&last_value,'\n');
		consume_line();
		*section=last_section;
		*key=last_key;
		*value=last_value;
		return 0; /* ok */
	}
	return 0;
}

