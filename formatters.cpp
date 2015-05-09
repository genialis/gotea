#ifdef TEST
#  include<stdio.h>
#  include<stdlib.h>
#  include<string.h>
#endif

/* handles only unsigned integers up to 1e9 */
static int size_to_string(char *where, size_t what)
{
	static const size_t lims[]={
		1,
		10,
		100,
		1000,
		10000,
		100000,
		1000000,
		10000000,
		100000000,
		1000000000
	};
	int lim,i;

	for (lim=0;lims[lim+1]<=what;lim++) ;
	for (i=0;lim>=0;i++)
	{
		where[i]=what/lims[lim]+'0';
		what%=lims[lim];
		lim--;
	}
	return i;
}

#include"dtoa_milo.h"
static inline int double_to_string(char *where, double what)
{
	return dtoa_milo(what,where);
}

#ifdef TEST
int main()
{
	char buf1[20],buf2[20];
	for (int i=0;i<=1000000000;i++)
	{
		sprintf(buf1,"%d",i);
		buf2[ size_to_string(buf2,(size_t)i) ]=0;
		if (strcmp(buf1,buf2)) { printf("%s -- %s\n",buf1,buf2); }
	}
}
#endif
