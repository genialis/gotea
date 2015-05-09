#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<sys/mman.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>
#include<fcntl.h>
#include<syslog.h>
#include<errno.h>

#include<map>
using std::map;

#include"preprocessor.h"

/* file opening */

struct mapped_file_t obo,gaf;
char *strings;
off_t strings_len;

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

	file->data=(char*)mmap(NULL,st.st_size,PROT_READ,MAP_FILE|MAP_PRIVATE,fd,0);
	if (!file->data)
	{
		syslog(LOG_ERR,"can't mmap file '%s' (%d - %s)",name,errno,strerror(errno));
		return 0;
	}
	madvise((void*)file->data,st.st_size,MADV_SEQUENTIAL);
	return 1;
}

/* obo and gaf processor implementations */

static int obo_processor(char **out_buf, size_t *buf_len)
{
	int ret;
	struct term_t *term=NULL;
	char *sect=NULL,*key,*val,*psect;
	for (;;)
	{
		psect=sect;
		ret=parse_line_obo(&sect,&key,&val);
		if (ret==1) { break; }
		if (ret==2) { continue; }
		if (ret==3)
		{
			syslog(LOG_ERR,"malformed input, sorry (cursor==%ld)",obo.cursor);
			return 3;
		}
		if (ret==0)
		{
			if (sect==NULL) { continue; }
			if (strcmp(sect,"Term") && strcmp(psect,"Term")) { continue; }
			if (key==NULL)
			{
				add_term(term);
				term=new_term();
			}
			else
			{
				if (!strcmp(key,"id"))
				{
					term->id=val;
				}
				if (!strcmp(key,"name"))
				{
					term->name=val;
				}
				if (!strcmp(key,"is_a"))
				{
					term->isa.insert(val);
				}
				if (!strcmp(key,"relationship"))
				{
					if (!strncmp(val,"part_of ",8))
					{
						term->isa.insert(val+8);
					}
					if (!strncmp(val,"has_part ",9))
					{
						term->isa.insert(val+9);
					}
					if (!strncmp(val,"regulates ",10))
					{
						term->isa.insert(val+10);
					}
					if (!strncmp(val,"positively_regulates ",21))
					{
						term->isa.insert(val+21);
					}
					if (!strncmp(val,"negatively_regulates ",21))
					{
						term->isa.insert(val+21);
					}
				}
				if (!strcmp(key,"is_obsolete") && !strcmp(val,"true"))
				{
					//term->obsolete=true; //XXX TODO odkomentiraj
				}
			}
		}
	}
	if (!strcmp(sect,"Term")) { add_term(term); }
	link_terms();

	return coagulate_terms(out_buf,buf_len)?0:4;
}

static int gaf_processor(char **out_buf, size_t *buf_len)
{
	char *fields[17];
	int ret;
	for (;;)
	{
		ret=parse_line_gaf(fields,17);
		if (ret==1) { break; }
		if (ret==2) { continue; }
		if (ret==3)
		{
			syslog(LOG_ERR,"malformed input, sorry (cursor==%ld)",gaf.cursor);
			return 3;
		}
		if (ret==0)
		{
			if (!fields[3] || (!strstr(fields[3],"NOT") && !strstr(fields[3],"colocalizes_with") && !strstr(fields[3],"contributes_to")))
			{
				add_link(fields[1],fields[4]);
			}
		}
	}
	return coagulate_links(out_buf,buf_len)?0:4;
}

/* driver function */

int main(int argc, char *argv[])
{
	int ret;
	char *out_buf;
	size_t out_buf_len;
	char prog_name[20];
	struct mapped_file_t *f=NULL;
	int (*proc)(char**, size_t*) = NULL;

	openlog("obo_gaf_processor",LOG_CONS|LOG_NDELAY|LOG_PID|LOG_PERROR,LOG_USER);

	if (argc!=4)
	{
		syslog(LOG_ERR,"wrong number of arguments");
		return 1; /* usage error */
	}
	if (strcmp(argv[1],"gaf") && strcmp(argv[1],"obo"))
	{
		syslog(LOG_ERR,"processor type should be either gaf or obo");
		return 1; /* usage error */
	}

	closelog();
	sprintf(prog_name,"%s_processor",argv[1]);
	openlog(prog_name,LOG_CONS|LOG_NDELAY|LOG_PID|LOG_PERROR,LOG_USER);

	if (!strcmp(argv[1],"obo"))
	{
		f=&obo;
		proc=obo_processor;
	}
	else if (!strcmp(argv[1],"gaf"))
	{
		f=&gaf;
		proc=gaf_processor;
	}

	if (!open_file(argv[2],f)) { return 2; /* system error */ }
	strings=(char*)malloc(f->size);
	if (!strings)
	{
		syslog(LOG_ERR,"can't allocate string parse buffer for input file (%d - %s)",errno,strerror(errno));
		return 2;
	}

	ret=proc(&out_buf,&out_buf_len);
	if (ret==0)
	{
		int fd=open(argv[3],O_WRONLY|O_CREAT|O_NOCTTY|O_TRUNC,S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH);
		write(fd,out_buf,out_buf_len);
		close(fd);
	}

	return ret;
}

